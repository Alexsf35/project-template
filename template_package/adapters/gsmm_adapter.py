# gsmm_adapter_minimal.py

from __future__ import annotations

import json
import re
from collections import defaultdict
from enum import Enum, auto
from typing import Optional, List, Union, Iterable, Iterator, Dict, Any

from cobra.io import read_sbml_model
from biocypher._logger import logger
from template_package.sanitize_utils import lowercase_bool_values

logger.debug(f"Loading module {__name__}.")


# ---------------- Enums ----------------

class NodeType(Enum):
    MODEL            = auto()
    ORGANISM         = auto()
    COMPARTMENT      = auto()
    MODEL_GENE       = auto()
    MODEL_REACTION   = auto()
    MODEL_METABOLITE = auto()


class EdgeType(Enum):
    MODEL_HAS_ORGANISM        = "model_has_organism"
    MODEL_HAS_GENE            = "model_has_gene"
    MODEL_HAS_REACTION        = "model_has_reaction"
    MODEL_HAS_METABOLITE      = "model_has_metabolite"
    MODEL_HAS_COMPARTMENT     = "model_has_compartment"
    GENE_REACTION             = "gene_reaction"
    REACTION_METABOLITE       = "reaction_metabolite"
    METABOLITE_IN_COMPARTMENT = "metabolite_in_compartment"
    REACTION_IN_COMPARTMENT   = "reaction_in_compartment"


# --------------- Adapter ---------------

class GSMMAdapter:
    # annotation keys we want to surface on nodes (normalized to snake_case)
    GENE_KEYS = [
        "ccds", "ecogene", "ensembl", "EnsemblGenomes-Gn", "EnsemblGenomes-Tr",
        "asap", "goa", "hgnc.symbol", "hprd", "interpro", "ncbigene", "ncbigi",
        "omim", "refseq_locus_tag", "refseq_old_locus_tag", "refseq_name",
        "refseq_synonym", "sbo", "sgd", "uniprot", "biocyc",
    ]
    RXN_KEYS = [
        "bigg.reaction", "biocyc", "ec-code", "kegg.reaction", "metanetx.reaction",
        "reactome", "reactome.reaction", "rhea", "sabiork", "sbo", "seed.reaction",
        "tcdb", "vmhreaction",
    ]
    MET_KEYS = [
        "bigg.metabolite", "biocyc", "chebi", "envipath", "hmdb", "inchi", "inchi_key",
        "kegg.compound", "kegg.drug", "kegg.glycan", "lipidmaps", "metanetx.chemical",
        "pubchem.compound", "reactome.compound", "sabiork", "sbo", "seed.compound", "slm",
        "vmhmetabolite",
    ]

    # mapping from annotation key -> CURIE prefix for preferred identifiers
    CURIE_PREFIX = {
        "biocyc": "EcoCyc",
        "bigg.reaction": "BiGG.R",
        "bigg.metabolite": "BiGG.M",
        "chebi": "CHEBI",
        "kegg.compound": "KEGG.C",
        "kegg.reaction": "KEGG.R",
        "hmdb": "HMDB",
        "rhea": "RHEA",
        "metanetx.reaction": "MNXR",
        "metanetx.chemical": "MNXC",
        "reactome.reaction": "REACTOME.R",
        "reactome.compound": "REACTOME.C",
        "uniprot": "UniProt",
        "ncbigene": "NCBIGene",
        "ecogene": "EcoGene",
        "EnsemblGenomes-Gn": "EnsemblGenomes.Gn",
        "EnsemblGenomes-Tr": "EnsemblGenomes.Tr",
        "kegg.compound": "KEGG.C",
        "bigg.metabolite": "BiGG.M",
        "bigg.reaction": "BiGG.R",
    }

    def __init__(
        self,
        sbml_paths: Union[str, Iterable[str]],
        node_types: Optional[List[NodeType]] = None,
        edge_types: Optional[List[EdgeType]] = None,
        *,
        curie_from_annotations: bool = True,
        provenance: Optional[Dict[str, str]] = None,
        organism_override: Optional[Dict[str, str]] = None,  # {'id': 'NCBITaxon:562', 'name': 'Escherichia coli K-12'}
        cache_in_memory: bool = False,
        pseudo_include_transport: bool = True,
    ):
        if isinstance(sbml_paths, str):
            sbml_paths = [sbml_paths]
        self.sbml_paths = list(sbml_paths)
        self.node_types = node_types or list(NodeType)
        self.edge_types = edge_types or list(EdgeType)
        self.curie_from_annotations = bool(curie_from_annotations)
        self.provenance = provenance or {}
        self.organism_override = organism_override or {}
        self.cache_in_memory = bool(cache_in_memory)
        self.pseudo_include_transport = bool(pseudo_include_transport)

        self.models = []
        self._pseudo_id_sets: Dict[str, Dict[str, set]] = {}
        for p in self.sbml_paths:
            logger.info(f"Loading SBML model from {p}")
            m = read_sbml_model(p)
            self.models.append(m)

            model_id = getattr(m, "id", "unknown_model")
            ex_ids = {rxn.id for rxn in getattr(m, "exchanges", []) or []}
            dm_ids = {rxn.id for rxn in getattr(m, "demands", []) or []}
            sk_ids = {rxn.id for rxn in getattr(m, "sinks", []) or []}
            self._pseudo_id_sets[model_id] = {
                "exchange": ex_ids,
                "demand": dm_ids,
                "sink": sk_ids,
            }

        # cache core iterables from cobra Model
        self._entries = [
            (getattr(m, "id", "unknown_model"), m, list(m.genes), list(m.reactions), list(m.metabolites))
            for m in self.models
        ]
        self._cached_nodes = None
        self._cached_edges = None

    # ---------- helpers ----------

    def _prov(self, d: Dict[str, Any]) -> Dict[str, Any]:
        """Merge global provenance (source/version/license) into a properties dict."""
        return {**d, **self.provenance} if self.provenance else d

    @staticmethod
    def _s(x: Optional[str]) -> str:
        """Safe string: coerce None -> '', escape quotes a little for CSV robustness."""
        return "" if x is None else str(x).replace('"', '""').replace("'", "")

    @staticmethod
    def _flat(x) -> str:
        """Flatten iterables into a '|' separated string; else str(x)."""
        if x is None:
            return ""
        if isinstance(x, (list, tuple, set)):
            return "|".join(sorted(map(str, x)))
        return str(x)

    @staticmethod
    def _norm_key(k: str) -> str:
        """Normalize annotation keys to snake_case: dots/hyphens -> underscores."""
        return re.sub(r"[\.\-]+", "_", k)

    @staticmethod
    def _first(v):
        """Stable pick of first item from a collection (sorted) or pass-through."""
        if isinstance(v, (list, tuple, set)):
            return sorted(map(str, v))[0] if v else ""
        return v

    def _curie(self, base_id: str, ann: Dict[str, Any], keys: List[str]) -> str:
        """Choose a CURIE id from annotations (preferred order), else base_id."""
        if not self.curie_from_annotations:
            return base_id
        for k in keys:
            v = ann.get(k)
            if not v:
                continue
            v = self._first(v)
            pref = self.CURIE_PREFIX.get(k)
            if pref and v:
                return f"{pref}:{v}"
        return base_id

    def _rxn_parts(self, r):
        """Yield stoichiometric parts for a reaction r."""
        for met, coeff in r.metabolites.items():
            yield dict(
                met_id=met.id,
                coefficient=float(coeff),
                abs_coefficient=float(abs(coeff)),
                role="reactant" if coeff < 0 else "product",
                ann=getattr(met, "annotation", {}) or {},
            )

    def _rxn_compartments(self, r):
        """Return list of compartments for a reaction (cobra may have r.compartments)."""
        comp = getattr(r, "compartments", None)
        if comp:
            return list(comp)
        # fallback: infer from metabolites
        return sorted({getattr(m, "compartment", "") for m in r.metabolites})

    # ---------- pseudo-reaction helpers ----------

    @staticmethod
    def _infer_comp_from_met_id(met_id: str) -> str:
        """Infer compartment suffix from metabolite identifier patterns."""
        if not met_id:
            return ""
        m = re.search(r"\[([A-Za-z0-9]+)\]$", met_id)
        if m:
            return m.group(1)
        if "_" in met_id:
            tail = met_id.rsplit("_", 1)[1]
            if re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", tail):
                return tail
        return ""

    @classmethod
    def _base_met_id_from_id(cls, met_id: str, comp: str) -> str:
        """Strip compartment suffix from metabolite id, with inference fallback."""
        if not met_id:
            return ""
        if comp and met_id.endswith(f"[{comp}]"):
            return met_id[: -(len(comp) + 2)]
        if comp and met_id.endswith(f"_{comp}"):
            return met_id[: -(len(comp) + 1)]

        inferred = cls._infer_comp_from_met_id(met_id)
        if inferred and met_id.endswith(f"_{inferred}"):
            return met_id[: -(len(inferred) + 1)]
        if inferred and met_id.endswith(f"[{inferred}]"):
            return met_id[: -(len(inferred) + 2)]
        return met_id

    def _is_transport_reaction(self, r) -> bool:
        """Detect transport by reused metabolites across multiple compartments."""
        comps = self._rxn_compartments(r)
        if len(comps) <= 1:
            return False

        met_comp_map: Dict[str, set] = {}
        for met in r.metabolites:
            met_id = getattr(met, "id", "") or ""
            comp = getattr(met, "compartment", "") or ""
            if not comp:
                comp = self._infer_comp_from_met_id(met_id)

            base = self._base_met_id_from_id(met_id, comp)
            if not base:
                continue

            met_comp_map.setdefault(base, set()).add(comp or "?")

        for _, cset in met_comp_map.items():
            known = {c for c in cset if c not in ("", "?", None)}
            if len(known) > 1:
                return True
        return False

    def _is_pseudo_reaction(self, r, pseudo_sets: Dict[str, set], rid_curie: Optional[str] = None) -> bool:
        rid_raw = getattr(r, "id", "") or ""

        # membership in exchange/demand/sink lists
        if rid_raw in pseudo_sets.get("exchange", set()):
            return True
        if rid_raw in pseudo_sets.get("demand", set()):
            return True
        if rid_raw in pseudo_sets.get("sink", set()):
            return True

        # explicit COBRA boundary flag
        if getattr(r, "boundary", False):
            return True

        # prefix-based fallback on raw id and CURIE tail
        prefixes = ("EX_", "DM_", "SK_", "BIOMASS")
        tail = ""
        if rid_curie and ":" in rid_curie:
            tail = rid_curie.split(":", 1)[1]
        elif rid_curie:
            tail = rid_curie
        if rid_raw.startswith(prefixes) or (tail and tail.startswith(prefixes)):
            return True

        # Optional transport detection
        if self.pseudo_include_transport and self._is_transport_reaction(r):
            return True

        return False

    # ---------- API ----------

    def get_nodes(self) -> Iterator:
        if self._cached_nodes is not None:
            return iter(self._cached_nodes)

        emitted_by_label = defaultdict(set)

        # model
        if NodeType.MODEL in self.node_types:
            label = "model"
            for model_id, m, *_ in self._entries:
                if model_id in emitted_by_label[label]:
                    continue
                emitted_by_label[label].add(model_id)
                yield (
                    model_id,
                    label,
                    self._prov({
                        "id": model_id,
                        "name": self._s(getattr(m, "name", "") or model_id),
                        "organism": self._s(getattr(m, "organism", "") or self.organism_override.get("name", "")),
                        "notes": json.dumps(getattr(m, "notes", {}) or {}),
                        "biolink_categories": "Model",
                    }),
                )

        # organism (one per model unless override provided)
        if NodeType.ORGANISM in self.node_types:
            label = "organism"
            for model_id, m, *_ in self._entries:
                org_name = self.organism_override.get("name") or self._s(getattr(m, "organism", "") or "")
                if not org_name:
                    continue
                org_id = self.organism_override.get("id") or org_name
                if org_id in emitted_by_label[label]:
                    continue
                emitted_by_label[label].add(org_id)
                yield (
                    org_id,
                    label,
                    self._prov({
                        "name": org_name,
                        "taxon_curie": self.organism_override.get("id", ""),
                        "strain": "",
                        "biolink_categories": "Organism",
                    }),
                )

        # compartments
        if NodeType.COMPARTMENT in self.node_types:
            label = "compartment"
            comp_seen = emitted_by_label[label]
            for _, _, _, reactions, metabolites in self._entries:
                for r in reactions:
                    for c in self._rxn_compartments(r):
                        if c and c not in comp_seen:
                            comp_seen.add(c)
                            yield (c, label, self._prov({"name": c, "biolink_categories": "CellularComponent"}))
                for m in metabolites:
                    c = getattr(m, "compartment", "")
                    if c and c not in comp_seen:
                        comp_seen.add(c)
                        yield (c, label, self._prov({"name": c, "biolink_categories": "CellularComponent"}))

        # genes
        if NodeType.MODEL_GENE in self.node_types:
            label = "model_gene"
            for _, _, genes, _, _ in self._entries:
                for g in genes:
                    ann = getattr(g, "annotation", {}) or {}
                    gid = self._curie(g.id, ann, ["uniprot", "ncbigene", "ecogene", "biocyc"])
                    if gid in emitted_by_label[label]:
                        continue
                    emitted_by_label[label].add(gid)
                    props = {
                        "name": self._s(getattr(g, "name", "") or ""),
                        "notes": json.dumps(getattr(g, "notes", {}) or {}),
                        "functional": bool(getattr(g, "functional", g.__dict__.get("_functional", False))),
                        "biolink_categories": "Gene",
                    }
                    for k in self.GENE_KEYS:
                        props[self._norm_key(k)] = self._s(self._flat(ann.get(k)))
                    lowercase_bool_values(props, in_place=True)
                    yield (gid, label, self._prov(props))

        # reactions
        if NodeType.MODEL_REACTION in self.node_types:
            label = "model_reaction"
            for model_id, _, _, reactions, _ in self._entries:
                for r in reactions:
                    ann = getattr(r, "annotation", {}) or {}
                    rid = self._curie(r.id, ann, ["biocyc", "bigg.reaction", "kegg.reaction", "rhea"])
                    if rid in emitted_by_label[label]:
                        continue
                    emitted_by_label[label].add(rid)
                    pseudo_sets = self._pseudo_id_sets.get(model_id, {})
                    props = {
                        "name": self._s(getattr(r, "name", "") or ""),
                        "notes": json.dumps(getattr(r, "notes", {}) or {}),
                        "subsystem": self._s(getattr(r, "subsystem", "") or ""),
                        "lower_bound": float(getattr(r, "lower_bound", 0.0)),
                        "upper_bound": float(getattr(r, "upper_bound", 0.0)),
                        "reversibility": bool(getattr(r, "reversibility", getattr(r, "lower_bound", 0.0) < 0)),
                        "gene_reaction_rule": self._s(getattr(r, "gene_reaction_rule", "") or ""),
                        "is_pseudo": self._is_pseudo_reaction(r, pseudo_sets, rid_curie=rid),
                        "biolink_categories": "MolecularActivity",
                    }
                    oc = getattr(r, "objective_coefficient", None)
                    if oc is not None:
                        props["objective_coefficient"] = float(oc)
                    for k in self.RXN_KEYS:
                        props[self._norm_key(k)] = self._s(self._flat(ann.get(k)))
                    lowercase_bool_values(props, in_place=True)
                    yield (rid, label, self._prov(props))

        # metabolites
        if NodeType.MODEL_METABOLITE in self.node_types:
            label = "model_metabolite"
            for _, _, _, _, mets in self._entries:
                for m in mets:
                    ann = getattr(m, "annotation", {}) or {}
                    mid = self._curie(m.id, ann, ["chebi", "bigg.metabolite", "kegg.compound", "metanetx.chemical", "biocyc"])
                    if mid in emitted_by_label[label]:
                        continue
                    emitted_by_label[label].add(mid)
                    props = {
                        "name": self._s(getattr(m, "name", "") or ""),
                        "notes": json.dumps(getattr(m, "notes", {}) or {}),
                        "formula": self._s(getattr(m, "formula", "") or ""),
                        "compartment": self._s(getattr(m, "compartment", "") or ""),
                        "charge": int(getattr(m, "charge", 0) or 0),
                        "bound": float(getattr(m, "bound", 0.0) or 0.0),
                        "is_boundary": bool(getattr(m, "boundary", False) or getattr(m, "is_boundary", False)),
                        "biolink_categories": "ChemicalEntity",
                    }
                    for k in self.MET_KEYS:
                        props[self._norm_key(k)] = self._s(self._flat(ann.get(k)))
                    lowercase_bool_values(props, in_place=True)
                    yield (mid, label, self._prov(props))

    def get_edges(self) -> Iterator:
        if self._cached_edges is not None:
            return iter(self._cached_edges)

        seen = set()

        for model_id, m, genes, reactions, metabolites in self._entries:
            # model → organism
            if EdgeType.MODEL_HAS_ORGANISM in self.edge_types:
                org_name = self.organism_override.get("name") or getattr(m, "organism", "") or ""
                if org_name:
                    org_id = self.organism_override.get("id") or org_name
                    eid = f"mdl_{model_id}_org_{org_id}"
                    if eid not in seen:
                        seen.add(eid)
                        yield (eid, model_id, org_id, "model_has_organism", self._prov({}))

            # model → gene
            if EdgeType.MODEL_HAS_GENE in self.edge_types:
                for g in genes:
                    gid = self._curie(g.id, getattr(g, "annotation", {}) or {}, ["uniprot", "ncbigene", "ecogene", "biocyc"])
                    eid = f"mdl_{model_id}_gene_{gid}"
                    if eid in seen:
                        continue
                    seen.add(eid)
                    yield (eid, model_id, gid, "model_has_gene", self._prov({}))

            # model → reaction
            if EdgeType.MODEL_HAS_REACTION in self.edge_types:
                for r in reactions:
                    rid = self._curie(r.id, getattr(r, "annotation", {}) or {}, ["biocyc", "bigg.reaction", "kegg.reaction", "rhea"])
                    eid = f"mdl_{model_id}_rxn_{rid}"
                    if eid in seen:
                        continue
                    seen.add(eid)
                    yield (eid, model_id, rid, "model_has_reaction", self._prov({}))

            # model → metabolite
            if EdgeType.MODEL_HAS_METABOLITE in self.edge_types:
                for met in metabolites:
                    mid = self._curie(met.id, getattr(met, "annotation", {}) or {}, ["chebi", "bigg.metabolite", "kegg.compound", "metanetx.chemical", "biocyc"])
                    eid = f"mdl_{model_id}_met_{mid}"
                    if eid in seen:
                        continue
                    seen.add(eid)
                    yield (eid, model_id, mid, "model_has_metabolite", self._prov({}))

            # model → compartment
            if EdgeType.MODEL_HAS_COMPARTMENT in self.edge_types:
                comp = set()
                for r in reactions:
                    comp.update(self._rxn_compartments(r))
                for met in metabolites:
                    c = getattr(met, "compartment", "")
                    if c:
                        comp.add(c)
                for c in comp:
                    if not c:
                        continue
                    eid = f"mdl_{model_id}_comp_{c}"
                    if eid in seen:
                        continue
                    seen.add(eid)
                    yield (eid, model_id, c, "model_has_compartment", self._prov({}))

            # gene → reaction (GPR)
            if EdgeType.GENE_REACTION in self.edge_types:
                for r in reactions:
                    rid = self._curie(r.id, getattr(r, "annotation", {}) or {}, ["biocyc", "bigg.reaction", "kegg.reaction", "rhea"])
                    for g in r.genes:
                        gid = self._curie(g.id, getattr(g, "annotation", {}) or {}, ["uniprot", "ncbigene", "ecogene", "biocyc"])
                        eid = f"gene_{gid}_to_rxn_{rid}"
                        if eid in seen:
                            continue
                        seen.add(eid)
                        yield (eid, gid, rid, "gene_reaction", self._prov({}))

            # reaction → metabolite (stoichiometry)
            if EdgeType.REACTION_METABOLITE in self.edge_types:
                for r in reactions:
                    rid = self._curie(r.id, getattr(r, "annotation", {}) or {}, ["biocyc", "bigg.reaction", "kegg.reaction", "rhea"])
                    for part in self._rxn_parts(r):
                        mid = self._curie(
                            part["met_id"], part["ann"], ["chebi", "bigg.metabolite", "kegg.compound", "metanetx.chemical", "biocyc"]
                        )
                        eid = f"rxn_{rid}_to_met_{mid}"
                        if eid in seen:
                            continue
                        seen.add(eid)
                        yield (
                            eid,
                            rid,
                            mid,
                            "reaction_metabolite",
                            self._prov({
                                "coefficient": part["coefficient"],
                                "abs_coefficient": part["abs_coefficient"],
                                "role": part["role"],
                            }),
                        )

            # metabolite → compartment
            if EdgeType.METABOLITE_IN_COMPARTMENT in self.edge_types:
                for met in metabolites:
                    c = getattr(met, "compartment", "")
                    if not c:
                        continue
                    mid = self._curie(met.id, getattr(met, "annotation", {}) or {}, ["chebi", "bigg.metabolite", "kegg.compound", "metanetx.chemical", "biocyc"])
                    eid = f"met_{mid}_in_comp_{c}"
                    if eid in seen:
                        continue
                    seen.add(eid)
                    yield (eid, mid, c, "metabolite_in_compartment", self._prov({}))

            # reaction → compartment
            if EdgeType.REACTION_IN_COMPARTMENT in self.edge_types:
                for r in reactions:
                    rid = self._curie(r.id, getattr(r, "annotation", {}) or {}, ["biocyc", "bigg.reaction", "kegg.reaction", "rhea"])
                    for c in self._rxn_compartments(r):
                        if not c:
                            continue
                        eid = f"rxn_{rid}_in_comp_{c}"
                        if eid in seen:
                            continue
                        seen.add(eid)
                        yield (eid, rid, c, "reaction_in_compartment", self._prov({}))

    # convenience
    def get_node_count(self) -> int:
        return sum(1 for _ in self.get_nodes())

    def get_edge_count(self) -> int:
        return sum(1 for _ in self.get_edges())
