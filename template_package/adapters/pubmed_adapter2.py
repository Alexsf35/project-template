import os
import uuid
import pandas as pd
from enum import Enum, auto
from typing import Optional, List, Iterator, Dict, Any
from biocypher._logger import logger

# Se tiveres o sanitize_utils no teu projeto, podes importar:
# from template_package.sanitize_utils import sanitize_row 

logger.debug(f"Loading module {__name__}.")

# -------------------------------
# Enums (Definição do Esquema)
# -------------------------------

class NodeType(Enum):
    PUBLICATION = auto()

class EdgeType(Enum):
    # Relação genérica: "Entidade (Gene/Reação) Mencionada na Publicação"
    MENTIONED_IN_PUBLICATION = "mentioned_in_publication"


# -------------------------------
# Node Class
# -------------------------------

class PubMedNode:
    """
    Representa um artigo do PubMed. 
    Estruturado para fornecer contexto rico em texto para o BioChatter (LLM).
    """
    def __init__(
        self,
        record: Dict[str, Any],
        provenance: Optional[Dict[str, Any]] = None,
    ):
        # O PMID é o identificador único. Usamos o formato CURIE padrão (PMID:XXXX)
        raw_pmid = str(record.get("pmid", "")).strip()
        self.id = f"PMID:{raw_pmid}" if raw_pmid else f"PMID:{uuid.uuid4()}"
        
        self.label = "publication"
        
        # Propriedades cruciais para a Task 3 (Explicação com LLM)
        self.properties = {
            "title": str(record.get("title", "")),
            "abstract": str(record.get("abstract", "")),
            "year": int(record.get("year", 0)) if pd.notna(record.get("year")) else None,
            "authors": str(record.get("authors", "")),
            "journal": str(record.get("journal", "")),
            "doi": str(record.get("doi", "")),                    # <--- ADICIONAR ISTO
            "pubmed_url": str(record.get("pubmed_url", "")),      # <--- ADICIONAR ISTO
            "biolink_categories": "Publication",
        }
        if provenance:
            self.properties.update(provenance)

    def get_id(self) -> str:
        return self.id

    def get_label(self) -> str:
        return self.label

    def get_properties(self) -> Dict[str, Any]:
        return self.properties


# -------------------------------
# Adapter Principal
# -------------------------------

class PubMedAdapter:
    """
    Lê dados de publicações e as suas ligações às entidades do GSMM,
    gerando nós e arestas para o BioCypher.
    """
    def __init__(
        self,
        data_dir: str,
        node_types: Optional[List[NodeType]] = None,
        edge_types: Optional[List[EdgeType]] = None,
        add_provenance: bool = True,
    ):
        self.data_dir = data_dir
        self.node_types = node_types if node_types else list(NodeType)
        self.edge_types = edge_types if edge_types else list(EdgeType)
        
        self.add_provenance = add_provenance
        self.provenance = {"source": "PubMed", "license": "NLM"} if add_provenance else {}

        # Ficheiros esperados
        self.publications_file = os.path.join(self.data_dir, "publications.csv")
        self.edges_file = os.path.join(self.data_dir, "entity_publication.csv")

    def _prov(self) -> Dict[str, Any]:
        return self.provenance

    def get_nodes(self) -> Iterator[tuple[str, str, Dict[str, Any]]]:
        """
        Lê o CSV de publicações e gera os nós.
        """
        logger.info("A gerar nós do PubMed...")
        
        if NodeType.PUBLICATION not in self.node_types:
            return

        if not os.path.exists(self.publications_file):
            logger.warning(f"Ficheiro não encontrado: {self.publications_file}")
            return

        df = pd.read_csv(self.publications_file)
        # Limpar NaN para evitar erros no Neo4j
        df = df.fillna("") 

        for rec in df.to_dict(orient="records"):
            node = PubMedNode(rec, provenance=self._prov())
            yield (node.get_id(), node.get_label(), node.get_properties())

    def get_edges(self) -> Iterator[tuple[str, str, str, str, Dict[str, Any]]]:
        """
        Lê o CSV de mapeamento e gera as ligações entre Entidades (ex: Reações) e Publicações.
        O CSV deve ter as colunas 'entity_id' (ex: BiGG.R:PFK) e 'pmid'.
        """
        logger.info("A gerar arestas de evidência bibliográfica...")
        
        if EdgeType.MENTIONED_IN_PUBLICATION not in self.edge_types:
            return

        if not os.path.exists(self.edges_file):
            logger.warning(f"Ficheiro não encontrado: {self.edges_file}")
            return

        df_edges = pd.read_csv(self.edges_file)
        
        for index, row in df_edges.iterrows():
            entity_id = str(row.get("entity_id", "")).strip()
            pmid = str(row.get("pmid", "")).strip()
            
            if not entity_id or not pmid:
                continue

            target_pmid_curie = f"PMID:{pmid}"
            edge_id = f"edge_{entity_id}_to_{pmid}"

            yield (
                edge_id,                 # ID da aresta
                entity_id,               # Nó de origem (Vem do GSMM Adapter)
                target_pmid_curie,       # Nó de destino (O artigo)
                EdgeType.MENTIONED_IN_PUBLICATION.value, # Label da aresta
                self._prov()             # Propriedades
            )

    def get_node_count(self) -> int:
        if os.path.exists(self.publications_file):
            return len(pd.read_csv(self.publications_file))
        return 0