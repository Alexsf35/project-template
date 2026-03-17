import os
import re
import time
import pandas as pd
from neo4j import GraphDatabase

class STRINGGraphEnricher:
    def __init__(self, data_dir, taxon_id, neo4j_uri, neo4j_user, neo4j_password, threshold_pct=0.05, interaction_threshold=700):
        self.data_dir = data_dir
        self.taxon_id = str(taxon_id)
        self.threshold_pct = threshold_pct
        self.interaction_threshold = interaction_threshold
        
        self.uri = neo4j_uri
        self.user = neo4j_user
        self.password = neo4j_password
        self.driver = None
        
        # Raw Input Files
        self.raw_terms_file = os.path.join(self.data_dir, f"{self.taxon_id}.protein.enrichment.terms.v12.0.txt")
        self.raw_links_file = os.path.join(self.data_dir, f"{self.taxon_id}.protein.links.full.v12.0.txt")
        self.raw_info_file = os.path.join(self.data_dir, f"{self.taxon_id}.protein.info.v12.0.txt")
        
        # Processed Output Directory and Files
        self.processed_dir = os.path.join(self.data_dir, "processed")
        os.makedirs(self.processed_dir, exist_ok=True)
        
        self.props_file = os.path.join(self.processed_dir, f"{self.taxon_id}_protein_properties.csv")
        self.info_cache_file = os.path.join(self.processed_dir, f"{self.taxon_id}_protein_info.csv") # New cache file
        self.terms_file = os.path.join(self.processed_dir, f"{self.taxon_id}_term_nodes.csv")
        self.edges_file = os.path.join(self.processed_dir, f"{self.taxon_id}_protein_term_edges.csv")
        self.filtered_links_file = os.path.join(self.processed_dir, f"{self.taxon_id}_filtered_links.csv")

    def connect(self, retries=10, delay=5):
        print("Attempting to connect to Neo4j...")
        for i in range(retries):
            try:
                self.driver = GraphDatabase.driver(self.uri, auth=(self.user, self.password))
                self.driver.verify_connectivity()
                print("Successfully connected to Neo4j.")
                return True
            except Exception as e:
                print(f"Neo4j is starting... waiting {delay}s (Attempt {i+1}/{retries})")
                time.sleep(delay)
        print("Failed to connect to Neo4j after multiple retries.")
        return False

    def close(self):
        if self.driver:
            self.driver.close()

    def get_valid_genes(self):
        query = "MATCH (g:ModelGene) WHERE g.refseq_locus_tag IS NOT NULL RETURN g.refseq_locus_tag AS id"
        valid_genes = set()
        with self.driver.session() as session:
            result = session.run(query)
            for record in result:
                valid_genes.add(str(record["id"]).strip())
        print(f"Found {len(valid_genes)} ModelGene nodes to act as anchor points.")
        return valid_genes

    def clean_property_name(self, text):
        clean = re.sub(r'[^a-zA-Z0-9]', '_', str(text)).lower()
        return re.sub(r'_+', '_', clean).strip('_')

    def process_enrichment_terms(self, valid_genes):
        print("\n--- Phase 1: Processing Enrichment Terms and Protein Info ---")
        
        # Check cache for all generated files
        if (os.path.exists(self.props_file) and os.path.exists(self.info_cache_file) and 
            os.path.exists(self.terms_file) and os.path.exists(self.edges_file)):
            print("Cache found! All Phase 1 files already processed. Skipping Phase 1...")
            return True

        # --- Sub-Phase A: Process Enrichment Terms ---
        if not os.path.exists(self.raw_terms_file):
            print(f"Error: Could not find {self.raw_terms_file}")
            return False

        print("Loading raw STRING enrichment file...")
        df = pd.read_csv(self.raw_terms_file, sep='\t')
        
        if '#string_protein_id' in df.columns:
            df.rename(columns={'#string_protein_id': 'string_protein_id'}, inplace=True)

        df['clean_id'] = df['string_protein_id'].apply(
            lambda x: str(x).split('.')[1].strip() if '.' in str(x) else str(x).strip()
        )

        original_row_count = len(df)
        df_enrich = df[df['clean_id'].isin(valid_genes)].copy()
        
        print(f"Filtered to model-specific annotations: kept {len(df_enrich)} out of {original_row_count} rows.")

        total_unique_proteins = df_enrich['clean_id'].nunique()
        cutoff_count = total_unique_proteins * self.threshold_pct
        
        term_counts = df_enrich.groupby('term')['clean_id'].nunique().reset_index()
        term_counts.rename(columns={'clean_id': 'protein_count'}, inplace=True)

        df_with_counts = pd.merge(df_enrich, term_counts, on='term')
        
        mask_generic = df_with_counts['protein_count'] >= cutoff_count
        mask_specific = (df_with_counts['protein_count'] < cutoff_count) & (df_with_counts['protein_count'] > 1)

        df_generic = df_with_counts[mask_generic].copy()
        df_specific = df_with_counts[mask_specific].copy()

        print(f"Identified {df_generic['term'].nunique()} generic terms (Properties).")
        print(f"Identified {df_specific['term'].nunique()} specific terms (Nodes).")

        df_generic['formatted_term'] = df_generic['term'] + ": " + df_generic['description']
        df_generic['prop_key'] = df_generic['category'].apply(self.clean_property_name)
        
        grouped_props = df_generic.groupby(['clean_id', 'prop_key'])['formatted_term'].apply(
            lambda x: '|'.join(set(x))
        ).reset_index()
        
        protein_properties = grouped_props.pivot(index='clean_id', columns='prop_key', values='formatted_term').reset_index()
        protein_properties.rename(columns={'clean_id': 'gene_id'}, inplace=True)
        protein_properties.fillna("", inplace=True)
        
        protein_properties.to_csv(self.props_file, index=False)

        nodes_df = df_specific[['term', 'category', 'description']].drop_duplicates()
        nodes_df.to_csv(self.terms_file, index=False)

        edges_df = df_specific[['clean_id', 'term']].drop_duplicates()
        edges_df.rename(columns={'clean_id': 'gene_id'}, inplace=True)
        edges_df.to_csv(self.edges_file, index=False)

        # --- Sub-Phase B: Process Protein Info ---
        if not os.path.exists(self.raw_info_file):
            print(f"Error: Could not find {self.raw_info_file}")
            return False

        print("Loading raw STRING protein info file...")
        df_info = pd.read_csv(self.raw_info_file, sep='\t')
        
        # Clean ID
        if '#string_protein_id' in df_info.columns:
             df_info.rename(columns={'#string_protein_id': 'string_protein_id'}, inplace=True)

        df_info['clean_id'] = df_info['string_protein_id'].apply(
            lambda x: str(x).split('.')[1].strip() if '.' in str(x) else str(x).strip()
        )
        
        # Filter for valid proteins
        df_info_filtered = df_info[df_info['clean_id'].isin(valid_genes)].copy()
        
        # Keep only required columns for cache
        clean_info = df_info_filtered[['clean_id', 'protein_size', 'annotation']].copy()
        clean_info.rename(columns={'clean_id': 'gene_id'}, inplace=True)
        
        # Save to cache CSV
        clean_info.to_csv(self.info_cache_file, index=False)
        print("Protein info cache file generated.")

        print("Enrichment processing complete. Files generated.")
        return True

    def process_protein_links(self, valid_genes):
        print("\n--- Phase 2: Processing Protein-Protein Links ---")

        if os.path.exists(self.filtered_links_file):
            print("Cache found! Link files already processed. Skipping Phase 2...")
            return True
        
        if not os.path.exists(self.raw_links_file):
            print(f"Error: Could not find {self.raw_links_file}")
            return False
            
        print("Loading raw full STRING links file...")
        df_links = pd.read_csv(self.raw_links_file, sep=' ')
        original_count = len(df_links)
        
        df_links['p1_clean'] = df_links['protein1'].apply(lambda x: str(x).split('.')[1].strip() if '.' in str(x) else str(x).strip())
        df_links['p2_clean'] = df_links['protein2'].apply(lambda x: str(x).split('.')[1].strip() if '.' in str(x) else str(x).strip())
        
        df_filtered = df_links[(df_links['p1_clean'].isin(valid_genes)) & (df_links['p2_clean'].isin(valid_genes))].copy()
        
        df_filtered = df_filtered[df_filtered['combined_score'] >= self.interaction_threshold]
        
        print(f"Filtered interactions from {original_count} down to {len(df_filtered)} relevant edges (Score >= {self.interaction_threshold}).")
        
        # Save all relevant columns
        columns_to_keep = [
            'p1_clean', 'p2_clean', 'neighborhood', 'neighborhood_transferred', 
            'fusion', 'cooccurence', 'homology', 'coexpression', 
            'coexpression_transferred', 'experiments', 'experiments_transferred', 
            'database', 'database_transferred', 'textmining', 'textmining_transferred', 
            'combined_score'
        ]
        
        clean_links = df_filtered[columns_to_keep].fillna(0)
        clean_links.to_csv(self.filtered_links_file, index=False)
        
        print("Links processing complete. File generated.")
        return True

    def inject_to_neo4j(self):
        print("\n--- Phase 3: Neo4j Data Injection ---")

        with self.driver.session() as session:
            # 1. Create Protein Nodes
            if os.path.exists(self.props_file) and os.path.exists(self.info_cache_file):
                print("Injecting Protein nodes, merging properties from enrichment and info files...")
                df_props = pd.read_csv(self.props_file)
                df_props.fillna("", inplace=True)
                
                # Load info data
                df_info = pd.read_csv(self.info_cache_file)
                df_info.fillna("", inplace=True)
                
                # Merge the dataframes on gene_id
                df_merged_props = pd.merge(df_props, df_info, on='gene_id', how='outer')
                df_merged_props.fillna("", inplace=True)
                
                batch = []
                for _, row in df_merged_props.iterrows():
                    row_dict = row.to_dict()
                    gene_id = str(row_dict.pop("gene_id")).strip()
                    
                    dynamic_properties = {}
                    
                    # Convert protein_size back to integer if it has data
                    protein_size = row_dict.pop("protein_size")
                    annotation = row_dict.pop("annotation")
                    
                    if protein_size != "":
                         dynamic_properties["size"] = int(protein_size)
                    if annotation != "":
                         dynamic_properties["anno"] = str(annotation)
                    
                    for key, value in row_dict.items():
                        if str(value).strip() != "":
                            dynamic_properties[key] = str(value)
                            
                    batch.append({
                        "gene_id": gene_id,
                        "properties": dynamic_properties
                    })
                
                # Modified Query: ON CREATE SET and ON MATCH SET to handle size and anno
                query_proteins = """
                UNWIND $batch AS data
                MATCH (g:ModelGene) WHERE g.refseq_locus_tag = data.gene_id
                MERGE (p:Protein {id: data.gene_id})
                ON CREATE SET p += data.properties
                ON MATCH SET p += data.properties, p.protein_size = data.size, p.annotation = data.anno
                SET p.name = g.name
                SET p.uniprot_id = g.id
                MERGE (g)-[:ENCODES]->(p)
                """
                result = session.run(query_proteins, batch=batch)
                summary = result.consume()
                print(f" -> Created {summary.counters.nodes_created} Protein nodes.")
                print(f" -> Updated {summary.counters.properties_set} properties.")
                print(f" -> Created {summary.counters.relationships_created} ENCODES relationships.")
            
            # 2. Create FunctionalTerm Nodes
            if os.path.exists(self.terms_file):
                print("\nInjecting FunctionalTerm nodes...")
                df_terms = pd.read_csv(self.terms_file)
                batch = []
                for _, row in df_terms.iterrows():
                    batch.append({
                        "term": str(row["term"]),
                        "category": str(row["category"]),
                        "description": str(row["description"])
                    })
                
                query_terms = """
                UNWIND $batch AS data
                MERGE (t:FunctionalTerm {id: data.term})
                ON CREATE SET t.category = data.category, t.description = data.description
                """
                result = session.run(query_terms, batch=batch)
                summary = result.consume()
                print(f" -> Created {summary.counters.nodes_created} FunctionalTerm nodes.")

            # 3. Create Protein -> Term Edges
            if os.path.exists(self.edges_file):
                print("\nInjecting Protein-Term relationships...")
                df_edges = pd.read_csv(self.edges_file)
                batch = []
                for _, row in df_edges.iterrows():
                    batch.append({
                        "gene_id": str(row["gene_id"]).strip(),
                        "term": str(row["term"])
                    })
                
                query_term_edges = """
                UNWIND $batch AS data
                MATCH (p:Protein {id: data.gene_id})
                MATCH (t:FunctionalTerm {id: data.term})
                MERGE (p)-[:HAS_FUNCTIONAL_ANNOTATION]->(t)
                """
                result = session.run(query_term_edges, batch=batch)
                summary = result.consume()
                print(f" -> Created {summary. counters.relationships_created} HAS_FUNCTIONAL_ANNOTATION relationships.")

            # 4. Create Protein -> Protein Interaction Edges
            if os.path.exists(self.filtered_links_file):
                print("\nInjecting Protein-Protein Interaction (PPI) network...")
                df_links = pd.read_csv(self.filtered_links_file)
                batch = []
                for _, row in df_links.iterrows():
                    batch.append({
                        "p1": str(row["p1_clean"]).strip(),
                        "p2": str(row["p2_clean"]).strip(),
                        "neighborhood": int(row["neighborhood"]),
                        "neighborhood_transferred": int(row["neighborhood_transferred"]),
                        "fusion": int(row["fusion"]),
                        "cooccurence": int(row["cooccurence"]),
                        "homology": int(row["homology"]),
                        "coexpression": int(row["coexpression"]),
                        "coexpression_transferred": int(row["coexpression_transferred"]),
                        "experiments": int(row["experiments"]),
                        "experiments_transferred": int(row["experiments_transferred"]),
                        "database": int(row["database"]),
                        "database_transferred": int(row["database_transferred"]),
                        "textmining": int(row["textmining"]),
                        "textmining_transferred": int(row["textmining_transferred"]),
                        "combined_score": int(row["combined_score"])
                    })
                
                query_ppi = """
                UNWIND $batch AS data
                MATCH (p1:Protein {id: data.p1})
                MATCH (p2:Protein {id: data.p2})
                MERGE (p1)-[r:INTERACTS_WITH]->(p2)
                ON CREATE SET 
                    r.neighborhood = data.neighborhood,
                    r.neighborhood_transferred = data.neighborhood_transferred,
                    r.fusion = data.fusion,
                    r.cooccurence = data.cooccurence,
                    r.homology = data.homology,
                    r.coexpression = data.coexpression,
                    r.coexpression_transferred = data.coexpression_transferred,
                    r.experiments = data.experiments,
                    r.experiments_transferred = data.experiments_transferred,
                    r.database = data.database,
                    r.database_transferred = data.database_transferred,
                    r.textmining = data.textmining,
                    r.textmining_transferred = data.textmining_transferred,
                    r.combined_score = data.combined_score
                """
                result = session.run(query_ppi, batch=batch)
                summary = result.consume()
                print(f" -> Created {summary.counters.relationships_created} INTERACTS_WITH relationships.")

        print("\nPipeline Complete: The Graph is now fully enriched with STRING data.")

def run_pipeline():
    import os
    
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIRECTORY = os.path.join(SCRIPT_DIR, "data", "string")
    ORGANISM_ID = "511145"
    
    NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
    NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
    NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "password")
    
    enricher = STRINGGraphEnricher(
        data_dir=DATA_DIRECTORY, 
        taxon_id=ORGANISM_ID,
        neo4j_uri=NEO4J_URI,
        neo4j_user=NEO4J_USER,
        neo4j_password=NEO4J_PASSWORD,
        threshold_pct=0.05,
        interaction_threshold=700 
    )
    
    if enricher.connect():
        valid_genes = enricher.get_valid_genes()
        if valid_genes:
            enricher.process_enrichment_terms(valid_genes)
            enricher.process_protein_links(valid_genes)
            enricher.inject_to_neo4j()
        enricher.close()

if __name__ == "__main__":
    run_pipeline()