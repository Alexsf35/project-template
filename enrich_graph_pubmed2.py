import os
import time
import urllib.request
import urllib.parse
import json
import pandas as pd
import functools
import xml.etree.ElementTree as ET
from neo4j import GraphDatabase

print = functools.partial(print, flush=True)

# --- CONFIGURATION ---
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "")

PUBMED_DATA_DIR = "/src/data/pubmed"
PUBS_FILE = os.path.join(PUBMED_DATA_DIR, "publications.csv")
EDGES_FILE = os.path.join(PUBMED_DATA_DIR, "entity_publication.csv")
PROCESSED_FILE = os.path.join(PUBMED_DATA_DIR, "processed_entities.csv")

class PubMedEnricher:
    def __init__(self, uri, user, password):
        self.uri = uri
        self.user = user
        self.password = password
        self.driver = None
        os.makedirs(PUBMED_DATA_DIR, exist_ok=True)

    def connect(self, retries=15, delay=5):
        for i in range(retries):
            try:
                self.driver = GraphDatabase.driver(self.uri, auth=(self.user, self.password))
                self.driver.verify_connectivity()
                print("✅ Connected to Neo4j successfully!")
                return True
            except Exception:
                time.sleep(delay)
        return False

    def get_all_biological_entities(self):
        query = """
        MATCH (n)
        WHERE n:ModelMetabolite OR n:ModelGene
        RETURN n.id AS id, split(n.name, ' [')[0] AS clean_name
        """
        entities = []
        with self.driver.session() as session:
            result = session.run(query)
            for record in result:
                entities.append({"id": record["id"], "name": record["clean_name"]})
        return entities

    def get_processed_entity_ids(self):
        if os.path.exists(PROCESSED_FILE):
            try:
                df = pd.read_csv(PROCESSED_FILE)
                return set(df['entity_id'].astype(str).tolist())
            except Exception:
                pass
        return set()

    def bulk_search_pmids(self, entities, max_results=3):
        """STEP 1: Find all article IDs as fast as possible."""
        entity_to_pmids = []
        all_pmids = set()

        print(f"🔍 Phase 1: Searching PubMed IDs for {len(entities)} entities...")
        for idx, entity in enumerate(entities):
            query_term = entity['name']
            full_query = f"{query_term}[Title/Abstract] AND Escherichia coli[Title/Abstract] AND metabolism" # dependente do modelo!!!! (Escherichia coli str. K-12 substr. MG1655)
            encoded_query = urllib.parse.quote(full_query)
            search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={encoded_query}&retmax={max_results}&retmode=json"

            try:
                with urllib.request.urlopen(search_url, timeout=10) as response:
                    data = json.loads(response.read().decode())
                    id_list = data.get("esearchresult", {}).get("idlist", [])
                    if id_list:
                        entity_to_pmids.append({"entity_id": entity['id'], "pmids": id_list})
                        all_pmids.update(id_list)
            except Exception:
                pass

            if (idx + 1) % 50 == 0:
                print(f"  ... searched {idx + 1}/{len(entities)}")
                
            time.sleep(0.35) # Max 3 requests per second to avoid NCBI ban

        return entity_to_pmids, list(all_pmids)

    def bulk_fetch_articles(self, pmids):
        """STEP 2: Download ALL abstracts in giant batches."""
        if not pmids:
            return []

        print(f"📥 Phase 2: Downloading abstracts for {len(pmids)} unique articles in bulk...")
        chunk_size = 200 # NCBI allows fetching hundreds of IDs at once
        articles = []

        for i in range(0, len(pmids), chunk_size):
            chunk = pmids[i:i + chunk_size]
            ids_string = ",".join(chunk)
            fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={ids_string}&retmode=xml"

            try:
                with urllib.request.urlopen(fetch_url, timeout=30) as response:
                    xml_data = response.read()
                    root = ET.fromstring(xml_data)

                    for article in root.findall('.//PubmedArticle'):
                        pmid = article.findtext('.//PMID')
                        if not pmid: continue

                        title = article.findtext('.//ArticleTitle') or ""

                        # THE FIX: itertext() reads through all HTML/Italic tags safely
                        abstract_parts = []
                        for a in article.findall('.//AbstractText'):
                            text = "".join(a.itertext())
                            if text: abstract_parts.append(text.strip())
                        abstract = " ".join(abstract_parts)

                        year = article.findtext('.//PubDate/Year') or article.findtext('.//ArticleDate/Year') or ""
                        authors = ", ".join([f"{a.findtext('LastName')} {a.findtext('Initials') or ''}".strip() for a in article.findall('.//Author') if a.findtext('LastName')])
                        journal = article.findtext('.//Journal/Title') or ""
                        
                        doi = ""
                        for el in article.findall('.//ArticleId'):
                            if el.get('IdType') == 'doi':
                                doi = el.text; break

                        articles.append({
                            "pmid": pmid, "title": title, "abstract": abstract,
                            "year": year, "authors": authors, "journal": journal,
                            "doi": doi, "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                        })
            except Exception as e:
                print(f"⚠️ Batch fetch failed: {e}")
            
            time.sleep(0.35)

        return articles

    def append_to_csv(self, file_path, data_list):
        if not data_list: return
        df = pd.DataFrame(data_list)
        df.to_csv(file_path, mode='a', header=not os.path.exists(file_path), index=False)

    def inject_from_csv(self):
        if not os.path.exists(PUBS_FILE) or not os.path.exists(EDGES_FILE):
            return

        print("📄 Injecting complete literature layer into Neo4j...")
        df_pubs = pd.read_csv(PUBS_FILE).drop_duplicates(subset=['pmid'])
        df_edges = pd.read_csv(EDGES_FILE)

        query = """
        UNWIND $batch AS data
        MERGE (p:Publication {id: data.pmid_curie})
        ON CREATE SET 
            p.title = data.title, p.abstract = data.abstract,
            p.year = data.year, p.authors = data.authors,
            p.journal = data.journal, p.doi = data.doi,
            p.pubmed_url = data.url, p.biolink_categories = 'Publication'
        WITH p, data
        MATCH (entity {id: data.entity_id})
        MERGE (entity)-[r:MENTIONED_IN_PUBLICATION {source: 'pubmed_cache'}]->(p)
        """
        
        batch = []
        for _, edge in df_edges.iterrows():
            pub_row = df_pubs[df_pubs['pmid'] == edge['pmid']]
            if pub_row.empty: continue
            
            pub = pub_row.iloc[0]
            try: safe_year = int(float(pub.get("year")))
            except: safe_year = None

            batch.append({
                "pmid_curie": f"PMID:{pub['pmid']}",
                "title": str(pub.get("title", "")),
                "abstract": str(pub.get("abstract", "")),
                "year": safe_year,
                "authors": str(pub.get("authors", "")),
                "journal": str(pub.get("journal", "")),
                "doi": str(pub.get("doi", "")),
                "url": str(pub.get("pubmed_url", "")),
                "entity_id": str(edge['entity_id'])
            })

        if batch:
            with self.driver.session() as session:
                session.run(query, batch=batch)
            print(f"🎉 INJECTION COMPLETE! Added {len(batch)} robust links.")

    def close(self):
        if self.driver: self.driver.close()

def run_workflow():
    print("🚀 STARTING BATCH PUBMED ENRICHMENT...")
    enricher = PubMedEnricher(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)
    if not enricher.connect(): return

    entities = enricher.get_all_biological_entities()
    processed_ids = enricher.get_processed_entity_ids()
    missing_entities = [e for e in entities if str(e['id']) not in processed_ids]
    
    if missing_entities:
        entity_pmid_map, all_pmids = enricher.bulk_search_pmids(missing_entities, max_results=3)
        articles = enricher.bulk_fetch_articles(all_pmids)
        
        print("💾 Saving clean data to local cache...")
        enricher.append_to_csv(PUBS_FILE, articles)
        
        edges_to_save = [{"entity_id": m["entity_id"], "pmid": p} for m in entity_pmid_map for p in m["pmids"]]
        enricher.append_to_csv(EDGES_FILE, edges_to_save)
        enricher.append_to_csv(PROCESSED_FILE, [{"entity_id": e['id'], "name": e['name']} for e in missing_entities])
    else:
        print("✅ Cache is up to date.")

    enricher.inject_from_csv()
    enricher.close()

if __name__ == "__main__":
    run_workflow()