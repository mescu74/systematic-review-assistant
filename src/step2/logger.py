# logger.py
import json
from datetime import datetime

def log_pubmed_search(query_string, pmid_list, log_file="search_log.json"):
    record = {
        "timestamp": datetime.utcnow().isoformat(),
        "query": query_string,
        "num_pmids": len(pmid_list),
        "pmids": pmid_list
    }
    try:
        with open(log_file, "r") as f:
            data = json.load(f)
    except FileNotFoundError:
        data = []
    
    data.append(record)
    
    with open(log_file, "w") as f:
        json.dump(data, f, indent=2)
