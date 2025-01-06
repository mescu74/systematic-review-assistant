# pubmed_integration.py

from Bio import Entrez
from dotenv import load_dotenv
import os
from typing import List, Dict, Union

# Load environment variables
load_dotenv()

# Get credentials from .env file
Entrez.email = os.getenv('NCBI_EMAIL')
Entrez.api_key = os.getenv('NCBI_API_KEY')

# Add error checking
if not Entrez.email or not Entrez.api_key:
    raise ValueError("Missing NCBI_EMAIL or NCBI_API_KEY in .env file")

def pubmed_search(query: str, max_results: int = 1000) -> List[str]:
    """
    Searches PubMed with the given query and returns a list of PMIDs.
    
    Args:
        query (str): Search query string
        max_results (int): Maximum number of results to return
    
    Returns:
        List[str]: List of PubMed IDs
    
    Raises:
        Exception: If the search fails
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        pmid_list = record.get("IdList", [])
        
        if not pmid_list:
            print(f"Warning: No results found for query: {query}")
            
        return pmid_list
    
    except Exception as e:
        raise Exception(f"PubMed search failed: {str(e)}")

def pubmed_fetch_details(pmids: List[str]) -> Dict:
    """
    Fetches article details for each PMID using EFetch.
    
    Args:
        pmids (List[str]): List of PubMed IDs
    
    Returns:
        Dict: Dictionary containing article details
    
    Raises:
        ValueError: If pmids is empty
        Exception: If fetch fails
    """
    if not pmids:
        raise ValueError("No PMIDs provided")
    
    try:
        batch_ids = ",".join(pmids)
        handle = Entrez.efetch(db="pubmed", id=batch_ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        # Validate records structure
        if not records or 'PubmedArticle' not in records:
            raise ValueError("Unexpected response format from PubMed")
            
        return records
        
    except Exception as e:
        raise Exception(f"Failed to fetch PubMed details: {str(e)}")

def extract_article_info(article: Dict) -> Dict[str, str]:
    """
    Extracts key information from a PubMed article record.
    
    Args:
        article (Dict): PubMed article record
    
    Returns:
        Dict[str, str]: Dictionary containing article information
    """
    try:
        medline = article['MedlineCitation']
        article_info = medline['Article']
        
        # Handle PMID correctly (it's a string element)
        pmid = str(medline['PMID'])
        
        # Handle abstract (might be a list of StringElements)
        abstract_text = article_info.get('Abstract', {}).get('AbstractText', ['No abstract'])
        if isinstance(abstract_text, list):
            abstract = ' '.join(str(text) for text in abstract_text)
        else:
            abstract = str(abstract_text)
        
        return {
            'pmid': pmid,
            'title': str(article_info.get('ArticleTitle', 'No title')),
            'abstract': abstract,
            'journal': str(article_info.get('Journal', {}).get('Title', 'No journal')),
            'year': str(article_info.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year', 'No year'))
        }
        
    except KeyError as e:
        print(f"Warning: Could not extract some article information: {e}")
        return {
            'pmid': 'Error',
            'title': 'Error extracting information',
            'abstract': 'Error extracting information',
            'journal': 'Error',
            'year': 'Error'
        }
