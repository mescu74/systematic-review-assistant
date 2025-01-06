from pubmed_integration import pubmed_search, pubmed_fetch_details, extract_article_info

def test_pubmed_functions():
    print("Testing PubMed search...")
    query = "systematic review AND machine learning"
    
    try:
        pmids = pubmed_search(query, max_results=5)
        print(f"✅ Search successful! Found {len(pmids)} results")
        print(f"First few PMIDs: {pmids[:5]}\n")
        
        if pmids:
            print("Testing article details fetch...")
            records = pubmed_fetch_details(pmids[:2])
            
            for article in records['PubmedArticle']:
                info = extract_article_info(article)
                print("\nArticle Details:")
                print(f"Title: {info['title']}")
                print(f"Journal: {info['journal']} ({info['year']})")
                print(f"PMID: {info['pmid']}")
                print("Abstract:", info['abstract'][:200] + "..." if len(info['abstract']) > 200 else info['abstract'])
            
            print("\n✅ Detail fetch successful!")
            
    except Exception as e:
        print(f"❌ Error: {str(e)}")

if __name__ == "__main__":
    test_pubmed_functions() 