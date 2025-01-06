import unittest
from .llama_index_integration import pubmed_records_to_documents, build_llama_index_from_pubmed
from llama_index import Document, VectorStoreIndex
import os
import shutil

class TestLlamaIndexIntegration(unittest.TestCase):
    def setUp(self):
        self.test_dir = "test_pubmed_index"

    def tearDown(self):
        # Clean up the test directory if it exists
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def test_pubmed_records_to_documents(self):
        # Test data
        pubmed_records = [
            {
                "MedlineCitation": {
                    "PMID": "123456",
                    "Article": {
                        "ArticleTitle": "Sample Title",
                        "Abstract": {
                            "AbstractText": ["This is a sample abstract."]
                        }
                    }
                }
            }
        ]
        
        # Test document creation
        documents = pubmed_records_to_documents(pubmed_records)
        
        # Assertions
        self.assertEqual(len(documents), 1)
        self.assertIsInstance(documents[0], Document)
        self.assertEqual(
            documents[0].text,
            "PMID: 123456\nTitle: Sample Title\nAbstract: This is a sample abstract."
        )

    def test_build_llama_index_from_pubmed(self):
        # Test data
        pubmed_records = [
            {
                "MedlineCitation": {
                    "PMID": "123456",
                    "Article": {
                        "ArticleTitle": "Sample Title",
                        "Abstract": {
                            "AbstractText": ["This is a sample abstract."]
                        }
                    }
                }
            }
        ]
        
        # Test index creation
        index = build_llama_index_from_pubmed(pubmed_records, self.test_dir)
        
        # Assertions
        self.assertIsNotNone(index)
        self.assertIsInstance(index, VectorStoreIndex)
        self.assertTrue(os.path.exists(self.test_dir))

if __name__ == '__main__':
    unittest.main() 