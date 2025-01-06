import unittest
from search_utils import build_pubmed_query

class TestSearchUtils(unittest.TestCase):
    def test_basic_keyword_search(self):
        keywords = ["diabetes", "obesity"]
        result = build_pubmed_query(
            keywords=keywords,
            inclusion_criteria={},
            exclusion_criteria={}
        )
        expected = "(diabetes) OR (obesity)"
        self.assertEqual(result, expected)
    
    def test_with_inclusion_criteria(self):
        keywords = ["cancer"]
        inclusion_criteria = {
            "language": "english",
            "study_type": "clinical trial"
        }
        result = build_pubmed_query(
            keywords=keywords,
            inclusion_criteria=inclusion_criteria,
            exclusion_criteria={}
        )
        expected = "(cancer) AND english[lang] AND clinical trial[pt]"
        self.assertEqual(result, expected)
    
    def test_with_exclusion_criteria(self):
        keywords = ["alzheimer"]
        exclusion_criteria = {
            "exclude_study_type": "review"
        }
        result = build_pubmed_query(
            keywords=keywords,
            inclusion_criteria={},
            exclusion_criteria=exclusion_criteria
        )
        expected = "(alzheimer) AND NOT review[pt]"
        self.assertEqual(result, expected)
    
    def test_with_additional_filters(self):
        keywords = ["covid"]
        additional_filters = {
            "date_range": "2020:2023"
        }
        result = build_pubmed_query(
            keywords=keywords,
            inclusion_criteria={},
            exclusion_criteria={},
            additional_filters=additional_filters
        )
        expected = "(covid) AND 2020:2023[dp]"
        self.assertEqual(result, expected)
    
    def test_empty_inputs(self):
        result = build_pubmed_query(
            keywords=[],
            inclusion_criteria={},
            exclusion_criteria={}
        )
        expected = ""
        self.assertEqual(result, expected)
    
    def test_complex_query(self):
        keywords = ["depression", "anxiety"]
        inclusion_criteria = {"language": "english"}
        exclusion_criteria = {"exclude_study_type": "review"}
        additional_filters = {"date_range": "2020:2023"}
        
        result = build_pubmed_query(
            keywords=keywords,
            inclusion_criteria=inclusion_criteria,
            exclusion_criteria=exclusion_criteria,
            additional_filters=additional_filters
        )
        expected = "(depression) OR (anxiety) AND english[lang] AND NOT review[pt] AND 2020:2023[dp]"
        self.assertEqual(result, expected)

if __name__ == '__main__':
    unittest.main() 