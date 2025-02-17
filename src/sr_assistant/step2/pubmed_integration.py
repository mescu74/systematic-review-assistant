"""PubMed integration.

The goal is to fetch relevant study metadata from PubMed.
We use Biopython and Entrez in the protoype,
see https://github.com/biopython/biopython/blob/fc085a532f015b502d0bdcff35bd63889c00651f/Doc/Tutorial/chapter_entrez.rst

Using Entrez esearch and efetch isn't really how this is supposed to be done, we should
use NCBI WebEnv history feature. To be fixed.

Notes:
    In [69]: import sr_assistant.step2.pubmed_integration as pmi

    In [70]: pmids = pmi.pubmed_search(query="(cancer) AND (immunotherapy) AND (clinical trial)", max_results=2)

    In [71]: pmids
    Out[71]: ['39844819', '39844281']

    In [72]: details = pmi.pubmed_fetch_details(pmids=pmids)

    In [73]: details.keys()
    Out[73]: dict_keys(['PubmedBookArticle', 'PubmedArticle'])

    In [74]: len(details['PubmedBookArticle'])
    Out[74]: 0

    In [75]: len(details['PubmedArticle'])
    Out[75]: 2

    In [76]: details['PubmedArticle'][0].keys()
    Out[76]: dict_keys(['MedlineCitation', 'PubmedData'])

    In [77]: details['PubmedArticle'][0]['MedlineCitation'].keys()
    Out[77]: dict_keys(['OtherAbstract', 'SpaceFlightMission', 'OtherID', 'InvestigatorList', 'KeywordList', 'GeneralNote', 'CitationSubset', 'PMID', 'DateCompleted', 'DateRevised', 'Article', 'MedlineJournalInfo', 'ChemicalList', 'MeshHeadingList', 'CoiStatement'])

    In [78]: details['PubmedArticle'][0]['PubmedData'].keys()
    Out[78]: dict_keys(['ReferenceList', 'History', 'PublicationStatus', 'ArticleIdList'])

"""

from __future__ import annotations

import os
from typing import Any, cast

from Bio import Entrez
from dotenv import load_dotenv
from loguru import logger

# Load environment variables
load_dotenv()

# Get credentials from .env file
email = os.getenv("NCBI_EMAIL")
api_key = os.getenv("NCBI_API_KEY")

if email is None or api_key is None:
    msg = "Missing NCBI_EMAIL or NCBI_API_KEY in .env file"
    raise ValueError(msg)

# Type ignore for Entrez attributes which are dynamically set
Entrez.email = email  # type: ignore
Entrez.api_key = api_key  # type: ignore


def pubmed_search(query: str, max_results: int = 1000) -> list[str]:
    """Searches PubMed with the given query and returns a list of PMIDs."""
    try:
        # Type ignore for untyped Entrez functions
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)  # type: ignore
        record = Entrez.read(handle)  # type: ignore
        handle.close()
        pmid_list = record.get("IdList", [])  # type: ignore

        if not pmid_list:
            msg = "No results found for query: {!r}"
            logger.warning(msg, query)

        return cast(list[str], pmid_list)

    except Exception as e:
        msg = f"PubMed search failed: {e!s}"
        logger.opt(exception=True).error(msg)
        raise


def pubmed_fetch_details(pmids: list[str]) -> dict[str, Any]:
    """Fetches article details for each PMID using EFetch.

    Args:
        pmids (list[str]): List of PMIDs to fetch details for.

    Returns:
        dict[str, Any]: Dictionary mapping PMID to article details.

    Example:
        In [18]: records.keys()
        Out[18]: dict_keys(['PubmedBookArticle', 'PubmedArticle'])
        The 'PubmedArticle' key contains a list of the article details.
        'PubmedBookArticle' is an empty list typically.

    """
    if not pmids:
        msg = "No PMIDs provided"
        raise ValueError(msg)

    try:
        # Process PMIDs in batches to avoid overwhelming the server
        batch_size = 100
        all_records = {}

        for i in range(0, len(pmids), batch_size):
            batch_ids = ",".join(pmids[i : i + batch_size])

            # Type ignore for untyped Entrez functions
            handle = Entrez.efetch(db="pubmed", id=batch_ids, retmode="xml")  # type: ignore
            records = Entrez.read(handle)  # type: ignore
            handle.close()
            # In [18]: records.keys()
            # Out[18]: dict_keys(['PubmedBookArticle', 'PubmedArticle'])

            # Add records to the collection
            if isinstance(records, dict):
                all_records.update(records)
            else:
                all_records.update(
                    {str(i): record for i, record in enumerate(records, start=i)}  # type: ignore
                )

        return all_records

    except Exception as e:
        msg = f"Failed to fetch PubMed details: {e!s}"
        logger.opt(exception=True).error(msg)
        raise


def extract_article_info(article: dict[str, Any]) -> dict[str, str]:
    """Extracts key information from a PubMed article record.

    Example:
        >>> article = {
        ...     "MedlineCitation": {
        ...         "PMID": "39826015",
        ...         "Article": {
        ...             "ArticleTitle": "Progress report on multiple endocrine neoplasia type 1.",
        ...             "Abstract": {
        ...                 "AbstractText": [
        ...                     "Multiple endocrine neoplasia type 1 (MEN1) syndrome is an "
        ...                     "autosomal dominant disorder caused by a germline pathogenic "
        ...                     "variant in the MEN1 tumor suppressor gene. Patients with MEN1 "
        ...                     "have a high risk for primary hyperparathyroidism (PHPT) with a "
        ...                     "penetrance of nearly 100%, pituitary adenomas (PitAd) in 40% of "
        ...                     "patients, and neuroendocrine neoplasms (NEN) of the pancreas "
        ...                     "with a mainline mortality of 40%..."
        ...                 ]
        ...             },
        ...             "Journal": {
        ...                 "Title": "Familial cancer",
        ...             },
        ...             "JournalIssue": {
        ...                 "PubDate": {
        ...                     "Year": "2024",
        ...                 },
        ...             },
        ...         },
        ...     },
        ...     "PubmedData": {
        ...         "ArticleIdList": [
        ...             {
        ...                 "Id": "PMC11753851",
        ...                 "IdType": "pmc",
        ...             },
        ...             {
        ...                 "Id": "10.1155/jimr/5845167",
        ...                 "IdType": "doi",
        ...             },
        ...         ],
        ...     },
        ... }
        >>> extract_article_info(article)
        {'abstract': 'Multiple endocrine neoplasia type 1 (MEN1) syndrome is an '
                    'autosomal dominant disorder caused by a germline pathogenic '
                    'variant in the MEN1 tumor suppressor gene. Patients with MEN1 '
                    'have a high risk for primary hyperparathyroidism (PHPT) with a '
                    'penetrance of nearly 100%, pituitary adenomas (PitAd) in 40% of '
                    'patients, and neuroendocrine neoplasms (NEN) of the pancreas '
                    '(40% of patients), duodenum, lung, and thymus. Increased '
                    'MEN1-related mortality is mainly related to duodenal-pancreatic '
                    'and thymic NEN. Management of PHPT differs from that of patients '
                    'with sporadic disease, as the surgical approach in MEN1-related '
                    'PHPT includes near-total or total parathyroidectomy because of '
                    'multigland hyperplasia in most patients and the consequent high '
                    'risk of recurrence. NEN management also differs from patients '
                    'with sporadic disease due to multiple synchronous and '
                    'metasynchronous neoplasms. In addition, the lifelong risk of '
                    'developing NEN requires special considerations to avoid '
                    "excessive surgeries and to minimize damage to the patient's "
                    'function and well-being. This progress report will outline '
                    'current insights into surveillance and management of the major '
                    'clinical manifestation of MEN1 syndrome in children and adults '
                    'with MEN1 diagnosis. In addition, we will discuss MEN1-like '
                    'clinical presentation with negative MEN1-genetic workup and '
                    'future clinical and research directions.',
        'journal': 'Familial cancer',
        'pmid': '39826015',
        'pmc': 'PMC11753851',
        'doi': '10.1155/jimr/5845167',
        'title': 'Progress report on multiple endocrine neoplasia type 1.',
        'year': '2025'}

    """
    try:
        pmc = doi = ""
        for v in article["PubmedData"]["ArticleIdList"]:
            if v.attributes.get("IdType") == "pmc":
                pmc = str(v)
            elif v.attributes.get("IdType") == "doi":
                doi = str(v)

        medline = article["MedlineCitation"]
        article_info = medline["Article"]

        # Handle PMID correctly (it's a string element)
        pmid = str(medline["PMID"])

        # Handle abstract (might be a list of StringElements)
        abstract_text = article_info.get("Abstract", {}).get(
            "AbstractText", ["No abstract"]
        )
        if isinstance(abstract_text, list):
            abstract = "".join(str(text) for text in abstract_text)
        else:
            abstract = str(abstract_text)

        return {
            "pmid": pmid,
            "pmc": pmc,
            "doi": doi,
            "title": str(article_info.get("ArticleTitle", "No title")),
            "abstract": abstract,
            "journal": str(article_info.get("Journal", {}).get("Title", "No journal")),
            "year": str(
                article_info.get("JournalIssue", {})
                .get("PubDate", {})
                .get("Year", "No year")
            ),
        }

    except KeyError:
        msg = "Could not extract article information"
        logger.opt(exception=True).error(msg)
        raise
