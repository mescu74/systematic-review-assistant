from __future__ import annotations

from enum import StrEnum, auto


class CitationFormat(StrEnum):
    """Standardized citation formats supported by academic databases.

    Different disciplines and journals prefer different citation formats. This
    classification helps ensure consistent citation export across databases.

    Examples:
        >>> ref = Reference(format=CitationFormat.VANCOUVER)
        >>> print(ref.format())
        'Smith J, et al. Cancer Research. 2023;12(3):45-50.'
        >>> ref.format = CitationFormat.APA
        >>> print(ref.format())
        'Smith, J., & Jones, M. (2023). Cancer research. Journal, 12(3), 45-50.'

    See Also:
        `Citation Guide <https://www.nlm.nih.gov/bsd/uniform_requirements.html>`_
    """

    VANCOUVER = auto()
    """International Committee of Medical Journal Editors (ICMJE) format.

    Standard format for medical literature, used by PubMed and most medical journals.
    """

    APA = auto()
    """American Psychological Association format.

    Used in psychology, education, and social sciences. Standard for PsycINFO.
    """

    HARVARD = auto()
    """Harvard citation format.

    Common in business and humanities. Used by many European journals.
    """

    IEEE = auto()
    """Institute of Electrical and Electronics Engineers format.

    Standard for engineering and computer science publications.
    """

    CHICAGO = auto()
    """Chicago Manual of Style format.

    Common in humanities and some social sciences.
    """

    MLA = auto()
    """Modern Language Association format.

    Standard format for humanities, especially literature and linguistics.
    """

    @property
    def description(self) -> str:
        """Returns detailed description of the citation format.

        Returns:
            str: Format description and typical use cases.
        """
        descriptions = {
            CitationFormat.VANCOUVER: (
                "Numerical citation format used in medical and scientific "
                "publications. Follows ICMJE recommendations."
            ),
            CitationFormat.APA: (
                "Author-date format used in social sciences. Currently in 7th "
                "edition. Features DOIs and URLs for electronic sources."
            ),
            CitationFormat.HARVARD: (
                "Author-date format popular in universities. Note that 'Harvard' "
                "has several variants based on institution."
            ),
            CitationFormat.IEEE: (
                "Numerical format with specific conventions for technical "
                "documents. Used in engineering fields."
            ),
            CitationFormat.CHICAGO: (
                "Offers both author-date and notes-bibliography systems. "
                "Currently in 17th edition."
            ),
            CitationFormat.MLA: (
                "Author-page format used in humanities. Currently in 9th "
                "edition. Emphasizes print sources."
            ),
        }
        return descriptions[self]

    @property
    def example(self) -> str:
        """Returns an example citation in this format.

        Returns:
            str: Example formatted citation.
        """
        examples = {
            CitationFormat.VANCOUVER: (
                "1. Smith J, Jones M, Liu X. Cancer patterns. N Engl J Med. "
                "2023;380(1):23-34. doi: 10.1056/nejm.2023.1234"
            ),
            CitationFormat.APA: (
                "Smith, J., Jones, M., & Liu, X. (2023). Cancer patterns in "
                "global populations. New England Journal of Medicine, 380(1), "
                "23-34. https://doi.org/10.1056/nejm.2023.1234"
            ),
            CitationFormat.HARVARD: (
                "Smith, J., Jones, M. and Liu, X., 2023. Cancer patterns in "
                "global populations. New England Journal of Medicine, 380(1), "
                "pp.23-34."
            ),
            CitationFormat.IEEE: (
                '[1] J. Smith, M. Jones, and X. Liu, "Cancer patterns in '
                'global populations," N. Engl. J. Med., vol. 380, no. 1, '
                "pp. 23-34, Jan. 2023."
            ),
            CitationFormat.CHICAGO: (
                'Smith, John, Mary Jones, and Xiao Liu. "Cancer Patterns in '
                'Global Populations." New England Journal of Medicine 380, '
                "no. 1 (January 2023): 23-34."
            ),
            CitationFormat.MLA: (
                'Smith, John, et al. "Cancer Patterns in Global Populations." '
                "New England Journal of Medicine, vol. 380, no. 1, Jan. 2023, "
                "pp. 23-34."
            ),
        }
        return examples[self]
