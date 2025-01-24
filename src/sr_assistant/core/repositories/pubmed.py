"""PubMed repository."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Any
from uuid import UUID

from loguru import logger

from sr_assistant.core.models.base import PubMedResult
from sr_assistant.step2.pubmed_integration import (
    extract_article_info,
)

if TYPE_CHECKING:
    from supabase import Client


class PubMedRepository:
    """PubMed repository."""

    def __init__(self, supabase: Client) -> None:
        # self.supabase: Client = st.session_state.supabase
        self.supabase = supabase

    def store_results(
        self, review_id: UUID, query: str, records: dict[str, Any]
    ) -> list[PubMedResult]:
        """Store PubMed search results."""
        try:
            results = []
            for record in records["PubmedArticle"]:
                try:
                    article_info = extract_article_info(record)
                    if not article_info:
                        continue
                    result = PubMedResult(
                        review_id=review_id,
                        query=query,
                        pmid=article_info["pmid"],
                        pmc=article_info["pmc"],
                        doi=article_info["doi"],
                        title=article_info["title"],
                        abstract=article_info["abstract"],
                        journal=article_info["journal"],
                        year=article_info["year"],
                    )
                    logger.info(f"Storing PubMed result: {result}")
                    results.append(result)
                except Exception as e:
                    # logger.warning(f"Failed to extract article info: {e!r}")
                    msg = f"Failed to extract article info from record: {record}"
                    logger.opt(exception=True).warning(msg)
                    continue

            if not results:
                return []

            # Convert list of models to JSON-compatible data
            data = [
                json.loads(result.model_dump_json(exclude={"created_at", "updated_at"}))
                for result in results
            ]

            ret = self.supabase.table("pubmed_results").insert(data).execute()
            logger.info(f"Stored {len(ret.data)} PubMed results for review {review_id}")

            # Supabase returns UUID as string
            return [
                # PubMedResult(
                #    **{**r, "review_id": UUID(r["review_id"]), "id": UUID(r["id"])}
                # )
                PubMedResult.model_validate(r)
                for r in ret.data
            ]

        except Exception as e:
            logger.error(f"Failed to store search results: {e!r}")
            raise

    def get_review_results(self, review_id: UUID) -> list[PubMedResult]:
        """Get all results for a review."""
        try:
            ret = (
                self.supabase.table("pubmed_results")
                .select("*")
                .eq("review_id", str(review_id))
                .execute()
            )
            return [
                PubMedResult(
                    **{**r, "review_id": UUID(r["review_id"]), "id": UUID(r["id"])}
                )
                for r in ret.data
            ]
        except Exception as e:
            logger.error(f"Failed to get results for review {review_id}: {e!r}")
            raise
