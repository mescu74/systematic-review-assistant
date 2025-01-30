"""PubMed repository."""

from __future__ import annotations

import json
from typing import TYPE_CHECKING, Any
from uuid import UUID

from loguru import logger

from sr_assistant.core.models import AbstractScreeningResult, PubMedResult
from sr_assistant.core.schemas import ScreeningResponse
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
                #PubMedResult(
                #    **{**r, "review_id": UUID(r["review_id"]), "id": UUID(r["id"])}
                #)
                PubMedResult.model_validate(r)
                for r in ret.data
            ]
        except Exception as e:
            logger.error(f"Failed to get results for review {review_id}: {e!r}")
            raise
"""Screening repositories.

These are the repositories for the screening decisions made by the reviewers.
"""

class ScreeningAbstractRepository:
    """Repository for managing abstract screening decisions."""

    def __init__(self, supabase: Client) -> None:
        """Initialize the repository with a Supabase client.

        Args:
            supabase: Initialized Supabase client
        """
        self.supabase = supabase

    def _map_screening_response_to_result(
        self, review_id: UUID, search_result_id: UUID, response: ScreeningResponse
    ) -> AbstractScreeningResult:
        """Map ScreeningResponse to AbstractScreeningResult.

        Args:
            review_id: ID of the systematic review
            search_result_id: ID of the PubMed search result
            response: Screening response from the model

        Returns:
            AbstractScreeningResult model instance
        """
        return AbstractScreeningResult(
            review_id=review_id,
            search_result_id=search_result_id,
            decision=response.decision,
            confidence_score=response.confidence_score,
            rationale=response.rationale,
            extracted_quotes=response.extracted_quotes,
            exclusion_reason_categories=response.exclusion_reason_categories,
        )

    def create_screening_result(
        self,
        review_id: UUID,
        search_result_id: UUID,
        response: ScreeningResponse,
    ) -> AbstractScreeningResult:
        """Create a new abstract screening result.

        Args:
            review_id: ID of the systematic review
            search_result_id: ID of the PubMed search result being screened
            response: Screening response from the model

        Returns:
            Created AbstractScreeningResult
        """
        result = self._map_screening_response_to_result(
            review_id=review_id,
            search_result_id=search_result_id,
            response=response,
        )

        try:
            data = json.loads(result.model_dump_json(exclude={"created_at", "updated_at"}))
            ret = self.supabase.table("abstract_screening_results").insert(data).execute()
            return AbstractScreeningResult.model_validate(ret.data[0])
        except Exception as e:
            logger.error(f"Failed to create screening result: {e}")
            raise

    def get_screening_results(self, review_id: UUID) -> list[AbstractScreeningResult]:
        """Get all screening results for a review.

        Args:
            review_id: ID of the systematic review

        Returns:
            List of AbstractScreeningResult for the review
        """
        try:
            ret = (
                self.supabase.table("abstract_screening_results")
                .select("*")
                .eq("review_id", str(review_id))
                .execute()
            )
            return [AbstractScreeningResult.model_validate(item) for item in ret.data]
        except Exception as e:
            logger.error(f"Failed to get screening results: {e}")
            raise

    def get_screening_result(
        self, review_id: UUID, search_result_id: UUID
    ) -> AbstractScreeningResult | None:
        """Get screening result for a specific search result.

        Args:
            review_id: ID of the systematic review
            search_result_id: ID of the PubMed search result

        Returns:
            AbstractScreeningResult if found, None otherwise
        """
        try:
            ret = (
                self.supabase.table("abstract_screening_results")
                .select("*")
                .eq("review_id", str(review_id))
                .eq("search_result_id", str(search_result_id))
                .execute()
            )
            return (
                AbstractScreeningResult.model_validate(ret.data[0])
                if ret.data
                else None
            )
        except Exception as e:
            logger.error(f"Failed to get screening result: {e}")
            raise

    def update_screening_result(
        self,
        review_id: UUID,
        search_result_id: UUID,
        response: ScreeningResponse,
    ) -> AbstractScreeningResult:
        """Update an existing screening result.

        Args:
            review_id: ID of the systematic review
            search_result_id: ID of the PubMed search result
            response: Updated screening response from the model

        Returns:
            Updated AbstractScreeningResult

        Raises:
            Exception: If the screening result doesn't exist
        """
        result = self._map_screening_response_to_result(
            review_id=review_id,
            search_result_id=search_result_id,
            response=response,
        )

        try:
            data = json.loads(result.model_dump_json(exclude={"created_at", "updated_at"}))
            ret = (
                self.supabase.table("abstract_screening_results")
                .update(data)
                .eq("review_id", str(review_id))
                .eq("search_result_id", str(search_result_id))
                .execute()
            )
            if not ret.data:
                raise Exception("Screening result not found")
            return AbstractScreeningResult.model_validate(ret.data[0])
        except Exception as e:
            logger.error(f"Failed to update screening result: {e}")
            raise

    def delete_screening_result(
        self, review_id: UUID, search_result_id: UUID
    ) -> None:
        """Delete a screening result.

        Args:
            review_id: ID of the systematic review
            search_result_id: ID of the PubMed search result
        """
        try:
            self.supabase.table("abstract_screening_results").delete().eq(
                "review_id", str(review_id)
            ).eq("search_result_id", str(search_result_id)).execute()
        except Exception as e:
            logger.error(f"Failed to delete screening result: {e}")
            raise
