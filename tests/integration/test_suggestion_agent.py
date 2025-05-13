import uuid

import pytest

from sr_assistant.core.models import SystematicReview
from sr_assistant.step1.suggestion_agent import SuggestionAgent


@pytest.fixture(scope="module")  # Use module scope if agent setup is expensive
def agent() -> SuggestionAgent:
    """Fixture to provide a SuggestionAgent instance."""
    # Add necessary setup here if needed, e.g., environment variables
    return SuggestionAgent()


@pytest.mark.integration
def test_agent_suggestions(agent: SuggestionAgent):
    """Test that the agent provides suggestions."""
    # Create a dummy review object
    mock_review = SystematicReview(
        id=uuid.uuid4(),  # Add required fields for model instantiation
        research_question="RQ...",
        exclusion_criteria="None",  # Add required fields for model instantiation
        background="Background info...",  # Add required fields for model instantiation
    )
    suggestions = agent.get_suggestions(review=mock_review)  # Pass the review object

    assert suggestions is not None
    assert isinstance(suggestions, dict)
    assert "pico" in suggestions
    assert "raw_response" in suggestions
    assert len(suggestions.get("raw_response", "")) > 10, (
        "Response should contain some text"
    )

    # History tests removed as agent has no memory
