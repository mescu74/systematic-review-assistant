"""Suggestion agent for systematic review criteria.

Needs a rewrite.
"""

from __future__ import annotations

from typing import cast

import uuid6
from langchain_core.messages import AIMessage, BaseMessage, HumanMessage
from langchain_core.runnables import RunnableConfig
from langchain_openai import ChatOpenAI
from langgraph.checkpoint.memory import MemorySaver
from langgraph.graph import START
from langgraph.graph.message import MessagesState
from langgraph.graph.state import CompiledStateGraph, StateGraph
from loguru import logger

from sr_assistant.core.models import SystematicReview


class SuggestionAgent:
    """Agent for providing suggestions on systematic review criteria."""

    def __init__(self, model: str = "gpt-4o", temperature: float = 0.0) -> None:
        """Initialize the agent with model configuration."""
        self.llm = ChatOpenAI(model=model, temperature=temperature)
        self.workflow = self._create_workflow()
        self.memory = MemorySaver()  # TODO: postgres checkpointer
        self.agent = self._compile_agent()
        # if st.session_state.criteria_messages:
        #     self.messages: list[BaseMessage] = st.session_state.criteria_messages
        # else:
        #     st.session_state.criteria_messages: list[BaseMessage] = []
        #     self.messages: list[BaseMessage] = []
        self.messages: list[BaseMessage] = []
        self.config = self._create_config()

    def _create_workflow(self) -> StateGraph:
        """Create the agent workflow."""
        workflow = StateGraph(state_schema=MessagesState)

        def call_model(state: MessagesState) -> dict[str, AIMessage]:
            """Call the model with the current state."""
            # FIXME: not AIMessage, MessaagesState?
            response = cast(AIMessage, self.llm.invoke(state["messages"]))
            return {"messages": response}

        workflow.add_node("model", call_model)
        workflow.add_edge(START, "model")
        return workflow

    def _compile_agent(self) -> CompiledStateGraph:
        """Compile the agent with memory."""
        # TODO: postgres
        return self.workflow.compile(checkpointer=self.memory)

    def _create_config(self) -> RunnableConfig:
        """Create agent configuration."""
        criteria_agent_thread_id = str(uuid6.uuid7())
        return RunnableConfig(
            recursion_limit=1000, configurable={"thread_id": criteria_agent_thread_id}
        )

    def get_suggestions(self, review: SystematicReview) -> str:
        """Get suggestions for the given criteria."""
        prompt = f"""\
        The user provided the following systematic review protocol:

        Background: {review.background}
        Research Question: {review.research_question}
        Inclusion Criteria: {review.inclusion_criteria}
        Exclusion Criteria: {review.exclusion_criteria}

        Please:
        1. Check if the criteria are clear and unambiguous.
        2. Suggest any improvements or additional details that might help during the screening process.
        3. Identify any contradictions or points needing clarification.
        """
        # Add message to history
        self.messages.append(HumanMessage(content=prompt))

        # Update agent state
        state = {"messages": self.messages}
        result = self.agent.invoke(state, config=self.config)

        # Store response
        ai_message = result["messages"][-1]
        if isinstance(ai_message, AIMessage):
            self.messages.append(ai_message)
            return cast(str, ai_message.content)
        msg = "Unexpected message type: {!r}"
        logger.error(msg, type(ai_message))
        raise ValueError(msg.format(type(ai_message)))

    def get_message_history(self) -> list[BaseMessage]:
        """Get the full message history."""
        return self.messages.copy()

    def clear_history(self) -> None:
        """Clear the message history."""
        self.messages = []
