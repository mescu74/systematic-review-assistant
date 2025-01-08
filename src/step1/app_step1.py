from __future__ import annotations

import datetime
import json
import os
import uuid

import streamlit as st
from dotenv import load_dotenv
from langchain_core.messages import HumanMessage
from langchain_core.runnables import RunnableConfig
from langchain_openai import ChatOpenAI
from langgraph.checkpoint.memory import MemorySaver
from langgraph.graph import START, MessagesState
from langgraph.graph.state import CompiledStateGraph, StateGraph

# Load environment variables from .env if present
load_dotenv()


def init_agent_config() -> RunnableConfig:
    if "agent_config" in st.session_state:
        return st.session_state.agent_config
    config = {"recursion_limit": 1000, "configurable": {"thread_id": uuid.uuid4().hex}}
    config = RunnableConfig(**config)
    st.session_state.agent_config = config
    return st.session_state.agent_config


def init_suggestion_agent() -> CompiledStateGraph:
    """Setup a simple LangGraph agent."""
    # Initialize LLM
    llm = ChatOpenAI(model="gpt-4o", temperature=0.0)

    # Define a new graph
    workflow = StateGraph(state_schema=MessagesState)

    # Define the function that calls the model
    def call_model(state: MessagesState) -> dict:
        response = llm.invoke(state["messages"])
        return {"messages": response}

    # Define the two nodes we will cycle between
    workflow.add_node("model", call_model)
    workflow.add_edge(START, "model")

    # Add memory
    memory = MemorySaver()

    # The thread id is a unique key that identifies
    # this particular conversation.
    # We'll just generate a random uuid here.
    st.session_state["suggestions_agent"] = workflow.compile(
        checkpointer=memory
    ).with_config(init_agent_config())
    return st.session_state["suggestions_agent"]


# Set up page
st.set_page_config(page_title="AI-Assisted Systematic Review - Step 1", layout="wide")
st.title("Step 1: Define Research Question & Criteria")

# Initialize session state variables
if "research_question" not in st.session_state:
    st.session_state["research_question"] = ""
if "inclusion_criteria" not in st.session_state:
    st.session_state["inclusion_criteria"] = ""
if "exclusion_criteria" not in st.session_state:
    st.session_state["exclusion_criteria"] = ""
if "llm_suggestions" not in st.session_state:
    st.session_state["llm_suggestions"] = ""
if "suggestions_agent" not in st.session_state:
    st.session_state["suggestions_agent"] = init_suggestion_agent()
# LangGraph agent state holding full message history
if "suggestions_agent_state" not in st.session_state:
    st.session_state.suggestions_agent_state = {}
    st.session_state.suggestions_agent_state["messages"] = []

# Layout for input fields and LLM feedback
col1, col2 = st.columns(2)

with col1:
    st.subheader("Enter Your Criteria")
    st.session_state["research_question"] = st.text_area(
        "Research Question",
        st.session_state["research_question"],
        placeholder="E.g. 'Does intervention X improve outcome Y in population Z?'",
    )
    st.session_state["inclusion_criteria"] = st.text_area(
        "Inclusion Criteria",
        st.session_state["inclusion_criteria"],
        placeholder="E.g. Adults >18, Intervention X, Outcome Y, RCTs, English",
    )
    st.session_state["exclusion_criteria"] = st.text_area(
        "Exclusion Criteria",
        st.session_state["exclusion_criteria"],
        placeholder="E.g. Children <18, Non-human studies, Case reports",
    )

    validate_button = st.button("Validate & Suggest Improvements")

with col2:
    st.subheader("LLM Suggestions")
    if st.session_state["llm_suggestions"]:
        st.markdown(st.session_state["llm_suggestions"])
    else:
        st.write("No suggestions yet. Click the 'Validate' button to get feedback.")

# LLM Interaction
if validate_button:
    user_input_summary = f"""\
    Research Question: {st.session_state['research_question']}
    Inclusion Criteria: {st.session_state['inclusion_criteria']}
    Exclusion Criteria: {st.session_state['exclusion_criteria']}
    """

    prompt = f"""\
    The user provided the following systematic review setup:

    {user_input_summary}

    Please:
    1. Check if the criteria are clear and unambiguous.
    2. Suggest any improvements or additional details that might help during the screening process.
    3. Identify any contradictions or points needing clarification.
    """

    st.session_state["suggestions_agent_state"]["messages"].append(
        HumanMessage(content=prompt)
    )
    st.session_state["suggestion_agent_state"] = st.session_state[
        "suggestions_agent"
    ].invoke(
        st.session_state["suggestions_agent_state"],
        config=st.session_state.agent_config,
    )
    st.session_state["llm_suggestions"] = st.session_state["suggestion_agent_state"][
        "messages"
    ][-1].content
    st.rerun()

st.subheader("Finalize & Save")
st.write("Once you are happy with the criteria, click the 'Save Criteria' button.")

save_button = st.button("Save Criteria")
if save_button:
    criteria_data = {
        "research_question": st.session_state["research_question"],
        "inclusion_criteria": st.session_state["inclusion_criteria"],
        "exclusion_criteria": st.session_state["exclusion_criteria"],
        "timestamp": datetime.datetime.now().isoformat(),
    }

    # Create directories if they don't exist
    os.makedirs("criteria", exist_ok=True)
    os.makedirs("logs", exist_ok=True)

    filename = (
        f"criteria/criteria_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    )
    with open(filename, "w") as f:
        json.dump(criteria_data, f, indent=2)

    # Log the event
    with open("logs/criteria_setup_log.txt", "a") as log_file:
        log_file.write(f"---\n{datetime.datetime.now().isoformat()}\n")
        log_file.write(f"Final Criteria:\n{json.dumps(criteria_data, indent=2)}\n")
        if st.session_state["llm_suggestions"]:
            log_file.write("LLM Suggestions:\n")
            log_file.write(st.session_state["llm_suggestions"] + "\n")

    st.success(f"Criteria saved to {filename}")
