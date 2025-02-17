# Chain callback handler streams tokens, cost, nof calls, run duration (on_xxx_end)
# base AsyncCallbackHandler to not run in a thread an be able to write to UI.
# OpenAI is sync and would run in executor no cotrol over threads.
#
# workers types
# - sync, ran with to_thread
# - async coro/task in same loop
# - async coro/task in different loop
# - sync threadpool executor
# - async taskggroup
"""
from langchain_community.callbacks.manager import get_openai_callback


with get_openai_callback() as cb:
    response = agent_executor.invoke(
        {
            "input": "What's a hummingbird's scientific name and what's the fastest bird species?"
        }
    )


cb.total_cost: float
cb.total_tokens: int
cb.prompt_tokens: int
cb.prompt_tokens_cached: int
cb.completion_tokens: int
cb.reasoning_tokens: int
cb.successful_requests: int
"""

import streamlit as st
import asyncio
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from typing import List, Dict, Any
from datetime import datetime
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession
from langchain.chat_models import ChatOpenAI
from langchain.callbacks.base import BaseCallbackHandler
from functools import partial
from contextlib import asynccontextmanager


class MetricsCallback(BaseCallbackHandler):
    """Callback handler for LLM metrics."""

    def __init__(self, metrics_container):
        self.metrics_container = metrics_container
        self.start_time = datetime.now()
        self.total_tokens = 0
        self.completed_calls = 0

    def on_llm_end(self, *args, **kwargs):
        self.completed_calls += 1
        self.total_tokens += kwargs.get("token_usage", {}).get("total_tokens", 0)

        with self.metrics_container:
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Completed Calls", self.completed_calls)
            with col2:
                st.metric("Total Tokens", self.total_tokens)
            with col3:
                elapsed = (datetime.now() - self.start_time).total_seconds()
                st.metric("Elapsed (s)", f"{elapsed:.1f}")


class BatchProcessor:
    """Handles batch processing using modern async patterns."""

    def __init__(self, db_url: str):
        self.db_url = db_url
        self.metrics = st.empty()
        self.logs = st.empty()
        self.results_df = st.empty()
        self.results: list[dict[str, Any]] = []

        # Set up thread pool for CPU-bound work
        self.executor = ThreadPoolExecutor(max_workers=4, thread_name_prefix="compute")

    @asynccontextmanager
    async def get_db(self):
        """Async context manager for DB sessions."""
        async with AsyncSession(engine) as session:  # Your engine setup
            yield session

    def process_dataframe(self, results: list[dict]) -> pd.DataFrame:
        """CPU-bound pandas processing."""
        df = pd.DataFrame(results)
        # Your pandas processing here
        return df

    async def update_ui(self, df: pd.DataFrame):
        """Update the UI with new results."""
        with self.results_df:
            edited_df = st.data_editor(df, num_rows="dynamic", key=f"df_{len(df)}")
            if not df.equals(edited_df):
                # Handle edits
                pass

    async def process_batch(self, items: list[str], batch_size: int = 5):
        """Process items in batches using TaskGroup."""

        # Initialize LLM with metrics callback
        llm = ChatOpenAI(callbacks=[MetricsCallback(self.metrics)])

        async def process_item_batch(batch: list[str]):
            """Process a single batch of items."""
            try:
                # LLM Processing
                responses = await llm.agenerate([batch])
                results = [r.text for r in responses.generations[0]]

                # DB Operations
                async with self.get_db() as session:
                    # Your DB save logic here
                    await session.commit()

                # Update results
                self.results.extend([{"text": r} for r in results])

                # CPU-bound processing in thread pool
                df = await asyncio.to_thread(self.process_dataframe, self.results)

                # Update UI
                await self.update_ui(df)

                with self.logs:
                    st.write(f"Processed batch of {len(batch)} items")

            except Exception as e:
                st.error(f"Error processing batch: {str(e)}")

        # Process all batches concurrently using TaskGroup
        async with asyncio.TaskGroup() as tg:
            for i in range(0, len(items), batch_size):
                batch = items[i : i + batch_size]
                tg.create_task(process_item_batch(batch))


def main():
    st.title("Modern Batch Processing")

    items = st.text_area(
        "Enter items to process (one per line)",
        value="Query 1\nQuery 2\nQuery 3\nQuery 4\nQuery 5",
        height=100,
    ).split("\n")
    items = [item.strip() for item in items if item.strip()]

    if st.button("Process Batch"):
        processor = BatchProcessor(db_url="your_db_url")

        try:
            asyncio.run(processor.process_batch(items))
            st.success("Processing complete!")

        except Exception as e:
            st.error(f"Error: {str(e)}")

        finally:
            processor.executor.shutdown()


if __name__ == "__main__":
    main()
