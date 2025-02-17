import streamlit as st
import asyncio
import time
from typing import List


async def simulate_processing(text: str, delay: float = 1.0) -> str:
    """Simulate async processing with delay."""
    await asyncio.sleep(delay)  # Non-blocking sleep
    return f"Processed: {text} at {time.strftime('%H:%M:%S')}"


async def process_batch(texts: list[str]) -> list[str]:
    """Process multiple items concurrently using TaskGroup."""
    results = []
    async with asyncio.TaskGroup() as tg:
        tasks = [tg.create_task(simulate_processing(text)) for text in texts]
        try:
            # Wait for tasks to complete and collect results
            results = [await task for task in tasks]
        except* Exception as e:
            st.error(f"Error in TaskGroup: {str(e)}")
    return results


def main():
    st.title("Test 2: asyncio.TaskGroup")
    st.write("""
    This test uses asyncio.TaskGroup for concurrent processing.
    TaskGroup provides structured concurrency and better error handling.
    Watch the timestamps to verify concurrent execution.
    """)

    # Input
    items = st.text_area(
        "Enter items to process (one per line)",
        value="Item 1\nItem 2\nItem 3",
        height=100,
    ).split("\n")
    items = [item.strip() for item in items if item.strip()]

    if st.button("Process Items"):
        try:
            with st.spinner("Processing..."):
                # Create progress placeholder
                progress_text = st.empty()
                results_container = st.container()

                # Process items
                results = asyncio.run(process_batch(items))

                # Display results
                for i, result in enumerate(results, 1):
                    with results_container:
                        st.write(f"{i}. {result}")

        except Exception as e:
            st.error(f"Error in main: {str(e)}")


if __name__ == "__main__":
    main()
