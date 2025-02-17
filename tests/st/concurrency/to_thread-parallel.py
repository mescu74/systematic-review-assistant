import streamlit as st
import asyncio
import time
from typing import List


def cpu_bound_work(text: str) -> str:
    """Simulate CPU-bound work that should run in a thread."""
    time.sleep(2)  # Blocking sleep to simulate work
    return f"Processed: {text} at {time.strftime('%H:%M:%S')}"


async def process_in_thread(text: str) -> str:
    """Process a single item in a separate thread using to_thread."""
    return await asyncio.to_thread(cpu_bound_work, text)


@st.fragment
def update_result(container, result: str):
    """Update individual result in the UI using a fragment."""
    with container:
        st.write(result)


async def process_batch(texts: list[str], results_container) -> list[str]:
    """Process multiple items concurrently using gather with to_thread."""
    # Create all thread tasks at once
    tasks = [asyncio.create_task(process_in_thread(text)) for text in texts]

    try:
        # Use gather to run all thread tasks concurrently
        results = await asyncio.gather(*tasks)

        # Update UI with results
        for i, result in enumerate(results, 1):
            update_result(results_container, f"{i}. {result}")

        return results
    except Exception as e:
        st.error(f"Error in process_batch: {str(e)}")
        return []


def main():
    st.title("Test 4: Concurrent to_thread with gather")
    st.write("""
    This test combines asyncio.gather with to_thread to achieve true concurrency for CPU-bound work.
    Watch the timestamps - they should be very close together if working correctly.
    """)

    # Input
    items = st.text_area(
        "Enter items to process (one per line)",
        value="Item 1\nItem 2\nItem 3\nItem 4\nItem 5",
        height=100,
    ).split("\n")
    items = [item.strip() for item in items if item.strip()]

    if st.button("Process Items"):
        try:
            with st.spinner("Processing..."):
                # Create container for results
                results_container = st.container()

                # Process items concurrently
                results = asyncio.run(process_batch(items, results_container))

                if results:
                    st.success(
                        f"Successfully processed {len(results)} items concurrently"
                    )

        except Exception as e:
            st.error(f"Error in main: {str(e)}")


if __name__ == "__main__":
    main()
