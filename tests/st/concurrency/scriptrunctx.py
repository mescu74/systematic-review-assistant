import asyncio
import threading
import time

import streamlit as st
from streamlit.runtime.scriptrunner import add_script_run_ctx, get_script_run_ctx


def cpu_bound_work(text: str) -> str:
    """Simulate CPU-bound work that needs access to Streamlit context.
    Now we can safely use st.* calls here because the context is properly attached.
    """
    # Access session state safely
    count = st.session_state.get("process_count", 0)
    st.session_state.process_count = count + 1

    time.sleep(2)  # Simulate work
    return f"Processed: {text} at {time.strftime('%H:%M:%S')} (Count: {count + 1})"


async def process_in_thread(text: str) -> str:
    """Process in thread while preserving Streamlit context."""
    # Get current script run context before creating thread
    ctx = get_script_run_ctx()
    result = None
    error = None

    def thread_worker():
        nonlocal result, error
        try:
            result = cpu_bound_work(text)
        except Exception as e:
            error = e

    # Create and configure thread with context
    thread = threading.Thread(target=thread_worker)
    add_script_run_ctx(thread, ctx)

    # Run thread via asyncio.to_thread to maintain concurrency
    await asyncio.to_thread(lambda: (thread.start(), thread.join()))

    if error:
        raise error
    return result


@st.fragment
def update_result(container, result: str):
    """Update individual result in the UI."""
    with container:
        st.write(result)


async def process_batch(texts: list[str], results_container) -> list[str]:
    """Process multiple items concurrently with proper thread context."""
    # Initialize counter in session state if not present
    if "process_count" not in st.session_state:
        st.session_state.process_count = 0

    # Create tasks for concurrent processing
    tasks = [asyncio.create_task(process_in_thread(text)) for text in texts]

    try:
        # Run all tasks concurrently
        results = await asyncio.gather(*tasks)

        # Update UI with results
        for i, result in enumerate(results, 1):
            update_result(results_container, f"{i}. {result}")

        return results
    except Exception as e:
        st.error(f"Error in process_batch: {e!s}")
        return []


def main():
    st.title("Test 5: Thread Context Management (Fixed)")
    st.write("""
    This test demonstrates proper Streamlit context management with threads.
    - Each thread maintains access to session state
    - UI updates work correctly from threads
    - Thread safety is maintained
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
                        f"Successfully processed {len(results)} items concurrently. "
                        f"Total processes: {st.session_state.process_count}"
                    )

        except Exception as e:
            st.error(f"Error in main: {e!s}")


if __name__ == "__main__":
    main()
