import streamlit as st
import asyncio
import time


def cpu_bound_work(text: str) -> str:
    """Simulate CPU-bound work that should run in a thread."""
    time.sleep(3)  # Blocking sleep to simulate work
    return f"Processed: {text} at {time.strftime('%H:%M:%S')}"


async def process_in_thread(text: str) -> str:
    """Process a single item in a separate thread using to_thread."""
    return await asyncio.to_thread(cpu_bound_work, text)


@st.fragment(run_every="0.5s")
def update_progress(pb_empty=None):
    if not st.session_state.pb_empty:
        st.write("No progress container found")
        return
    """Update progress in the UI using a fragment."""
    total = st.session_state.get("total", 0)
    current = st.session_state.get("current", 0)

    with st.session_state.pb_empty.container():
        st.progress(current / total)
        st.write(f"Processed {current} of {total} items")


async def process_batch(
    texts: list[str], progress_container, results_container
) -> list[str]:
    """Process multiple items using to_thread with progress updates."""
    results = []
    total = len(texts)
    st.session_state.total = total
    st.session_state.current = 0

    for i, text in enumerate(texts, 1):
        try:
            result = await process_in_thread(text)
            results.append(result)

            # Update progress
            st.session_state.current = i
            update_progress()

            # Display result
            with results_container:
                st.write(f"{i}. {result}")

        except Exception as e:
            st.error(f"Error processing item {i}: {str(e)}")

    return results


def main():
    st.title("Test 3: asyncio.to_thread with Fragments")
    st.write("""
    This test uses asyncio.to_thread to handle CPU-bound work in separate threads,
    combined with Streamlit fragments for UI updates.
    Watch the timestamps and progress updates.
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
                # Create containers for progress and results
                progress_container = st.empty()
                st.session_state.pb_empty = progress_container
                results_container = st.container()

                # Process items
                results = asyncio.run(
                    process_batch(items, progress_container, results_container)
                )

        except Exception as e:
            st.error(f"Error in main: {str(e)}")


if __name__ == "__main__":
    main()
