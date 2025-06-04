import asyncio
import time

import streamlit as st


async def simulate_processing(text: str, delay: float = 1.0) -> str:
    """Simulate async processing with delay."""
    await asyncio.sleep(delay)  # Non-blocking sleep
    return f"Processed: {text} at {time.strftime('%H:%M:%S')}"


async def process_batch(texts: list[str]) -> list[str]:
    """Process multiple items concurrently using asyncio.gather."""
    tasks = [asyncio.create_task(simulate_processing(text)) for text in texts]
    try:
        results = await asyncio.gather(*tasks)
        return results
    except Exception as e:
        st.error(f"Error in process_batch: {e!s}")
        return []


def main():
    st.title("Test 1: asyncio.gather with Tasks")
    st.write("""
    This test uses asyncio.gather with create_task to process items concurrently.
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

                # Note: This is where we need to be careful with the event loop
                results = asyncio.run(process_batch(items))

                # Display results
                for i, result in enumerate(results, 1):
                    with results_container:
                        st.write(f"{i}. {result}")

        except Exception as e:
            st.error(f"Error in main: {e!s}")


if __name__ == "__main__":
    main()
