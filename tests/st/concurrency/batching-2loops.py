import streamlit as st
import asyncio
import threading
import time
from typing import List
from streamlit.runtime.scriptrunner import add_script_run_ctx, get_script_run_ctx


class WorkerThread:
    """Manages work in a separate thread with its own event loop."""

    def __init__(self, name: str):
        self.name = name
        self.loop = asyncio.new_event_loop()
        self.thread = threading.Thread(
            target=self._run_loop, name=f"{self.name}_worker"
        )
        self.ctx = get_script_run_ctx()

    def _run_loop(self):
        """Sets up and runs the event loop in this thread."""
        asyncio.set_event_loop(self.loop)
        self.loop.run_forever()

    def start(self):
        """Start the worker thread with proper context."""
        add_script_run_ctx(self.thread, self.ctx)
        self.thread.start()

    def stop(self):
        """Stop the worker thread cleanly."""
        self.loop.call_soon_threadsafe(self.loop.stop)
        self.thread.join()

    async def run_task(self, coro):
        """Run a coroutine in this worker's loop."""
        return await asyncio.wrap_future(
            asyncio.run_coroutine_threadsafe(coro, self.loop)
        )


# Simulate different types of work
async def db_operation(item: str, status_container) -> str:
    """Simulate a database operation."""
    with status_container:
        st.write(f"Starting DB operation for {item}")
    await asyncio.sleep(1)  # Simulate IO
    result = f"DB processed {item} at {time.strftime('%H:%M:%S')}"
    with status_container:
        st.write(f"Completed: {result}")
    return result


async def cpu_work(item: str, status_container) -> str:
    """Simulate CPU-intensive work."""
    with status_container:
        st.write(f"Starting CPU work for {item}")
    await asyncio.sleep(1)  # Using sleep instead of blocking for demo
    result = f"CPU processed {item} at {time.strftime('%H:%M:%S')}"
    with status_container:
        st.write(f"Completed: {result}")
    return result


async def process_batches(
    db_items: list[str], cpu_items: list[str], status_container, results_container
) -> None:
    """Process items using different workers."""

    # Create and start workers
    db_worker = WorkerThread("db")
    cpu_worker = WorkerThread("cpu")

    db_worker.start()
    cpu_worker.start()

    try:
        # Create tasks for both types of work
        db_tasks = [
            db_worker.run_task(db_operation(item, status_container))
            for item in db_items
        ]
        cpu_tasks = [
            cpu_worker.run_task(cpu_work(item, status_container)) for item in cpu_items
        ]

        # Wait for all tasks to complete
        db_results = await asyncio.gather(*db_tasks)
        cpu_results = await asyncio.gather(*cpu_tasks)

        # Show final results
        with results_container:
            st.write("### DB Results")
            for result in db_results:
                st.write(result)

            st.write("### CPU Results")
            for result in cpu_results:
                st.write(result)

    finally:
        # Clean shutdown
        db_worker.stop()
        cpu_worker.stop()


def main():
    st.title("Test 6: Multiple Event Loops (Fixed)")
    st.write("""
    This test demonstrates running multiple event loops in separate threads:
    - One loop for DB/IO operations
    - One loop for CPU-bound work
    - Real-time status updates
    """)

    # Input
    col1, col2 = st.columns(2)
    with col1:
        db_items = st.text_area(
            "DB Items (one per line)", value="Query 1\nQuery 2\nQuery 3", height=100
        ).split("\n")
        db_items = [item.strip() for item in db_items if item.strip()]

    with col2:
        cpu_items = st.text_area(
            "CPU Items (one per line)", value="Task 1\nTask 2\nTask 3", height=100
        ).split("\n")
        cpu_items = [item.strip() for item in cpu_items if item.strip()]

    if st.button("Process All"):
        try:
            with st.spinner("Processing..."):
                status_container = st.empty()
                results_container = st.container()

                asyncio.run(
                    process_batches(
                        db_items, cpu_items, status_container, results_container
                    )
                )

                st.success("All processing complete!")

        except Exception as e:
            st.error(f"Error in main: {str(e)}")


if __name__ == "__main__":
    main()
