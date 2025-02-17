import streamlit as st
import asyncio
import threading
import time
from typing import List, Any
from streamlit.runtime.scriptrunner import add_script_run_ctx, get_script_run_ctx
from queue import Queue
from concurrent.futures import ThreadPoolExecutor


class AsyncWorkerPool:
    """Manages a dedicated event loop in a separate thread for specific work types."""

    def __init__(self, name: str):
        self.name = name
        self.loop: asyncio.AbstractEventLoop | None = None
        self.thread: threading.Thread | None = None
        self.queue: asyncio.Queue | None = None
        self.running = False
        self.ctx = get_script_run_ctx()  # Capture context from creating thread

    async def worker(self):
        """Process tasks from the queue."""
        while self.running:
            try:
                task = await self.queue.get()
                if task is None:  # Shutdown signal
                    break

                func, args = task
                try:
                    result = await func(*args)
                    st.session_state[f"last_{self.name}_result"] = result
                except Exception as e:
                    st.error(f"Error in {self.name} worker: {str(e)}")
                finally:
                    self.queue.task_done()

            except asyncio.CancelledError:
                break

    def run_loop(self):
        """Run the event loop in this thread."""
        # Set up the loop
        self.loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self.loop)
        self.queue = asyncio.Queue()

        # Start worker task
        self.running = True
        worker_task = self.loop.create_task(self.worker())

        try:
            self.loop.run_forever()
        finally:
            # Clean up
            worker_task.cancel()
            self.loop.run_until_complete(worker_task)
            self.loop.close()

    def start(self):
        """Start the worker thread with proper context."""
        self.thread = threading.Thread(target=self.run_loop, name=f"{self.name}_worker")
        add_script_run_ctx(self.thread, self.ctx)  # Attach Streamlit context
        self.thread.start()

    async def stop(self):
        """Stop the worker thread cleanly."""
        if self.running and self.loop:
            self.running = False
            await self.queue.put(None)  # Signal shutdown
            self.loop.call_soon_threadsafe(self.loop.stop)
            self.thread.join()

    async def submit(self, func, *args):
        """Submit a task to this worker's queue."""
        if self.queue:
            await self.queue.put((func, args))


# Simulate different types of blocking work
async def db_operation(item: str) -> str:
    """Simulate a database operation."""
    await asyncio.sleep(1)  # Simulate IO
    return f"DB processed {item} at {time.strftime('%H:%M:%S')}"


async def cpu_work(item: str) -> str:
    """Simulate CPU-intensive work."""
    time.sleep(1)  # Actually block
    return f"CPU processed {item} at {time.strftime('%H:%M:%S')}"


async def process_batches(
    db_items: list[str], cpu_items: list[str], results_container
) -> None:
    """Process items using different worker pools."""

    # Create worker pools
    db_pool = AsyncWorkerPool("db")
    cpu_pool = AsyncWorkerPool("cpu")

    # Start workers
    db_pool.start()
    cpu_pool.start()

    try:
        # Submit work to both pools
        db_tasks = [db_pool.submit(db_operation, item) for item in db_items]
        cpu_tasks = [cpu_pool.submit(cpu_work, item) for item in cpu_items]

        # Wait for all tasks to be submitted
        await asyncio.gather(*(db_tasks + cpu_tasks))

        # Wait for queues to be processed
        if db_pool.queue:
            await db_pool.queue.join()
        if cpu_pool.queue:
            await cpu_pool.queue.join()

        # Update results
        with results_container:
            st.write("### DB Results")
            if "last_db_result" in st.session_state:
                st.write(st.session_state.last_db_result)

            st.write("### CPU Results")
            if "last_cpu_result" in st.session_state:
                st.write(st.session_state.last_cpu_result)

    finally:
        # Clean shutdown
        await db_pool.stop()
        await cpu_pool.stop()


def main():
    st.title("Test 6: Multiple Event Loops")
    st.write("""
    This test demonstrates running multiple event loops in separate threads:
    - One loop for DB/IO operations
    - One loop for CPU-bound work
    - Each with proper Streamlit context
    Each loop processes work independently but can update the UI.
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
                results_container = st.container()
                asyncio.run(process_batches(db_items, cpu_items, results_container))

        except Exception as e:
            st.error(f"Error in main: {str(e)}")


if __name__ == "__main__":
    main()
