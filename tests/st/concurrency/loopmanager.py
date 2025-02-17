import asyncio
import threading
import time
from typing import Set
from dataclasses import dataclass, field
import streamlit as st


@dataclass
class LoopManager:
    """Manages an event loop running in a separate thread."""

    name: str
    loop: asyncio.AbstractEventLoop = field(default_factory=asyncio.new_event_loop)
    thread: threading.Thread = None
    running: bool = False
    tasks: set[asyncio.Task] = field(default_factory=set)

    def start(self):
        """Start the loop in a new thread."""
        self.running = True
        self.thread = threading.Thread(target=self._run_loop, name=f"{self.name}_loop")
        self.thread.start()

    def _run_loop(self):
        """Run the event loop and handle cleanup."""
        print(f"Starting loop in {threading.current_thread().name}")
        asyncio.set_event_loop(self.loop)

        try:
            self.loop.run_forever()
        finally:
            try:
                # Cancel pending tasks
                self._cancel_tasks()

                # Run loop one last time to complete cancellation
                self.loop.run_until_complete(self._cleanup_tasks())
            finally:
                print(f"Closing loop in {threading.current_thread().name}")
                self.loop.close()

    async def _cleanup_tasks(self):
        """Wait for tasks to complete their cancellation."""
        tasks = list(self.tasks)
        for task in tasks:
            if not task.cancelled():
                try:
                    await task
                except asyncio.CancelledError:
                    pass
                except Exception as e:
                    print(f"Task error during cleanup: {e}")

    def _cancel_tasks(self):
        """Cancel all pending tasks."""
        for task in self.tasks:
            task.cancel()
        self.tasks.clear()

    def stop(self):
        """Stop the loop and wait for thread to finish."""
        if self.running:
            self.running = False
            self.loop.call_soon_threadsafe(self._cancel_tasks)
            self.loop.call_soon_threadsafe(self.loop.stop)
            self.thread.join()

    async def long_task(self, task_id: int, duration: int):
        """A task that takes time and can be cancelled."""
        try:
            print(f"Task {task_id} starting")
            for i in range(duration):
                await asyncio.sleep(1)
                print(f"Task {task_id}: {i+1}/{duration} seconds")
            print(f"Task {task_id} completed")
            return f"Task {task_id} result"

        except asyncio.CancelledError:
            print(f"Task {task_id} cancelled")
            raise

        finally:
            print(f"Task {task_id} cleanup")

    def submit_task(self, task_id: int, duration: int):
        """Submit a new task to the loop."""
        if not self.running:
            raise RuntimeError("Loop is not running")

        coro = self.long_task(task_id, duration)
        future = asyncio.run_coroutine_threadsafe(coro, self.loop)

        # Keep track of the task
        task = asyncio.Task(coro, loop=self.loop)
        self.tasks.add(task)
        task.add_done_callback(self.tasks.discard)

        return future

    def submit_task(self, task_id: int, duration: int):
        """Submit a new task to the loop."""
        if not self.running:
            raise RuntimeError("Loop is not running")

        coro = self.long_task(task_id, duration)

        # Create and track the task
        def create_task():
            task = asyncio.Task(coro, loop=self.loop)
            self.tasks.add(task)
            task.add_done_callback(self.tasks.discard)
            return task

        return self.loop.call_soon_threadsafe(create_task)


def main():
    st.title("Event Loop Lifecycle Demo")

    if "loop_manager" not in st.session_state:
        st.session_state.loop_manager = None

    # Start/Stop loop
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Start Loop"):
            if st.session_state.loop_manager is None:
                st.session_state.loop_manager = LoopManager("demo")
                st.session_state.loop_manager.start()
                st.success("Loop started!")

    with col2:
        if st.button("Stop Loop"):
            if st.session_state.loop_manager is not None:
                st.session_state.loop_manager.stop()
                st.session_state.loop_manager = None
                st.success("Loop stopped!")

    # Submit tasks
    if st.session_state.loop_manager is not None:
        duration = st.slider("Task Duration (seconds)", 1, 10, 3)
        if st.button("Submit Task"):
            try:
                task_id = int(time.time())
                future = st.session_state.loop_manager.submit_task(task_id, duration)
                st.info(f"Task {task_id} submitted!")

                # Show task result
                try:
                    result = future.result(timeout=0.1)  # Non-blocking check
                    st.success(f"Task completed: {result}")
                except TimeoutError:
                    st.info("Task still running...")
                except Exception as e:
                    st.error(f"Task error: {e}")

            except Exception as e:
                st.error(f"Error submitting task: {e}")


if __name__ == "__main__":
    main()
