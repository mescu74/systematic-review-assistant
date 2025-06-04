"""Pytest fixtures for integration testing."""

import os
from collections.abc import Generator
from typing import Any

import pytest
import sqlalchemy as sa
from dotenv import load_dotenv
from sqlmodel import Session, create_engine


# FIXME: THIS EXPOSES CREDENTIALS IN LOGS/CI! See alembic/env.py for how to mask the conn string.
@pytest.fixture(scope="session")
def db_engine() -> Generator[sa.Engine, None, None]:
    """Yields a SQLAlchemy engine scoped to the test session."""
    load_dotenv(".env.test", override=True)
    # WARN: Contains plaintext credentials!
    db_url = os.getenv("SRA_DATABASE_URL", os.getenv("DATABASE_URL", ""))

    if not db_url:
        pytest.fail("DATABASE_URL or SRA_DATABASE_URL not set in environment/.env.test")

    if "sra_integration_test" not in db_url:
        pytest.fail(
            f"DATABASE_URL '{db_url}' does not seem to be the integration test DB ('sra_integration_test')"
        )

    print(f"\nCreating session engine for: {db_url}")
    engine = create_engine(db_url, echo=False)
    yield engine
    print("\nDisposing session engine.")
    engine.dispose()


@pytest.fixture(scope="function", autouse=True)
def clean_db(db_engine: sa.Engine, request: pytest.FixtureRequest):
    """Cleans (drops/recreates public schema) and recreates tables before each integration test."""
    marker = request.node.get_closest_marker("integration")
    if marker is None:
        yield
        return

    print("\nRunning DB cleanup for integration test...")
    # Drop and recreate public schema directly
    with db_engine.connect() as conn:
        conn.execute(sa.text("DROP SCHEMA public CASCADE;"))
        conn.execute(sa.text("CREATE SCHEMA public;"))
        # Grant permissions if necessary for your test user/role
        conn.execute(
            sa.text("GRANT ALL ON SCHEMA public TO postgres;")
        )  # Usually the owner in Supabase
        conn.execute(
            sa.text("GRANT ALL ON SCHEMA public TO public;")
        )  # Or specific test roles
        conn.commit()
    print("Public schema dropped and recreated.")

    # Recreate schema using Alembic
    # print("Applying Alembic migrations to test database...")
    # alembic_env = os.environ.copy()
    # alembic_env["ENVIRONMENT"] = "test"
    # # Ensure SRA_DATABASE_URL is set for alembic in this specific env context if not already inherited perfectly
    # # db_url is already fetched and validated above, so we can reuse it.
    # alembic_env["SRA_DATABASE_URL"] = db_engine.url.render_as_string(
    #     hide_password=False
    # )

    # alembic_command = ["uv", "run", "alembic", "upgrade", "head"]
    # # Using subprocess.run to execute alembic
    # # Need to capture output to see if alembic itself had issues.
    # result = subprocess.run(
    #     alembic_command, env=alembic_env, capture_output=True, text=True, check=False
    # )
    # if result.returncode != 0:
    #     print("Alembic upgrade stdout:", result.stdout)
    #     print("Alembic upgrade stderr:", result.stderr)
    #     pytest.fail(f"Alembic upgrade head failed for test DB setup: {result.stderr}")
    # else:
    #     print("Alembic upgrade head successful for test DB.")
    #     if result.stdout:
    #         print("Alembic stdout:", result.stdout)
    #     if result.stderr:  # Log stderr even on success, might contain warnings
    #         print("Alembic stderr:", result.stderr)

    # # SQLModel.metadata.create_all(db_engine) # No longer using this for schema creation
    # print("Test database schema setup via Alembic complete.")

    # Create schema using SQLModel definitions
    print("Creating database schema using SQLModel.metadata.create_all()...")
    from sr_assistant.core import (
        models,  # Import models to ensure metadata is populated
    )

    models.SQLModel.metadata.create_all(db_engine)
    print("Test database schema setup via SQLModel.metadata.create_all() complete.")

    # Attempt to reset connection state
    # with db_engine.connect().execution_options(isolation_level="AUTOCOMMIT") as conn:
    #     conn.execute(sa.text("DISCARD ALL;"))
    #     # No conn.commit() needed due to autocommit
    # print("PostgreSQL session state discarded using DISCARD ALL.")

    yield  # Test runs here
    print("DB cleanup fixture finished.")


@pytest.fixture
def db_session(db_engine: sa.Engine, clean_db: Any) -> Generator[Session, None, None]:
    """Yields a DB session with transaction management for a test function."""
    with Session(db_engine) as session:
        print("\nYielding DB session for test...")
        try:
            yield session
        finally:
            print("Closing DB session after test.")
