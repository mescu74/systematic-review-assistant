"""Database setup.

See `async SQLAlchemy <https://docs.sqlalchemy.org/en/20/orm/extensions/asyncio.html>`_

.. code-block:: python

    # refresh a collection
    # force the collection to load by naming it in attribute_names
    await async_session.refresh(a_obj, ["bs"])
"""

import typing as t

import streamlit as st
from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool
from sqlmodel import create_engine
from sqlmodel.ext.asyncio.session import AsyncSession as AsyncSQLModelSession
from sqlmodel.orm.session import Session as SQLModelSession
from sqlmodel.sql.expression import select

import sr_assistant.app.utils as ut
from sr_assistant.app.config import get_settings
from sr_assistant.core.models import SQLModelBase

settings = get_settings()

# sync
engine = create_engine(
    url=str(settings.DATABASE_URL),
    echo=False,
    pool_size=10,  # Number of connections to maintain in the pool
    max_overflow=20,  # Additional connections beyond pool_size
    pool_timeout=30,  # Seconds to wait for a connection from the pool
    pool_recycle=3600,  # Recycle connections after 1 hour
    pool_pre_ping=True,  # Validate connections before use
)

session_factory = sessionmaker(
    bind=engine, class_=SQLModelSession, expire_on_commit=False
)
"""`sessionmaker` context manager for sync sessions.

Usage: ``with session_factory.begin() as session:`` to auto-commit and rollback
on exit.
"""

if ut.in_streamlit() and "session_factory" not in st.session_state:
    st.session_state.session_factory = session_factory


def create_tables() -> None:
    with engine.begin() as conn:
        conn.run_sync(SQLModelBase.metadata.create_all)


## ---------------------------------------------------- queries


def persist_model(model: SQLModelBase) -> SQLModelBase:
    session_factory = t.cast(
        "sessionmaker[SQLModelSession]", st.session_state.session_factory
    )
    # with session_factory() as session:
    #    model = session.merge(model)
    #    session.commit()
    #    session.refresh(model)
    #    return model
    with session_factory.begin() as session:
        return session.merge(model)


def is_persisted(model: SQLModelBase) -> bool:
    with st.session_state.session_factory.begin() as session:
        stmt = select(type(model)).where(model.id == type(model).id)
        response = session.exec(stmt).first()
        return isinstance(response, type(model))


## ---------------------------------------------------- queries

# async
async_engine = create_async_engine(
    url=str(settings.DATABASE_URL),
    echo=False,
    pool_size=10,
    max_overflow=20,
    pool_timeout=30,
    pool_recycle=3600,
    pool_pre_ping=True,
)

asession_factory = async_sessionmaker(
    bind=async_engine,
    expire_on_commit=False,
    class_=AsyncSQLModelSession,
    sync_session_class=session_factory,
    # https://docs.sqlalchemy.org/en/20/orm/extensions/asyncio.html#using-multiple-asyncio-event-loops
    poolclass=NullPool,
)
"""`async_sessionmaker` context manager for async sessions.

Usage: ``async with asession_factory.begin() as session:`` to auto-commit and rollback
on exit.
"""

if ut.in_streamlit() and "asession_factory" not in st.session_state:
    st.session_state.asession_factory = asession_factory
"""Session state is thread-safe (RLock) and can be shared across pools if you've the ctx holding it."""


async def acreate_tables() -> None:
    async with async_engine.begin() as conn:
        await conn.run_sync(SQLModelBase.metadata.create_all)
