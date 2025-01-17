import os
import re
from logging.config import fileConfig
from typing import NoReturn

from alembic import context
from sqlalchemy import engine_from_config, pool
from sqlmodel import SQLModel


def die(msg: str) -> NoReturn:
    msg = f"ERROR: {msg}"
    raise SystemExit(msg)


# Safety checks
env = os.environ.get("ENVIRONMENT")
# note that prototype code uses "prod" for hosted Supabase, while GHA has a "prototype"
# environment.
valid_envs = ["CI", "local", "dev", "staging", "prod"]
direct_url = os.environ.get("SUPABASE_SR_DIRECT_URL")

if not direct_url:
    die("SUPABASE_SR_DIRECT_URL environment variable is required")

if env not in valid_envs:
    die(f"Unknown ENVIRONMENT: got {env!r}, want {valid_envs!r}")

if env != "prod":
    # Matches Supabase production database URL pattern
    prod_pattern = r"db\.[a-z0-9-]+\.supabase\.co"
    if re.search(prod_pattern, direct_url):
        die(f"Production database URL detected but ENVIRONMENT is {env!r}")

# this is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config

config.set_main_option("sqlalchemy.url", direct_url)

# Interpret the config file for Python logging.
# This line sets up loggers basically.
if config.config_file_name is not None:
    fileConfig(config.config_file_name)


# add your model's MetaData object here
# for 'autogenerate' support
# from myapp import mymodel
# target_metadata = mymodel.Base.metadata
# Note: SQlModel.metadata has access to all tables inheriting from SQLModel
target_metadata = SQLModel.metadata

# other values from the config, defined by the needs of env.py,
# can be acquired:
# my_important_option = config.get_main_option("my_important_option")
# ... etc.


def run_migrations_offline() -> None:
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.

    """
    url = config.get_main_option("sqlalchemy.url")
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """Run migrations in 'online' mode.

    In this scenario we need to create an Engine
    and associate a connection with the context.

    """
    connectable = engine_from_config(
        config.get_section(config.config_ini_section, {}),
        prefix="sqlalchemy.",
        poolclass=pool.NullPool,
    )

    with connectable.connect() as connection:
        context.configure(connection=connection, target_metadata=target_metadata)

        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
