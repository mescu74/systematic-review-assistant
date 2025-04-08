"""tz aware timestamps

Revision ID: 0ab256b04824
Revises: edd72f8c1a16
Create Date: 2025-01-31 07:11:51.408417+00:00

"""

from typing import Union
from collections.abc import Sequence

from alembic import op
import sqlalchemy as sa
import sqlmodel


# revision identifiers, used by Alembic.
revision: str = "0ab256b04824"
down_revision: str | None = "edd72f8c1a16"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    # Abstract screening results
    for column in ("created_at", "updated_at"):
        op.alter_column(
            table_name="abstract_screening_results",
            column_name=column,
            nullable=True,
            existing_type=sa.DateTime(),
            type_=sa.DateTime(timezone=True),
            postgresql_using=f"{column} AT TIME ZONE 'UTC'",
        )

    # Reviews
    for column in ("created_at", "updated_at"):
        op.alter_column(
            table_name="reviews",
            column_name=column,
            nullable=True,
            existing_type=sa.DateTime(),
            type_=sa.DateTime(timezone=True),
            postgresql_using=f"{column} AT TIME ZONE 'UTC'",
        )

    # Pubmed results
    op.alter_column(
        table_name="pubmed_results",
        column_name="created_at",
        nullable=True,
        existing_type=sa.DateTime(),
        type_=sa.DateTime(timezone=True),
        postgresql_using="created_at AT TIME ZONE 'UTC'",
    )


def downgrade() -> None:
    # Abstract screening results
    for column in ("created_at", "updated_at"):
        op.alter_column(
            table_name="abstract_screening_results",
            column_name=column,
            nullable=True,
            existing_type=sa.DateTime(timezone=True),  # Note the timezone=True here
            type_=sa.DateTime(),
            postgresql_using=f"{column} AT TIME ZONE 'UTC'",
        )

    # Reviews
    for column in ("created_at", "updated_at"):
        op.alter_column(
            table_name="reviews",
            column_name=column,
            nullable=True,
            existing_type=sa.DateTime(timezone=True),  # Note the timezone=True here
            type_=sa.DateTime(),
            postgresql_using=f"{column} AT TIME ZONE 'UTC'",
        )

    # Pubmed results
    op.alter_column(
        table_name="pubmed_results",
        column_name="created_at",
        nullable=True,
        existing_type=sa.DateTime(timezone=True),  # Note the timezone=True here
        type_=sa.DateTime(),
        postgresql_using="created_at AT TIME ZONE 'UTC'",
    )
