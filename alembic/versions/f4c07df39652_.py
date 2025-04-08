"""empty message

Revision ID: f4c07df39652
Revises: 2c7373437db5
Create Date: 2025-04-08 02:43:54.247798+00:00

"""
from typing import Union
from collections.abc import Sequence

from alembic import op
import sqlalchemy as sa
import sqlmodel


# revision identifiers, used by Alembic.
revision: str = 'f4c07df39652'
down_revision: str | None = ('2c7373437db5')
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    pass


def downgrade() -> None:
    pass
