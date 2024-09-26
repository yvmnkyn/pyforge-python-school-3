"""Initial migration

Revision ID: 48a578e55f64
Revises: c1836ed8b472
Create Date: 2024-08-30 15:11:29.387518

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = '48a578e55f64'
down_revision: Union[str, None] = 'c1836ed8b472'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    pass


def downgrade() -> None:
    pass
