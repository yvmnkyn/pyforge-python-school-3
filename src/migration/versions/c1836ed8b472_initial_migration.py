"""Initial migration

Revision ID: c1836ed8b472
Revises: 19b4077419f8
Create Date: 2024-08-30 14:50:30.895384

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = 'c1836ed8b472'
down_revision: Union[str, None] = '19b4077419f8'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    pass


def downgrade() -> None:
    pass
