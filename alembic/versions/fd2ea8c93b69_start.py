"""start

Revision ID: fd2ea8c93b69
Revises:
Create Date: 2018-03-15 12:00:31.726213

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = 'fd2ea8c93b69'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.add_column('organism', sa.Column('ncbi_id', sa.Integer()))


def downgrade():
    pass
