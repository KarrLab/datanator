from flask_script import Manager
from flask_migrate import Migrate, MigrateCommand

from kinetic_datanator import create_app, db

app = create_app()
manager = Manager(app)

@manager.command
def create_db():
    """Creates the db tables."""
    db.create_all()

@manager.command
def drop_db():
    """Drops the db tables."""
    db.drop_all()

# @manager.command
# def create_admin():
#     """Creates the admin user."""
#     db.session.add(User(username='admin', email='ad@min.com', password='admin', admin=True))
#     db.session.commit()


if __name__ == '__main__':
    manager.run()
