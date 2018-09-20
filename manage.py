from flask_script import Manager, Server
from flask_migrate import Migrate, MigrateCommand
from kinetic_datanator import app, db
from kinetic_datanator.core import common_schema
# from app.server.model import User

migrate = Migrate(app, db)
manager = Manager(app)

# migrations
manager.add_command('migrate', MigrateCommand)
manager.add_command("runserver", Server(host='localhost'))


@manager.command
def create_db():
    """Create the database tables."""
    db.create_all()


@manager.command
def drop_db():
    """Drop the database tables."""
    db.drop_all()


@manager.option('--do-not-restore-data', dest='restore_data', action='store_false',
                help='If set, do not restore the data')
@manager.option('--restore-schema', dest='restore_schema', action='store_true',
                help='If set, restore the schema')
@manager.option('--do-not-exit-on-error', dest='exit_on_error', action='store_false',
                help='If set, do not exit on errors')
def restore_db(restore_data=True, restore_schema=False, exit_on_error=True):
    """ Restore the content of the database 
    Args:
        restore_data (:obj:`bool`, optional): if :obj:`True`, restore the data
        restore_schema (:obj:`bool`, optional): if :obj:`True`, restore the schema
        exit_on_error (:obj:`bool`, optional): if :obj:`True`, exit on error
    """
    if not restore_data and not restore_schema:
        raise Exception('One or more of the data and schema must be restored')

    common_schema.CommonSchema(clear_content=True,
                               restore_backup_data=restore_data,
                               restore_backup_schema=restore_schema,
                               restore_backup_exit_on_error=exit_on_error,
                               load_content=False,
                               verbose=True)


# @manager.command
# def create_admin():
#     """Creates the admin user."""
#     db.session.add(User(username='admin', email='ad@min.com', password='admin', admin=True))
#     db.session.commit()


if __name__ == '__main__':
    manager.run()
