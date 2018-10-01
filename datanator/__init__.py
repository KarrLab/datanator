import os

from flask import Flask, render_template
from flask_login import LoginManager
from flask_bcrypt import Bcrypt
from flask_debugtoolbar import DebugToolbarExtension
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask_marshmallow import Marshmallow
from flask_cors import CORS

def register_blueprints(app):
    # register blueprints
    from datanator.api.urls import api_blueprint
    app.register_blueprint(api_blueprint)

#TODO: Include API Templates
app = Flask(
    __name__,
    template_folder='api/client/templates',
    static_folder='api/client/static'
)

# set config
app_settings = os.getenv(
    'APP_SETTINGS', 'datanator.config.config.CircleTestingConfig')

app.config.from_object(app_settings)

login_manager = LoginManager(app)
bcrypt = Bcrypt(app)
toolbar = DebugToolbarExtension(app)

db = SQLAlchemy(app)
ma = Marshmallow(app)
cors = CORS(app)
migrate = Migrate(app, db)
register_blueprints(app)

# # flask login
# from app.server.model import User
# login_manager.login_view = 'user.login'
# login_manager.login_message_category = 'danger'
#
# @login_manager.user_loader
# def load_user(user_id):
#     return User.query.filter(User.id == int(user_id)).first()

@app.errorhandler(401)
def unauthorized_page(error):
    return render_template('errors/401.html'), 401
#
# @app.errorhandler(403)
# def forbidden_page(error):
#     return render_template('errors/403.html'), 403
#
@app.errorhandler(404)
def page_not_found(error):
    return render_template('errors/404.html'), 404
#
# @app.errorhandler(500)
# def server_error_page(error):
#     return render_template('errors/500.html'), 500






import pkg_resources

with open(pkg_resources.resource_filename('datanator', 'VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from . import config
from . import core
from . import datanator
from . import data_source
from . import io
from . import util
