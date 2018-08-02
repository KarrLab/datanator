import os

from flask import Flask, render_template
from flask_login import LoginManager
from flask_bcrypt import Bcrypt
from flask_debugtoolbar import DebugToolbarExtension
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate

import re


# instantiate the extensions
login_manager = LoginManager()
bcrypt = Bcrypt()
toolbar = DebugToolbarExtension()
migrate = Migrate()
db = SQLAlchemy()


def create_app():

    #TODO: Include API Templates

    app = Flask(__name__)

    # set config
    app_settings = os.getenv(
        'APP_SETTINGS', 'kinetic_datanator.config.config.CircleTestingConfig')
    app.config.from_object(app_settings)

    # set up extensions
    login_manager.init_app(app)
    bcrypt.init_app(app)
    toolbar.init_app(app)
    db.init_app(app)
    migrate.init_app(app, db)

    # register blueprints
    from kinetic_datanator.api.urls import api_blueprint
    app.register_blueprint(api_blueprint)


    # # flask login
    # from app.server.model import User
    # login_manager.login_view = 'user.login'
    # login_manager.login_message_category = 'danger'
    #
    # @login_manager.user_loader
    # def load_user(user_id):
    #     return User.query.filter(User.id == int(user_id)).first()

    # @app.errorhandler(401)
    # def unauthorized_page(error):
    #     return render_template('errors/401.html'), 401
    #
    # @app.errorhandler(403)
    # def forbidden_page(error):
    #     return render_template('errors/403.html'), 403
    #
    # @app.errorhandler(404)
    # def page_not_found(error):
    #     return render_template('errors/404.html'), 404
    #
    # @app.errorhandler(500)
    # def server_error_page(error):
    #     return render_template('errors/500.html'), 500

    return app
