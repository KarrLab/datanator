import os
basedir = os.path.abspath(os.path.dirname(__file__))


class BaseConfig(object):
    """Base configuration."""
    SECRET_KEY = 'my_precious'
    DEBUG = False
    BCRYPT_LOG_ROUNDS = 13
    WTF_CSRF_ENABLED = False
    DEBUG_TB_ENABLED = False
    DEBUG_TB_INTERCEPT_REDIRECTS = False
    SQLALCHEMY_TRACK_MODIFICATIONS = False


class LocalDevelopmentConfig(BaseConfig):
    """Development configuration."""
    DEBUG = True
    TESTING = False
    BCRYPT_LOG_ROUNDS = 4
    WTF_CSRF_ENABLED = False
    SQLALCHEMY_DATABASE_URI = 'postgres://localhost/CommonSchema'
    # SQLALCHEMY_BINDS = {'data': 'postgres://localhost/User'}
    DEBUG_TB_ENABLED = True



class CircleTestingConfig(BaseConfig):
    """Testing configuration."""
    DEBUG = True
    TESTING = True
    BCRYPT_LOG_ROUNDS = 4
    WTF_CSRF_ENABLED = False
    SQLALCHEMY_DATABASE_URI = 'postgresql://ubuntu@localhost/CommonSchema?sslmode=disable'
    DEBUG_TB_ENABLED = False
    PRESERVE_CONTEXT_ON_EXCEPTION = False


class BuildConfig(BaseConfig):
    """Testing configuration."""
    DEBUG = True
    TESTING = True
    BCRYPT_LOG_ROUNDS = 4
    WTF_CSRF_ENABLED = False
    #TODO: Need to determine the URI
    SQLALCHEMY_DATABASE_URI = 'postgres:///localhost/CommonSchema'
    SQLALCHEMY_BINDS = {'data': 'postgres:///localhost/User'}
    DEBUG_TB_ENABLED = False
    PRESERVE_CONTEXT_ON_EXCEPTION = False

class ProductionConfig(BaseConfig):
    """Production configuration."""
    SECRET_KEY = 'my_precious'
    DEBUG = False
    SQLALCHEMY_DATABASE_URI = 'postgresql://localhost/example'
    SQLALCHEMY_BINDS = {
                        'data': 'sqlite:///{0}'.format(os.path.join(basedir,'cache', 'data.db'))
                        }
    DEBUG_TB_ENABLED = False
    WTF_CSRF_ENABLED = True
    TESTING = False
