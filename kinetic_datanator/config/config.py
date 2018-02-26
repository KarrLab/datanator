
class BaseConfig(object):
    """Base configuration."""
    SECRET_KEY = 'my_precious'
    DEBUG = False
    DEBUG_TB_ENABLED = False
    DEBUG_TB_INTERCEPT_REDIRECTS = False
    SQLALCHEMY_TRACK_MODIFICATIONS = False


class Config(BaseConfig):
    """Development configuration."""
    DEBUG = True
    SQLALCHEMY_DATABASE_URI = 'sqlite:///:memory:'
    ## The below URI is meant for migrations
    #SQLALCHEMY_DATABASE_URI = 'sqlite:////Users/pochis01/Desktop/GitHub/kinetic_datanator/kinetic_datanator/data_source/cache/FlaskCommonSchema.sqlite'
    DEBUG_TB_ENABLED = True
    SQLALCHEMY_TRACK_MODIFICATIONS = True
    TESTING = True
