""" Configuration

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

import configobj
import os
import pkg_resources
import wc_utils.config.core
import wc_utils.debug_logs.config


def get_config(extra=None):
    """ Get configuration

    Args:
        extra (:obj:`dict`, optional): additional configuration to override

    Returns:
        :obj:`configobj.ConfigObj`: nested dictionary with the configuration settings loaded from the configuration source(s).
    """
    paths = wc_utils.config.core.ConfigPaths(
        default=pkg_resources.resource_filename('datanator', 'config/core.default.cfg'),
        schema=pkg_resources.resource_filename('datanator', 'config/core.schema.cfg'),
        user=(
            'datanator.ini',
            os.path.expanduser('~/.wc/datanator.cfg'),
        )
    )

    return wc_utils.config.core.ConfigManager(paths).get_config(extra=extra)
def get_mongo_config():
    """ Get a configuration to pass directly into the mongo client
    Args:
        extra (:obj: 'dict', optional): override the Mongo information loaded from the config file
    Returns: 
        :obj:'dict': dictionary containing parameters to pass into the MongoDB util constructor
    """
    config=get_config()
    username =config['datanator']['mongodb']['user']
    password =config['datanator']['mongodb']['password']
    port = config['datanator']['mongodb']['port']
    MongoDB = config['datanator']['mongodb']['server']
    replSet = config['datanator']['mongodb']['replSet']
    mongo_config = {"MongoDB":MongoDB,"username":username, "password": password, "replicaSet": replSet}
    return mongo_config




def get_debug_logs_config(extra=None):
    """ Get debug logs configuration

    Args:
        extra (:obj:`dict`, optional): additional configuration to override

    Returns:
        :obj:`configobj.ConfigObj`: nested dictionary with the configuration settings loaded from the configuration source(s).
    """
    paths = wc_utils.debug_logs.config.paths.deepcopy()
    paths.user = (
        'datanator.debug.cfg',
        os.path.expanduser('~/.wc/datanator.debug.cfg'),
    )
    return wc_utils.config.core.ConfigManager(paths).get_config(extra=extra)
