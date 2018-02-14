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
        default=pkg_resources.resource_filename('kinetic_datanator', 'config/default.cfg'),
        schema=pkg_resources.resource_filename('kinetic_datanator', 'config/schema.cfg'),
        user=(
            'kinetic_datanator.cfg',
            os.path.expanduser('~/.wc/kinetic_datanator.cfg'),
        ),
    )

    return wc_utils.config.core.ConfigManager(paths).get_config(extra=extra)


def get_debug_logs_config(extra=None):
    """ Get debug logs configuration

    Args:
        extra (:obj:`dict`, optional): additional configuration to override

    Returns:
        :obj:`configobj.ConfigObj`: nested dictionary with the configuration settings loaded from the configuration source(s).
    """
    paths = wc_utils.debug_logs.config.paths.deepcopy()
    paths.user = (
        'kinetic_datanator.debug.cfg',
        os.path.expanduser('~/.wc/kinetic_datanator.debug.cfg'),
    )
    return wc_utils.config.core.ConfigManager(paths).get_config(extra=extra)
