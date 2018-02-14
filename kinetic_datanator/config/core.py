""" Configuration

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

from pkg_resources import resource_filename
from wc_utils.config.core import ConfigPaths
from wc_utils.debug_logs.config import paths as debug_logs_default_paths
import os


core = ConfigPaths(
    default=resource_filename('kinetic_datanator', 'config/default.cfg'),
    schema=resource_filename('kinetic_datanator', 'config/schema.cfg'),
    user=(
        'kinetic_datanator.cfg',
        os.path.expanduser('~/.wc/kinetic_datanator.cfg'),
    ),
)

debug_logs = debug_logs_default_paths.deepcopy()
debug_logs.user = (
    'kinetic_datanator.debug.cfg',
    os.path.expanduser('~/.wc/kinetic_datanator.debug.cfg'),
)
