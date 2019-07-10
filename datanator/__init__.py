import os

import pkg_resources

with open(pkg_resources.resource_filename('VERSION'), 'r') as file:
    __version__ = file.read().strip()
# :obj:`str`: version

# API
from . import config
from . import core
from . import data_source
from . import util