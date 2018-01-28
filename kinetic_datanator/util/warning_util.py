""" Warning utilities

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

import openbabel
import requests.packages.urllib3


def disable_warnings():
    """ Disable warning messages from openbabel and urllib """
    openbabel.obErrorLog.SetOutputLevel(openbabel.obError)
    requests.packages.urllib3.disable_warnings(requests.packages.urllib3.exceptions.InsecureRequestWarning)


def enable_warnings():
    """ Enable warning messages from openbabel and urllib """
    openbabel.obErrorLog.SetOutputLevel(openbabel.obWarning)
    requests.packages.urllib3.warnings.resetwarnings()
