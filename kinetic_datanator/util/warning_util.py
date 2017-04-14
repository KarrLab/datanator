""" Warning utilities

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

from requests.packages import urllib3
import openbabel


def set_warnings():
    """ Suppress warning messages from openbabel and urllib"""
    openbabel.obErrorLog.SetOutputLevel(openbabel.obError)
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
