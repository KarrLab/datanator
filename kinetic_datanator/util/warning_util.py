""" Warning utilities

:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

from requests.packages import urllib3


def set_warnings():
    requests.packages.urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
