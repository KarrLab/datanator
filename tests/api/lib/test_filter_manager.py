""" Test of text metabolite manager

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-08
:Copyright: 2018, Karr Lab
:License: MIT
"""

from kinetic_datanator.api.lib.filter.manager import FilterManager
from kinetic_datanator.api.lib.metabolite.manager import metabolite_manager
from kinetic_datanator.core import common_schema, models
import unittest
import tempfile
import shutil

class TestFilterManager(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.adenine = metabolite_manager._search('adenine')

    def test_run(self):
        results = FilterManager(data=self.adenine)
