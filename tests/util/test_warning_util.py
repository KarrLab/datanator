""" Tests of the warning utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from capturer import CaptureOutput
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import warning_util
import requests
import unittest


class TestWarningUtil(unittest.TestCase):

    def test_set_warnings_openbabel(self):
        adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'

        with CaptureOutput() as capturer:
            molecule_util.Molecule(adp).to_inchi()
            self.assertNotEqual(capturer.get_text(), '')

        warning_util.set_warnings()

        with CaptureOutput() as capturer:
            molecule_util.Molecule(adp).to_inchi()
            self.assertEqual(capturer.get_text(), '')

    @unittest.skip('todo: implement')
    def test_set_warnings_urllib3(self):
        pass
