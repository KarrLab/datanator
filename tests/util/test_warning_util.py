""" Tests of the warning utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from capturer import CaptureOutput
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import warning_util
import time
import unittest


class TestWarningUtil(unittest.TestCase):
    adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'

    def test_enable_warnings_openbabel(self):
        warning_util.enable_warnings()
        with CaptureOutput() as capturer:
            molecule_util.Molecule(self.adp).to_inchi()
            time.sleep(0.1)
            self.assertNotEqual(capturer.get_text(), '')

    def test_disable_warnings_openbabel(self):
        warning_util.disable_warnings()
        with CaptureOutput() as capturer:
            molecule_util.Molecule(self.adp).to_inchi()
            time.sleep(0.1)
            self.assertEqual(capturer.get_text(), '')

    @unittest.skip('todo: implement')
    def test_disable_warnings_urllib3(self):
        pass
