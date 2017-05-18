""" Tests of the command line utility

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from capturer import CaptureOutput
from kinetic_datanator.__main__ import App
from kinetic_datanator.util import warning_util
import os
import re
import shutil
import tempfile
import time
import unittest

warning_util.disable_warnings()


class TestWithTempFile(unittest.TestCase):

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_GetDataController(self):
        input_filename = os.path.join(os.path.dirname(__file__), 'fixtures', 'five_reactions.xlsx')
        output_filename = os.path.join(self.dirname, 'output.xlsx')
        argv = ['get-data', input_filename, output_filename,
                '--max-taxon-dist', '8',
                '--taxon-dist-scale', '1.6',
                '--include-variants',
                '--temperature', '37',
                '--temperature-std', '1',
                '--ph', '7.5',
                '--ph-std', '0.3',
                ]
        with App(argv=argv) as app:
            app.run()

        # todo: add assertions
        # self.assertTrue(os.path.isfile(output_filename))

    def test_GenerateTemplateController(self):
        filename = os.path.join(self.dirname, 'template.xlsx')
        with App(argv=['generate-template', filename]) as app:
            app.run()
        self.assertTrue(os.path.isfile(filename))


class TestWithoutTempFile(unittest.TestCase):

    def test_taxonomy_get_rank(self):
        with App(argv=['taxonomy', 'get-rank', '2097']) as app:
            with CaptureOutput() as capturer:
                app.run()
                time.sleep(0.1)
                self.assertEqual(capturer.get_text(), "species")

        with App(argv=['taxonomy', 'get-rank', 'Mycoplasma genitalium']) as app:
            with CaptureOutput() as capturer:
                app.run()
                time.sleep(0.1)
                self.assertEqual(capturer.get_text(), "species")

        with App(argv=['taxonomy', 'get-rank', 'Mycoplasma genitalium XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_parents(self):
        with App(argv=['taxonomy', 'get-parents', 'bacteria']) as app:
            with CaptureOutput() as capturer:
                app.run()
                time.sleep(0.1)
                self.assertEqual(capturer.get_text(), "root\ncellular organisms")

        with App(argv=['taxonomy', 'get-parents', 'XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_common_ancestor(self):
        with App(argv=['taxonomy', 'get-common-ancestor', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae']) as app:
            with CaptureOutput() as capturer:
                app.run()
                time.sleep(0.1)
                self.assertEqual(capturer.get_text(), 'Mycoplasma')

        with App(argv=['taxonomy', 'get-common-ancestor', 'Mycoplasma genitalium', 'XXX']) as app:
            with CaptureOutput() as capturer:
                self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_distance_to_common_ancestor(self):
        with App(argv=['taxonomy', 'get-distance-to-common-ancestor', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae']) as app:
            with CaptureOutput() as capturer:
                app.run()
                time.sleep(0.1)
                self.assertEqual(float(capturer.get_text()), 1.)

        with App(argv=['taxonomy', 'get-distance-to-common-ancestor', 'Mycoplasma genitalium', 'XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_distance_to_root(self):
        with App(argv=['taxonomy', 'get-distance-to-root', 'bacteria']) as app:
            with CaptureOutput() as capturer:
                app.run()
                time.sleep(0.1)
                self.assertEqual(float(capturer.get_text()), 2.)

        with App(argv=['taxonomy', 'get-distance-to-root', 'XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_molecule_get_structure(self):
        name = 'water'
        with App(argv=['molecule', 'get-structure', '--by-name', name]) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(re.split('  +', capturer.get_text().split('\n')
                                          [-1]), ['water', 'pubchem.compound', '962', 'InChI=1S/H2O/h1H2'])

        namespace = 'chebi'
        id = '15377'
        with App(argv=['molecule', 'get-structure', '--by-id', '--namespace', namespace, id]) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(re.split('  +', capturer.get_text().split('\n')[-1]), ['', 'chebi', '15377', 'InChI=1S/H2O/h1H2'])

        namespace = 'chebi'
        id = '0'
        with App(argv=['molecule', 'get-structure', '--by-id', '--namespace', namespace, id]) as app:
            with CaptureOutput(merged=False) as capturer:
                app.run()
                self.assertEqual(capturer.stderr.get_text(), 'Unable to find structure')

    def test_molecule_convert_structure(self):
        with App(argv=['molecule', 'convert-structure', 'O', 'can']) as app:
            with CaptureOutput() as capturer:
                app.run()
                time.sleep(0.1)
                self.assertEqual(capturer.get_text(), 'O')

    def test_get_ec_number(self):
        dr1p = 'OCC1OC(CC1O)OP([O-])([O-])=O'
        dr5p = 'OC1CC(O)C(COP([O-])([O-])=O)O1'

        reaction = 'Deoxyribose 1-phosphate --> Deoxyribose-5-P'
        with App(argv=['reaction', 'get-ec-number', reaction]) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(re.split('  +', capturer.get_text().split('\n')[2]), ['5.4.2', '16.00'])

        reaction = '{} --> {}'.format(dr1p, dr5p)
        with App(argv=['reaction', 'get-ec-number', reaction]) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(re.split('  +', capturer.get_text().split('\n')[2]), ['5.4.2', '16.00'])

        reaction = '{} > {}'.format(dr1p, dr5p)
        with App(argv=['reaction', 'get-ec-number', reaction]) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), 'The reaction is ill-formed')

        reaction = 'xxxxxxxxx --> Deoxyribose-5-P'
        with App(argv=['reaction', 'get-ec-number', reaction]) as app:
            with CaptureOutput() as capturer:
                app.run()
                self.assertTrue(capturer.get_text().startswith('Unable to interpret participants:\n'))
