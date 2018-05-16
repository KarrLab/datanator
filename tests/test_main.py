""" Tests of the command line utility

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from capturer import CaptureOutput
from cement.utils import test
from kinetic_datanator.__main__ import App
from kinetic_datanator.util import warning_util
import kinetic_datanator
import os
import re
import shutil
import tempfile
import unittest
from os import path

warning_util.disable_warnings()


class BaseControllerTestCase(unittest.TestCase):
    def test_get_version(self):
        with CaptureOutput() as capture_output:
            with App(argv=['-v']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(capture_output.get_text(), kinetic_datanator.__version__)

        with CaptureOutput() as capture_output:
            with App(argv=['--version']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(capture_output.get_text(), kinetic_datanator.__version__)


class TestUploadData(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.dirname = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.dirname)

    def test_upload_ref_seq(self):
        with App(argv=['upload',
                       'reference-genome', "{}/data_source/test_mpn_sequence.gb".format(path.dirname(__file__)),
                       '--path_to_database='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Refseq.sqlite'))


class TestBuildController(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.dirname = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.dirname)

    def test_build_corum(self):
        with App(argv=['build', 'corum', '--path='+self.dirname, '--max-entries=1', '--verbose=True']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Corum.sqlite'))

    def test_build_intact(self):
        with App(argv=['build', 'intact', '--path='+self.dirname, '--max-entries=1', '--verbose=True']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/IntAct.sqlite'))

    @unittest.skip('skip')
    def test_build_sabio(self):
        with App(argv=['build', 'sabio', '--path='+self.dirname, '--max-entries=1', '--verbose=True']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/SabioRk.sqlite'))

    def test_build_jaspar(self):
        with App(argv=['build', 'jaspar', '--path='+self.dirname, '--max-entries=1', '--verbose=True']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Jaspar.sqlite'))

    def test_build_ecmdb(self):
        with App(argv=['build', 'ecmdb', '--path='+self.dirname, '--max-entries=1', '--verbose=True']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Ecmdb.sqlite'))

    def test_build_pax(self):
        with App(argv=['build', 'pax', '--path='+self.dirname, '--max-entries=1', '--verbose=True']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Pax.sqlite'))


class TestDownloadController(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.dirname = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.dirname)

    def test_download_corum(self):
        with App(argv=['download', 'corum', '--path='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Corum.sqlite'))

    def test_download_intact(self):
        with App(argv=['download', 'intact', '--path='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/IntAct.sqlite'))

    def test_download_sabio(self):
        with App(argv=['download', 'sabio', '--path='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/SabioRk.sqlite'))

    def test_download_jaspar(self):
        with App(argv=['download', 'jaspar', '--path='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Jaspar.sqlite'))

    def test_download_ecmdb(self):
        with App(argv=['download', 'ecmdb', '--path='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Ecmdb.sqlite'))

    def test_download_pax(self):
        with App(argv=['download', 'pax', '--path='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Pax.sqlite'))

    def test_download_aggregate(self):
        with App(argv=['download', 'aggregate', '--path='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/FlaskCommonSchema.sqlite'))

    def test_download_array_express(self):
        with App(argv=['download', 'array-express', '--path='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/ArrayExpress.sqlite'))

    def test_download_uniprot(self):
        with App(argv=['download', 'uniprot', '--path='+self.dirname]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Uniprot.sqlite'))


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

        # todo: add assertions after implementing CLI method
        # self.assertTrue(os.path.isfile(output_filename))

    def test_GenerateTemplateController(self):
        filename = os.path.join(self.dirname, 'template.xlsx')
        with App(argv=['generate-template', filename]) as app:
            app.run()
        self.assertTrue(os.path.isfile(filename))


class TestWithoutTempFile(unittest.TestCase):

    def test_taxonomy_get_rank(self):
        with App(argv=['taxonomy', 'get-rank', '2097']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), "species")

        with App(argv=['taxonomy', 'get-rank', 'Mycoplasma genitalium']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), "species")

        with App(argv=['taxonomy', 'get-rank', 'Mycoplasma genitalium XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_parents(self):
        with App(argv=['taxonomy', 'get-parents', 'bacteria']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), "root\ncellular organisms")

        with App(argv=['taxonomy', 'get-parents', 'XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_common_ancestor(self):
        with App(argv=['taxonomy', 'get-common-ancestor', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), 'Mycoplasma')

        with App(argv=['taxonomy', 'get-common-ancestor', 'Mycoplasma genitalium', 'XXX']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_distance_to_common_ancestor(self):
        with App(argv=['taxonomy', 'get-distance-to-common-ancestor', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(float(capturer.get_text()), 1.)

        with App(argv=['taxonomy', 'get-distance-to-common-ancestor', 'Mycoplasma genitalium', 'XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_taxonomy_get_distance_to_root(self):
        with App(argv=['taxonomy', 'get-distance-to-root', 'bacteria']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(float(capturer.get_text()), 2.)

        with App(argv=['taxonomy', 'get-distance-to-root', 'XXX']) as app:
            self.assertRaises(ValueError, lambda: app.run())

    def test_molecule_get_structure(self):
        name = 'water'
        with App(argv=['molecule', 'get-structure', '--by-name', name]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(re.split('  +', capturer.get_text().split('\n')
                                          [-1]), ['water', 'pubchem.compound', '962', 'InChI=1S/H2O/h1H2'])

        namespace = 'chebi'
        id = '15377'
        with App(argv=['molecule', 'get-structure', '--by-id', '--namespace', namespace, id]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
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
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), 'O')

    def test_get_ec_number(self):
        dr1p = 'OCC1OC(CC1O)OP([O-])([O-])=O'
        dr5p = 'OC1CC(O)C(COP([O-])([O-])=O)O1'

        reaction = 'Deoxyribose 1-phosphate --> Deoxyribose-5-P'
        with App(argv=['reaction', 'get-ec-number', reaction]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(re.split('  +', capturer.get_text().split('\n')[2]), ['5.4.2', '16.00'])

        reaction = '{} --> {}'.format(dr1p, dr5p)
        with App(argv=['reaction', 'get-ec-number', reaction]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(re.split('  +', capturer.get_text().split('\n')[2]), ['5.4.2', '16.00'])

        reaction = '{} > {}'.format(dr1p, dr5p)
        with App(argv=['reaction', 'get-ec-number', reaction]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), 'The reaction is ill-formed')

        reaction = 'xxxxxxxxx --> Deoxyribose-5-P'
        with App(argv=['reaction', 'get-ec-number', reaction]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(capturer.get_text().startswith('Unable to interpret participants:\n'))


class HelpTestCase(unittest.TestCase):
    def test(self):
        with App(argv=[]) as app:
            app.run()
        with App(argv=['upload']) as app:
            app.run()
        with App(argv=['build']) as app:
            app.run()
        with App(argv=['download']) as app:
            app.run()
        with App(argv=['taxonomy']) as app:
            app.run()
        with App(argv=['molecule']) as app:
            app.run()
        with App(argv=['reaction']) as app:
            app.run()
