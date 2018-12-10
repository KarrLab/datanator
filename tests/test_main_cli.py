""" Tests of the command line utility

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from capturer import CaptureOutput
from cement.utils import test
from datanator.__main__ import App
from datanator.util import warning_util
import datanator
import mock
import os
import re
import shutil
import sqlalchemy.orm
import sqlalchemy_utils
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
                self.assertEqual(capture_output.get_text(), datanator.__version__)

        with CaptureOutput() as capture_output:
            with App(argv=['--version']) as app:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(capture_output.get_text(), datanator.__version__)


class TestUploadData(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dirname = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dirname)

    def test_upload_ref_seq(self):
        with App(argv=['upload', 'reference-genome',
                       os.path.join(path.dirname(__file__), 'data_source', 'test_mpn_sequence.gb'),
                       '--db-path', self.dirname]) as app:
            app.run()
        self.assertTrue(os.path.exists(os.path.join(self.dirname, 'Refseq.sqlite')))


class TestBuildController(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dirname = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dirname)

    def test_build_corum(self):
        with App(argv=['build', '--path='+self.dirname, '--max-entries=1', '--verbose=True', 'corum']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Corum.sqlite'))

    def test_build_intact(self):
        with App(argv=['build', '--path='+self.dirname, '--max-entries=1', '--verbose=True', 'intact']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/IntAct.sqlite'))

    @unittest.skip('skip')
    def test_build_sabio(self):
        with App(argv=['build', '--path='+self.dirname, '--max-entries=1', '--verbose=True', 'sabio']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/SabioRk.sqlite'))

    def test_build_jaspar(self):
        with App(argv=['build', '--path='+self.dirname, '--max-entries=1', '--verbose=True', 'jaspar']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Jaspar.sqlite'))

    def test_build_ecmdb(self):
        with App(argv=['build', '--path='+self.dirname, '--max-entries=1', '--verbose=True', 'ecmdb']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Ecmdb.sqlite'))

    def test_build_pax(self):
        with App(argv=['build', '--path='+self.dirname, '--max-entries=1', '--verbose=True', 'pax']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Pax.sqlite'))


@unittest.skip('skip because too long')
class TestDownloadController(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dirname = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dirname)

    def test_download_corum(self):
        with App(argv=['download', '--path='+self.dirname, 'corum']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Corum.sqlite'))

    def test_download_intact(self):
        with App(argv=['download', '--path='+self.dirname, 'intact']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/IntAct.sqlite'))

    def test_download_sabio(self):
        with App(argv=['download', '--path='+self.dirname, 'sabio']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/SabioRk.sqlite'))

    def test_download_jaspar(self):
        with App(argv=['download', '--path='+self.dirname, 'jaspar']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Jaspar.sqlite'))

    def test_download_ecmdb(self):
        with App(argv=['download', '--path='+self.dirname, 'ecmdb']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Ecmdb.sqlite'))

    def test_download_pax(self):
        with App(argv=['download', '--path='+self.dirname, 'pax']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/Pax.sqlite'))

    def test_download_array_express(self):
        with App(argv=['download', '--path='+self.dirname, 'array-express']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertTrue(os.path.exists(self.dirname+'/ArrayExpress.sqlite'))

    def test_download_uniprot(self):
        with App(argv=['download', '--path='+self.dirname, 'uniprot']) as app:
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
            self.assertRaises(SystemExit, lambda: app.run())

    def test_taxonomy_get_parents(self):
        with App(argv=['taxonomy', 'get-parents', 'bacteria']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), "root\ncellular organisms")

        with App(argv=['taxonomy', 'get-parents', 'XXX']) as app:
            self.assertRaises(SystemExit, lambda: app.run())

    def test_taxonomy_get_common_ancestor(self):
        with App(argv=['taxonomy', 'get-common-ancestor', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(capturer.get_text(), 'Mycoplasma')

        with App(argv=['taxonomy', 'get-common-ancestor', 'Mycoplasma genitalium', 'XXX']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                self.assertRaises(SystemExit, lambda: app.run())

    def test_taxonomy_get_distance_to_common_ancestor(self):
        with App(argv=['taxonomy', 'get-distance-to-common-ancestor', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(float(capturer.get_text()), 1.)

        with App(argv=['taxonomy', 'get-distance-to-common-ancestor', 'Mycoplasma genitalium', 'XXX']) as app:
            self.assertRaises(SystemExit, lambda: app.run())

    def test_taxonomy_get_distance_to_root(self):
        with App(argv=['taxonomy', 'get-distance-to-root', 'bacteria']) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                self.assertEqual(float(capturer.get_text()), 2.)

        with App(argv=['taxonomy', 'get-distance-to-root', 'XXX']) as app:
            self.assertRaises(SystemExit, lambda: app.run())

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

    @unittest.skip('skipping for development purposes')
    def test_get_ec_number(self):
        dr1p = 'OCC1OC(CC1O)OP([O-])([O-])=O'
        dr5p = 'OC1CC(O)C(COP([O-])([O-])=O)O1'

        reaction = 'Deoxyribose 1-phosphate --> Deoxyribose-5-P'
        with App(argv=['reaction', 'get-ec-number', reaction]) as app:
            with CaptureOutput(termination_delay=0.1) as capturer:
                app.run()
                print(re.split('  +', capturer.get_text().split('\n')))
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


class DbControllerTestCase(unittest.TestCase):
    def setUp(self):
        self.url = datanator.db.engine.url
        if sqlalchemy_utils.functions.database_exists(self.url):
            sqlalchemy_utils.functions.drop_database(self.url)
            datanator.db.engine.dispose()
        self.assertFalse(sqlalchemy_utils.functions.database_exists(self.url))

        if os.path.isdir('migrations'):
            shutil.rmtree('migrations')
        self.assertFalse(os.path.isdir('migrations'))

    def tearDown(self):
        if sqlalchemy_utils.functions.database_exists(self.url):
            sqlalchemy_utils.functions.drop_database(self.url)
            datanator.db.engine.dispose()
        self.assertFalse(sqlalchemy_utils.functions.database_exists(self.url))

        if os.path.isdir('migrations'):
            shutil.rmtree('migrations')
        self.assertFalse(os.path.isdir('migrations'))

    def test_create(self):
        with App(argv=['db', 'create']) as app:
            app.run()
        self.assertTrue(sqlalchemy_utils.functions.database_exists(self.url))
        self.assertNotEqual(datanator.db.engine.table_names(), [])

    def test_create_table_only(self):
        sqlalchemy_utils.functions.create_database(self.url)
        self.assertTrue(sqlalchemy_utils.functions.database_exists(self.url))
        self.assertEqual(datanator.db.engine.table_names(), [])
        with App(argv=['db', 'create']) as app:
            app.run()
        self.assertNotEqual(datanator.db.engine.table_names(), [])

    def test_migrate(self):
        with App(argv=['db', 'create']) as app:
            app.run()

        with App(argv=['db', 'migrate']) as app:
            app.run()
        self.assertTrue(os.path.isdir('migrations'))
        self.assertNotEqual(os.listdir('migrations'), [])

    def test_drop(self):
        with App(argv=['db', 'create']) as app:
            app.run()
        self.assertTrue(sqlalchemy_utils.functions.database_exists(self.url))

        with App(argv=['db', 'drop']) as app:
            app.run()
        self.assertFalse(sqlalchemy_utils.functions.database_exists(self.url))

    def test_restore(self):
        with App(argv=['db', 'create']) as app:
            app.run()
        self.assertTrue(sqlalchemy_utils.functions.database_exists(self.url))

        session = sqlalchemy.orm.sessionmaker(bind=datanator.db.engine)()
        query = session.query(datanator.core.models.Observation)
        self.assertEqual(query.count(), 0)
        session.close()

        with App(argv=['db', 'restore', '--restore-schema', '--do-not-exit-on-error']) as app:
            # todo: remove --restore-schema and --do-not-exit-on-error after fixing Alembic issue with migrations
            app.run()

        session = sqlalchemy.orm.sessionmaker(bind=datanator.db.engine)()
        query = session.query(datanator.core.models.Observation)
        self.assertGreater(query.count(), 0)
        session.close()
