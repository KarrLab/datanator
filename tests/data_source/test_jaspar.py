"""
This code tests all aspects of jaspar.py

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-24
:Copyright: 2017, Karr Lab
:License: MIT

"""

import unittest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from kinetic_datanator.data_source import jaspar
import random
import tempfile
import shutil
import requests


class TestStructure(unittest.TestCase):

    def setUp(self):
        engine = self.engine = create_engine('sqlite:///:memory:')
        session = self.session = sessionmaker(bind=engine)()
        jaspar.Base.metadata.create_all(engine)
        # Create instance

        tf = self.tf = jaspar.TranscriptionFactor(id='MA0035', name='Gata1')
        session.add(tf)

        tf.species = jaspar.Species(id=12342345)

        tf.subunits.append(jaspar.Subunit(uniprot_id='P17679'))
        tf.subunits.append(jaspar.Subunit(uniprot_id='P17679-2'))

        tf.classes.append(jaspar.Class(name='Zinc-coordinating'))
        tf.classes.append(jaspar.Class(name='Zinc-coordinating-2'))

        tf.families.append(jaspar.Family(name='GATA'))
        tf.families.append(jaspar.Family(name='GATA-2'))

        tf.collection = jaspar.Collection(name='CORE')

        matrix = jaspar.Matrix(transcription_factor=tf, version=1)
        matrix.type = jaspar.Type(name='SELEX')

        matrix.references.append(jaspar.Resource(id=8321207))
        matrix.references.append(jaspar.Resource(id=83212072))

        matrix.positions.append(jaspar.MatrixPosition(
            position=1,
            frequency_a=0,
            frequency_c=42,
            frequency_g=12,
            frequency_t=1,
        ))

    def tearDown(self):
        jaspar.Base.metadata.drop_all(self.engine)

    def test_structure_matrix_observation(self):
        # test whatever you need to test by comparing expected to result and assert ifequal
        results = self.session.query(jaspar.TranscriptionFactor).all()

        self.assertEqual(results, [self.tf])

        tf = results[0]
        self.assertEqual(tf.id, 'MA0035')
        self.assertEqual(tf.name, 'Gata1')

        self.assertEqual(set([s.uniprot_id for s in tf.subunits]), set(['P17679', 'P17679-2']))
        self.assertEqual(tf.subunits[0].transcription_factors, [tf])

    # Run other tests for Querying such as taking a piece of real DB nand seeing if it works


class TestParseFunctions(unittest.TestCase):

    def test_make_first_col_ints(self):
        self.data1 = [['1', 'a'], ['2', 'b'], ['3', 'c']]
        result = jaspar.make_first_col_ints(self.data1)
        for i in range(0, len(result)):
            self.assertIsInstance(result[i][0], int)

    def test_make_data_int(self):
        self.data2 = [['1', 'a', '4', '7'], ['2', 'b', '5', '8'], ['3', 'c', '6', '9'], ['4', 'd', '7', '10']]
        result = jaspar.make_data_int(self.data2)
        for i in range(0, len(result)):
            self.assertIsInstance(result[i][0], int)
            self.assertIsInstance(result[i][2], int)
            self.assertIsInstance(result[i][3], float) ##change to int

    def test_sort(self):
        random.seed()
        self.data3 = [[random.random()*100, random.random()*100], [random.random()*100, random.random()*100]]
        result = jaspar.sort(self.data3)
        self.assertLess(result[0][0], result[1][0])

    def test_group(self):
        self.data4 = [[4, 'this', 2], [4, 'is', 3], [4, 'a', 1], [4, 'test',7]]
        returned = jaspar.group(self.data4,1)
        self.assertEqual(returned[4][0], ['this', 2])
        self.assertEqual(returned[4][3], ['test',7])

        self.data5 = [[4, 'A', 1, 100], [4, 'C' , 1, 123], [4, 'A', 2, 423], [4, 'C', 2, 999]]
        returned = jaspar.group(self.data5,1)
        returned[4] = jaspar.group(returned[4] , 2)
        self.assertEqual(returned[4][1][0][1] , 100)
        self.assertEqual(returned[4][1][0][0] , 'A')

        self.assertRaises(TypeError , jaspar.group(returned[4],3))


class TestQuery(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_query(self):
        src = jaspar.Jaspar(cache_dirname=self.cache_dirname, clear_content = False, load_content=False, download_backup=False, verbose = True)
        src.load_content()
        session = src.session

        z = session.query(jaspar.Matrix).get(9436)
        self.assertEqual(str(z.transcription_factor_id), 'MA0193')
        self.assertEqual(z.version, 1)
        self.assertEqual(z.type_id, 5)

        y = session.query(jaspar.TranscriptionFactor).get(3)
        self.assertEqual(str(y.id), 'MA0003')
        self.assertEqual(y.name, 'TFAP2A')
        self.assertEqual(y.species_id, 9606)
        self.assertEqual(y.collection_id, 1)

        x = session.query(jaspar.MatrixPosition).get(7)
        self.assertEqual(x.position, 7)
        self.assertEqual(x.frequency_a, 65)
        self.assertEqual(x.frequency_c, 5)
        self.assertEqual(x.frequency_g, 5)
        self.assertEqual(x.frequency_t, 22)
        self.assertEqual(x.matrix_id, 9229)

        type = session.query(jaspar.Type).filter(jaspar.Type.id == 5).first()
        self.assertEqual(type.name, 'bacterial 1-hybrid')

        family = session.query(jaspar.Family).filter(jaspar.Family.id == 44).first()
        self.assertEqual(family.name, 'Paired-related HD factors')

        class_ = session.query(jaspar.Class).filter(jaspar.Class.id == 16).first()
        self.assertEqual(class_.name, 'Homeo domain factors')

        subunit = session.query(jaspar.Subunit).filter(jaspar.Subunit.id == 5).first()
        self.assertEqual(subunit.uniprot_id, 'P17839')
