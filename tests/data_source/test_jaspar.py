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

    def test_make_jaspar_int(self):
        self.data1 = [['1', 'a'], ['2', 'b'], ['3', 'c']]
        result = jaspar.make_jaspar_int(self.data1)
        for i in range(0, len(result)):
            self.assertIsInstance(result[i][0], int)

    def test_make_data_int(self):
        self.data2 = [['1', 'a', '4', '7'], ['2', 'b', '5', '8'], ['3', 'c', '6', '9'], ['4', 'd', '7', '10']]
        result = jaspar.make_data_int(self.data2)
        for i in range(0, len(result)):
            self.assertIsInstance(result[i][0], int)
            self.assertIsInstance(result[i][2], int)
            self.assertIsInstance(result[i][3], float)

    def test_sort(self):
        random.seed()
        self.data3 = [[random.random()*100, random.random()*100], [random.random()*100, random.random()*100]]
        result = jaspar.sort(self.data3)
        self.assertLess(result[0][0], result[1][0])

    def test_group_by_jaspar_id(self):
        self.data4 = [[4, 'this'], [4, 'is'], [4, 'a'], [4, 'test']]
        group_result, key_result = jaspar.group_by_jaspar_id(self.data4)
        self.assertEqual(group_result[0][0][0], key_result[0])


class TestCallFunctions(unittest.TestCase):

    def test_call_Attribute(self):
        self.key_attribute = [1, 2, 3]
        self.data_attribute = [[[1, 'at1']], [[2, 'at2']], [[3, 'at3']]]
        result_answer, result_key_list, result_data_list = jaspar.call_Attribute(self.key_attribute, self.data_attribute, 2)
        self.assertEqual(result_answer, ['at2'])


class TestQuery(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_query(self):
        src = jaspar.Jaspar(cache_dirname=self.cache_dirname, clear_content=True, load_content=False, download_backup=False)
        src.parse_db()
        session = src.session

        z = session.query(jaspar.MatrixObservation).get(9436)
        self.assertEqual(z.tf_name, 'schlank')
        self.assertEqual(z.collection.name, 'CORE')
        self.assertEqual(z.tf_id, 'MA0193')
        self.assertEqual(z.version, 1)
        self.assertEqual(z.matrix_id, 9436)

        type = session.query(jaspar.Type).filter(jaspar.Type.matrix_id == 9436).first()
        self.assertEqual(type.name, 'bacterial 1-hybrid')

        family = session.query(jaspar.Family).filter(jaspar.Family.matrix_id == 9433).first()
        self.assertEqual(family.name, 'Paired-related HD factors')

        class_ = session.query(jaspar.Class).filter(jaspar.Class.matrix_id == 9436).first()
        self.assertEqual(class_.name, 'Homeo domain factors')

        resource = session.query(jaspar.Resource).filter(jaspar.Resource.matrix_id == 9436).first()
        self.assertEqual(resource.medline_id, '18332042')

        species = session.query(jaspar.Species).filter(jaspar.Species.matrix_id == 9436).first()
        self.assertEqual(species.ncbi_id, '7227')

        tf = session.query(jaspar.TranscriptionFactorSubunit).filter(jaspar.TranscriptionFactorSubunit.matrix_id == 9436).first()
        self.assertEqual(tf.uniprot_id, 'Q9W423')

        matrix = session.query(jaspar.BindingMatrix).filter(jaspar.BindingMatrix.matrix_id == 9436).first()
        self.assertEqual(matrix.position, 1)
        self.assertEqual(matrix.frequency_A, 0)
        self.assertEqual(matrix.frequency_C, 19)
        self.assertEqual(matrix.frequency_G, 0)
        self.assertEqual(matrix.frequency_T, 0)
