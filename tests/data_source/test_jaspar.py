"""
This module tests all aspects of jaspar.py

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Jonathan Karr <jonrkar@gmail.com>
:Date: 2017-07-24
:Copyright: 2017, Karr Lab
:License: MIT

"""

from kinetic_datanator.data_source import jaspar
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import random
import requests
import shutil
import tempfile
import unittest


class TestStructure(unittest.TestCase):

    def setUp(self):
        engine = create_engine('sqlite:///:memory:')
        self.session = sessionmaker(bind=engine)()
        jaspar.Base.metadata.create_all(engine)

    def test_structure_matrix_observation(self):
        session = self.session

        # create object instances
        tf = jaspar.TranscriptionFactor(id='MA0035', name='Gata1')
        session.add(tf)

        tf.species.append(jaspar.Species(ncbi_id=12342345))

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

        # test that the relationships among the classes are implemented correctly
        results = session.query(jaspar.TranscriptionFactor).all()

        self.assertEqual(results, [tf])

        tf = results[0]
        self.assertEqual(tf.id, 'MA0035')
        self.assertEqual(tf.name, 'Gata1')

        self.assertEqual(set([s.uniprot_id for s in tf.subunits]), set(['P17679', 'P17679-2']))
        self.assertEqual(tf.subunits[0].transcription_factors, [tf])


class TestParseFunctions(unittest.TestCase):

    def test_type_cast_matrix_ids_to_ints(self):
        data = [
            ['1', 'a'],
            ['2', 'b'],
            ['3', 'c'],
        ]
        table = jaspar.type_cast_matrix_ids_to_ints(data)
        for row in table:
            self.assertIsInstance(row[0], int)

    def test_type_cast_matrix_positions_and_frequencies(self):
        data = [
            ['1', 'a', '4', '7'],
            ['2', 'b', '5', '8'],
            ['3', 'c', '6', '9'],
            ['4', 'd', '7', '10']
        ]
        table = jaspar.type_cast_matrix_positions_and_frequencies(data)
        for row in table:
            self.assertIsInstance(row[2], int)
            self.assertIsInstance(row[3], float)

    def test_group_by_matrix_ids(self):
        table = [
            [4, 'this', 2],
            [4, 'is', 3],
            [4, 'a', 1],
            [4, 'test', 7],
            [3, 'this', 2],
            [3, 'is', 3],
            [5, 'a', 1],
            [5, 'test', 7],
        ]
        result = jaspar.group_by_matrix_ids(table)
        self.assertEqual(len(result), 3)
        self.assertEqual(set([tuple(row) for row in result[4]]), set([('this', 2), ('is', 3), ('a', 1), ('test', 7)]))
        self.assertEqual(set([tuple(row) for row in result[3]]), set([('this', 2), ('is', 3)]))
        self.assertEqual(set([tuple(row) for row in result[5]]), set([('a', 1), ('test', 7)]))

    def test_group_by_position(self):
        table = [
            ['A', 1, 100],
            ['C', 1, 123],
            ['A', 2, 423],
            ['C', 2, 999],
        ]
        result = jaspar.group_by_position(table)
        self.assertEqual(len(result), 2)
        self.assertEqual(set([tuple(base) for base in result[1]]), set([('A', 100), ('C', 123)]))
        self.assertEqual(set([tuple(base) for base in result[2]]), set([('A', 423), ('C', 999)]))

    def test_group_by_matrix_ids_and_position(self):
        table = [
            [4, 'A', 1, 410],
            [4, 'C', 1, 411],
            [4, 'A', 2, 420],
            [4, 'C', 2, 421],
            [3, 'A', 1, 310],
            [3, 'C', 1, 311],
            [3, 'A', 2, 320],
            [3, 'C', 2, 321],
            [5, 'A', 1, 510],
            [5, 'C', 1, 511],
            [5, 'A', 2, 520],
            [5, 'C', 2, 521],
        ]
        result = jaspar.group_by_matrix_ids(table)
        self.assertEqual(len(result), 3)

        result2 = jaspar.group_by_position(result[4])
        self.assertEqual(len(result2), 2)
        self.assertEqual(set([tuple(base) for base in result2[1]]), set([('A', 410), ('C', 411)]))
        self.assertEqual(set([tuple(base) for base in result2[2]]), set([('A', 420), ('C', 421)]))


class TestLoadContent(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test(self):
        src = jaspar.Jaspar(cache_dirname=self.cache_dirname, clear_content=False, load_content=False, download_backup=False, verbose=True)
        src.load_content()
        session = src.session

        q = session.query(jaspar.Matrix).get(9436)
        self.assertEqual(str(q.transcription_factor_id), 'MA0193')
        self.assertEqual(q.version, 1)
        self.assertEqual(q.type_id, 5)

        q = session.query(jaspar.TranscriptionFactor).get(3)
        self.assertEqual(str(q.id), 'MA0003')
        self.assertEqual(q.name, 'TFAP2A')
        self.assertEqual(q.collection_id, 1)

        q = session.query(jaspar.MatrixPosition).get(7)
        self.assertEqual(q.position, 7)
        self.assertEqual(q.frequency_a, 65)
        self.assertEqual(q.frequency_c, 5)
        self.assertEqual(q.frequency_g, 5)
        self.assertEqual(q.frequency_t, 22)
        self.assertEqual(q.matrix_id, 9229)

        type = session.query(jaspar.Type).filter(jaspar.Type.id == 5).first()
        self.assertEqual(type.name, 'bacterial 1-hybrid')

        family = session.query(jaspar.Family).filter(jaspar.Family.id == 44).first()
        self.assertEqual(family.name, 'Paired-related HD factors')

        class_ = session.query(jaspar.Class).filter(jaspar.Class.id == 16).first()
        self.assertEqual(class_.name, 'Homeo domain factors')

        subunit = session.query(jaspar.Subunit).filter(jaspar.Subunit.id == 5).first()
        self.assertEqual(subunit.uniprot_id, 'P17839')
