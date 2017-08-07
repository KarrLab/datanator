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

        matrix.references.append(jaspar.Resource(pubmed_id=8321207))
        matrix.references.append(jaspar.Resource(pubmed_id=83212072))

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

    def test_import_a_few_tfs(self):
        src = jaspar.Jaspar(cache_dirname=self.cache_dirname, load_content=False, download_backup=False, max_entries=10)
        src.load_content()
        session = src.session

        matrix = session.query(jaspar.Matrix).get(9234)
        self.assertEqual(matrix.version, 1)
        self.assertEqual(matrix.type.name, 'SELEX')
        self.assertEqual(len(matrix.references), 1)
        self.assertEqual(matrix.references[0].pubmed_id, 7592839)

        tf = matrix.transcription_factor
        self.assertEqual(tf.id, 'MA0006')
        self.assertEqual(tf.name, 'Ahr::Arnt')
        self.assertEqual(tf.collection.name, 'CORE')

        self.assertEqual(set([subunit.uniprot_id for subunit in tf.subunits]), set(['P30561', 'P53762']))

        self.assertEqual(len(tf.classes), 1)
        self.assertEqual(tf.classes[0].name, 'Basic helix-loop-helix factors (bHLH)')

        self.assertEqual(len(tf.families), 1)
        self.assertEqual(tf.families[0].name, 'PAS domain factors')

        self.assertEqual(len(tf.species), 1)
        self.assertEqual(tf.species[0].ncbi_id, 10090)

        self.assertEqual(len(matrix.positions), 6)
        pos = next(pos for pos in matrix.positions if pos.position == 3)
        self.assertEqual(pos.frequency_a, 0)
        self.assertEqual(pos.frequency_c, 23)
        self.assertEqual(pos.frequency_g, 0)
        self.assertEqual(pos.frequency_t, 1)

    def test_import_all(self):
        src = jaspar.Jaspar(cache_dirname=self.cache_dirname, load_content=False, download_backup=False)
        src.load_content()
        session = src.session

        matrix = session.query(jaspar.Matrix).get(9436)
        self.assertEqual(matrix.version, 1)
        self.assertEqual(matrix.type.name, 'bacterial 1-hybrid')
        self.assertEqual(len(matrix.references), 1)
        self.assertEqual(matrix.references[0].pubmed_id, 18332042)

        tf = matrix.transcription_factor
        self.assertEqual(tf.id, 'MA0193')
        self.assertEqual(tf.name, 'schlank')
        self.assertEqual(tf.collection.name, 'CORE')

        self.assertEqual(len(tf.subunits), 1)
        self.assertEqual(tf.subunits[0].uniprot_id, 'Q9W423')

        self.assertEqual(len(tf.classes), 1)
        self.assertEqual(tf.classes[0].name, 'Homeo domain factors')

        self.assertEqual(len(tf.families), 0)

        self.assertEqual(len(tf.species), 1)
        self.assertEqual(tf.species[0].ncbi_id, 7227)

        self.assertEqual(len(matrix.positions), 7)
        pos = next(pos for pos in matrix.positions if pos.position == 7)
        self.assertEqual(pos.frequency_a, 12)
        self.assertEqual(pos.frequency_c, 0)
        self.assertEqual(pos.frequency_g, 4)
        self.assertEqual(pos.frequency_t, 3)
