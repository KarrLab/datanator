"""
This code tests all aspects of jaspar.py

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017 July 24th
:Copyright: 2017, Karr Lab
:License: MIT

"""

import unittest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import jaspar
from jaspar import Base
import random


class TestStructure(unittest.TestCase):

    def setUp(self):
        self.engine = create_engine('sqlite:///:memory:')
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        Base.metadata.create_all(self.engine)
        #Create instance
        self.type = jaspar.Type(type_name = 'SELEX', jaspar_matrix_ID = 1)
        self.session.add(self.type)
        self.family = jaspar.Family(family_name = u'pochiraju', jaspar_matrix_ID = 1)
        self.session.add(self.family)
        self.classes = jaspar.Class(class_name = 'Zinc-Fingers', jaspar_matrix_ID = 1)
        self.session.add(self.classes)
        self.resources = jaspar.Resources(medline_id = '0123456789', jaspar_matrix_ID = 1)
        self.session.add(self.resources)
        self.species = jaspar.Species(NCBI_id = '12342345', jaspar_matrix_ID = 1)
        self.session.add(self.species)
        self.tf = jaspar.TranscriptionFactor(uniprot_id = 'P612034', jaspar_matrix_ID = 1)
        self.session.add(self.tf)
        self.data = jaspar.BindingMatrix(
                    position = 1,
                    frequency_A = 0.0,
                    frequency_C = 42.0,
                    frequency_G = 12.0,
                    frequency_T = 1.0,
                    jaspar_matrix_ID = 1
                )
        self.session.add(self.data)
        self.matrixobservation = jaspar.MatrixObservation(
                    jaspar_matrix_ID = 1,
                    jaspar_tf_ID = 'MA0437',
                    version = 1,
                    tf_name = 'YPR196W',
                    jaspar_collection = 'CORE',
                    data = [self.data],
                    references = [self.resources],
                    type_info = [self.type],
                    class_info = [self.classes],
                    family_info = [self.family],
                    tf = [self.tf],
                    species = [self.species]
                )
        self.session.add(self.matrixobservation)
        self.session.commit()

    def tearDown(self):
        Base.metadata.drop_all(self.engine)

    def test_structure_matrixobservation(self):
        #test whatever you need to test by comparing expected to result and assert ifequal
        expected = [self.matrixobservation.jaspar_matrix_ID]
        result = self.session.query(jaspar.MatrixObservation.jaspar_matrix_ID).all()
        ## Outputs list of tuple. Convert to int ##
        result = [i[0] for i in result]
        self.assertEqual(result, expected)

    ## Run other tests for Querying such as taking a piece of real DB nand seeing if it works

class TestParseFunctions(unittest.TestCase):

    def test_make_jaspar_int(self):
        self.data1 = [['1','a'],['2','b'],['3','c']]
        result = jaspar.make_jaspar_int(self.data1)
        for i in range(0,len(result)):
            self.assertIsInstance(result[i][0], int)

    def test_make_data_int(self):
        self.data2 = [['1','a','4','7'],['2','b','5','8'],['3','c','6','9'],['4','d','7','10']]
        result = jaspar.make_data_int(self.data2)
        for i in range(0,len(result)):
            self.assertIsInstance(result[i][0],int)
            self.assertIsInstance(result[i][2],int)
            self.assertIsInstance(result[i][3],float)

    def test_sort(self):
        random.seed()
        self.data3 = [[random.random()*100,random.random()*100], [random.random()*100,random.random()*100]]
        result = jaspar.sort(self.data3)
        self.assertLess(result[0][0],result[1][0])

    def test_group_by_jaspar_id(self):
        self.data4 = [[4,'this'],[4,'is'],[4,'a'],[4,'test']]
        group_result, key_result = jaspar.group_by_jaspar_id(self.data4)
        self.assertEqual(group_result[0][0][0],key_result[0])

class TestCallFunctions(unittest.TestCase):
    def test_call_Attribute(self):
        self.key_attribute = [1,2,3]
        self.data_attribute = [[[1,'at1']],[[2,'at2']],[[3,'at3']]]
        result_answer, result_key_list, result_data_list = jaspar.call_Attribute(self.key_attribute, self.data_attribute, 2)
        self.assertEqual(result_answer, ['at2'])


class TestQuery(unittest.TestCase):
    def setUp(self):
        self.engine = create_engine('sqlite:///:memory:')
        ## Change from uniprot to strings
        self.engine.raw_connection().connection.text_factory = str
        self.Session = sessionmaker(bind=self.engine)
        self.session = self.Session()
        Base.metadata.drop_all(self.engine)
        Base.metadata.create_all(self.engine)
        database_url = 'http://jaspar.genereg.net/html/DOWNLOAD/database/'
        jaspar.parse_Jaspar_db(self.session, database_url)
        self.session.commit()
        print 'Succesfully Committed'



    def test_query(self):
        self.z = self.session.query(jaspar.MatrixObservation).get(9436)
        self.assertEqual(self.z.tf_name,'schlank\n')
        self.assertEqual(self.z.jaspar_collection,'CORE')
        self.assertEqual(self.z.jaspar_tf_ID, 'MA0193')
        self.assertEqual(self.z.version, 1)
        self.assertEqual(self.z.jaspar_matrix_ID, 9436)

        self.type = self.session.query(jaspar.Type).filter(jaspar.Type.jaspar_matrix_ID == 9436).first()
        self.assertEqual(self.type.type_name, 'bacterial 1-hybrid')

        self.family = self.session.query(jaspar.Family).filter(jaspar.Family.jaspar_matrix_ID == 9433).first()
        self.assertEqual(self.family.family_name, 'Paired-related HD factors')

        self.class_ = self.session.query(jaspar.Class).filter(jaspar.Class.jaspar_matrix_ID == 9436).first()
        self.assertEqual(self.class_.class_name, 'Homeo domain factors')

        self.resource = self.session.query(jaspar.Resources).filter(jaspar.Resources.jaspar_matrix_ID == 9436).first()
        self.assertEqual(self.resource.medline_id, '18332042')

        self.species = self.session.query(jaspar.Species).filter(jaspar.Species.jaspar_matrix_ID == 9436).first()
        self.assertEqual(self.species.NCBI_id,'7227')

        self.tf = self.session.query(jaspar.TranscriptionFactor).filter(jaspar.TranscriptionFactor.jaspar_matrix_ID == 9436).first()
        self.assertEqual(self.tf.uniprot_id, 'Q9W423')

        self.data = self.session.query(jaspar.BindingMatrix).filter(jaspar.BindingMatrix.jaspar_matrix_ID == 9436).first()
        self.assertEqual(self.data.position, 1)
        self.assertEqual(self.data.frequency_A, 0)
        self.assertEqual(self.data.frequency_C, 19)
        self.assertEqual(self.data.frequency_G, 0)
        self.assertEqual(self.data.frequency_T, 0)

    def tearDown(self):
        Base.metadata.drop_all(self.engine)


if __name__ == '__main__':
    unittest.main()
