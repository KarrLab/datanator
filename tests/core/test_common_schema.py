# -*- coding: utf-8 -*-

""" Test of common schema

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-07-31
:Copyright: 2017, Karr Lab
:License: MIT
"""
import unittest
from sqlalchemy.orm import sessionmaker
from kinetic_datanator.core import common_schema
import tempfile
import shutil


class ShortTestCommonSchema(unittest.TestCase):
    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()

    # def tearDown(self):
        # shutil.rmtree(self.cache_dirname)

    def test_working(self):
        cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname , clear_content = True,
            load_content= False, download_backup=False, max_entries = 5)
        # cs = common_schema.CommonSchema(name = 'aggregate', clear_content = True,
        #                                 load_content=False, download_backup=False,
        #                                 max_entries = 5)
        cs.load_content()
        session = cs.session

        #TODO: ADD Tests
        # deoxyuridine = session.query(cs.PhysicalEntity).filter_by(name = 'Deoxyuridine').first()
        # print deoxyuridine.observation_id
        # print deoxyuridine.name
        # print deoxyuridine.type
        # _compound = session.query(cs.Compound).get(deoxyuridine.observation_id)
        # print _compound.compound_id
        # self.assertEqual(deoxyuridine.type, 'Compound')
        # compound = session.query(common_schema.Compound).get(deoxyuridine.observation_id)
        # self.assertEqual(compound.description, "'2'-Deoxyuridine is a naturally occurring nucleoside. \
        #     It is similar in chemical structure to uridine, but without the 2'-hydroxyl group.  \
        #     It is considered to be an antimetabolite that is converted to deoxyuridine triphosphate during DNA synthesis.")
        # structure = session.query(common_schema.Structure).get(compound.structure_id)
        # self.assertEqual(structure._value_inchi, 'InChI=1S/C9H12N2O5/c12-4-6-5(13)3-8(16-6)11-2-1-7(14)10-9(11)15/h1-2,5-6,8,12-13H,3-4H2,(H,10,14,15)/t5-,6+,8+/m0/s1')
        # observe = session.query(common_schema.Observation).get(deoxyuridine.observation_id)
        # metadata = session.query(common_schema.Observation).join(observe._metadata)




# class LongTestCommonSchema(unittest.TestCase):
#     def setUp(self):
#         self.cache_dirname = tempfile.mkdtemp()
#
#     def tearDown(self):
#         shutil.rmtree(self.cache_dirname)
#
#     def test_working(self):
#         # cs = common_schema.CommonSchema(cache_dirname = self.cache_dirname , clear_content = True, load_content=False, download_backup=False)
#         cs = common_schema.CommonSchema(name = 'aggregate', clear_content = True,
#                                         load_content=False, download_backup=False
#                                         max_entries = 20)
#         cs.load_content()
#         session = cs.session
#
#
#         subunit = session.query(common_schema.Subunit).filter_by(uniprot_id = 'P41182').first()
#         self.assertEqual(subunit.entrez_id, 604)
#
#         taxon = session.query(common_schema.Taxon).get(882)
#         self.assertEqual(taxon.species_name, 'D.vulgaris')
