# -*- coding: utf-8 -*-

""" Test of CORUM database

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-08-13
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

import unittest
from datanator.data_source import corum_nosql
import tempfile
import shutil
import pymongo
import datanator.config.core

class TestCorumNoSQL(unittest.TestCase):
 
    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.db = 'test'
        self.username = datanator.config.core.get_config()['datanator']['mongodb']['user']
        self.password = datanator.config.core.get_config()['datanator']['mongodb']['password']
        self.MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
        self.port = datanator.config.core.get_config()['datanator']['mongodb']['port']
        self.replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    # @unittest.skip("loading everything")
    def test_load_some_content(self):
        src = corum_nosql.CorumNoSQL(self.MongoDB, self.db, replicaSet=self.replSet, cache_dirname = self.cache_dirname,
            verbose = True, max_entries=20, username = self.username, password = self.password)
        collection = src.load_content()
        self.assertTrue(collection.count_documents({}) in [18, 19, 20])
        cursor = collection.find({'subunits_uniprot_id': 'P41182'}).limit(3)
        self.assertEqual(cursor.count(), 3)
        self.assertEqual(cursor[1]['complex_id'], 2)
        self.assertEqual(cursor[2]['subunits_protein_name'], ['B-cell lymphoma 6 protein', 'Histone deacetylase 7'])
        src.client.close()

    @unittest.skip(" loading all contents")
    def test_load_all_content(self):
        db = 'datanator'
        src = corum_nosql.CorumNoSQL(self.MongoDB, db, verbose = True, username = self.username, password = self.password,
                                    cache_dirname = self.cache_dirname)
        collection = src.load_content()
        collection = src.load_content(endpoint='splice')

        c = collection.find_one({'ComplexID':'6620'})
        self.assertEqual(c['ComplexName'], 'MGAT1-SLC35A2(splice variant 2) complex')


    def test_parse_list(self):
        self.assertEqual(corum_nosql.parse_list(None), [None])
        self.assertEqual(corum_nosql.parse_list('protein-A'), ['protein-A'])
        self.assertEqual(corum_nosql.parse_list('protein-A;protein-B'), ['protein-A', 'protein-B'])
        self.assertEqual(corum_nosql.parse_list('protein-A [A;B;C];protein-B'), ['protein-A [A;B;C]', 'protein-B'])

    def test_correct_protein_name_list(self):
        self.assertEqual(len(corum_nosql.parse_list(corum_nosql.correct_protein_name_list(
            'Nuclear pore complex protein Nup98-Nup96 [Cleaved into: Nuclear pore complex protein Nup98;'
            'Protein SEC13 homolog;'
            'Nuclear pore complex protein Nup107;'
            'Nuclear pore complex protein Nup160;'
            'Nucleoporin Nup43 ;'
            'Nucleoporin Nup37 ;'
            'Nuclear pore complex protein Nup133;'
            'Nucleoporin SEH1;'
            'Nuclear pore complex protein Nup85'))), 9)
        self.assertEqual(len(corum_nosql.parse_list(corum_nosql.correct_protein_name_list(
            'Importin subunit alpha-4;'
            ' Prelamin-A/C[Cleaved into: Lamin-A/C;'
            ' Nucleophosmin;'
            'Poly [ADP-ribose] polymerase 1;'
            'Histone H2A.Z ;'
            'DNA topoisomerase 2-alpha;'
            'Transcriptional repressor CTCF ;'
            'Importin subunit alpha-5;'
            'Histone H2A type 2-A'))), 9)
        self.assertEqual(len(corum_nosql.parse_list(corum_nosql.correct_protein_name_list(
            'Prelamin-A/C [Cleaved into: Lamin-A/C ;'
            'Nuclear transcription factor Y subunit alpha'))), 2)

        self.assertEqual(len(corum_nosql.parse_list(corum_nosql.correct_protein_name_list(
            'Transcription initiation factor TFIID subunit 4;'
            ' Maltase-glucoamylase, intestinal [Includes: Maltase ;'
            'Transcription factor E2F6;'
            'Transcription initiation factor TFIID subunit 1;'
            'Heat shock 70 kDa protein 4;'
            'Transcription initiation factor TFIID subunit 6;'
            'Host cell factor 1;'
            'WD repeat-containing protein 5;'
            'Histone-lysine N-methyltransferase 2A;'
            'Retinoblastoma-binding protein 5;'
            'Transcription initiation factor TFIID subunit 7;'
            'Transcription initiation factor TFIID subunit 9;'
            'INO80 complex subunit C ;'
            'KAT8 regulatory NSL complex subunit 1;'
            'Chromatin complexes subunit BAP18 ;'
            'Proline-, glutamic acid- and leucine-rich protein 1;'
            'U4/U6 small nuclear ribonucleoprotein Prp31 ;'
            'Microspherule protein 1;'
            'E3 ubiquitin-protein ligase RING2;'
            'PHD finger protein 20;'
            'Sentrin-specific protease 3;'
            'Histone acetyltransferase KAT8;'
            'Chromodomain-helicase-DNA-binding protein 8 ;'
            'Testis-expressed sequence 10 protein;'
            'Set1/Ash2 histone methyltransferase complex subunit ASH2;'
            'RuvB-like 1;'
            'Ribosomal biogenesis protein LAS1L'))), 27)

    def test_parse_subunits(self):
        subunits = ['P78381-2', "P01234"]
        result = corum_nosql.parse_subunits(subunits)
        exp = ['P78381', 'P01234']
        self.assertEqual(result, exp)
        subunits = [None]
        result = corum_nosql.parse_subunits(subunits)
        self.assertEqual(result, [None])