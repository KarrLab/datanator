# -*- coding: utf-8 -*-

""" Test of CORUM database

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Saahith Pochiraju <saahith116@gmail.com>
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


class TestCorumNoSQL(unittest.TestCase):
 
    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.db = 'test'
        self.MongoDB = 'mongodb://mongo:27017/'

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_con_db(self):
        src = corum_nosql.CorumNoSQL(
            self.cache_dirname, self.MongoDB, self.db, replicaSet=None, verbose = True, max_entries = 20)
        client, db, collection = src.con_db('corum')
        self.assertNotEqual(collection, 'Server not available')
        client.close()

    #@unittest.skip("loading everything")
    def test_load_some_content(self):
        src = corum_nosql.CorumNoSQL(
            self.cache_dirname, self.MongoDB, self.db, replicaSet=None, verbose = True, max_entries = 20)
        client, _, collection = src.load_content()
        self.assertEqual(collection.find().count(), 20)
        cursor = collection.find({'subunits(UniProt IDs)': 'P41182'}).limit(3)
        self.assertEqual(cursor.count(), 3)
        self.assertEqual(cursor[1]['ComplexName'], 'BCL6-HDAC5 complex')
        self.assertEqual(cursor[2]['subunits(Protein name)'], ['B-cell lymphoma 6 protein', 'Histone deacetylase 7'])
        collection.drop()
        
    @unittest.skip("will not work on circle ci due to file directory setting")
    def test_load_all_content(self):
        db = 'datanator'
        cache_dirname = '../../datanator/data_source/cache'
        src = corum_nosql.CorumNoSQL(
            cache_dirname, self.MongoDB, db, verbose = True)
        client, _, collection = src.load_content()

        c = collection.find_one({'ComplexID':80})
        self.assertEqual(c['ComplexName'], 'Ubiquitin E3 ligase (SKP1A, SKP2, CUL1, RBX1)')

        s = collection.find_one({'subunits(UniProt IDs)':'Q9UQL6'})
        self.assertEqual(s['subunits(Protein name)'][1], 'Histone deacetylase 5')
        collection.drop()
        client.close()

        # t = session.query(corum.Taxon).filter(corum.Taxon.ncbi_id == 9606).first()
        # self.assertEqual(t.swissprot_id, 'Homo sapiens (Human)')

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
