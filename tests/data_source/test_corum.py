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
from datanator.data_source import corum
from sqlalchemy.orm import sessionmaker
import tempfile
import shutil


class TestCorumDBCreation(unittest.TestCase):
    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    def test_load_some_content(self):
        src = corum.Corum(cache_dirname=self.cache_dirname, load_content=False, download_backups=False, max_entries=10)
        src.load_content()
        session = src.session

        subunit = session.query(corum.Subunit).get(3)
        self.assertEqual(subunit.su_uniprot, 'P41182')
        self.assertEqual(str(subunit.protein_name), 'B-cell lymphoma 6 protein')

        complx = subunit.complex

        self.assertEqual(complx.complex_id, 2)
        self.assertEqual(complx.complex_name, 'BCL6-HDAC5 complex')

        obs = complx.observation
        self.assertEqual(obs.id, 2)
        self.assertEqual(obs.pubmed_id, 11929873)

        tax = obs.taxon

        self.assertEqual(tax.swissprot_id, 'Homo sapiens (Human)')

    def test_load_all_content(self):
        src = corum.Corum(cache_dirname=self.cache_dirname, load_content=False, download_backups=False, verbose=False)
        src.load_content()
        session = src.session

        c = session.query(corum.Complex).get(80)
        self.assertEqual(c.complex_name, 'Ubiquitin E3 ligase (SKP1A, SKP2, CUL1, RBX1)')

        s = session.query(corum.Subunit).filter(corum.Subunit.su_uniprot == 'Q9UQL6').first()
        self.assertEqual(s.protein_name, 'Histone deacetylase 5')

        t = session.query(corum.Taxon).filter(corum.Taxon.ncbi_id == 9606).first()
        self.assertEqual(t.swissprot_id, 'Homo sapiens (Human)')

    def test_parse_list(self):
        self.assertEqual(corum.parse_list(None), [None])
        self.assertEqual(corum.parse_list('protein-A'), ['protein-A'])
        self.assertEqual(corum.parse_list('protein-A;protein-B'), ['protein-A', 'protein-B'])
        self.assertEqual(corum.parse_list('protein-A [A;B;C];protein-B'), ['protein-A [A;B;C]', 'protein-B'])

    def test_correct_protein_name_list(self):
        self.assertEqual(len(corum.parse_list(corum.correct_protein_name_list(
            'Nuclear pore complex protein Nup98-Nup96 [Cleaved into: Nuclear pore complex protein Nup98;'
            'Protein SEC13 homolog;'
            'Nuclear pore complex protein Nup107;'
            'Nuclear pore complex protein Nup160;'
            'Nucleoporin Nup43 ;'
            'Nucleoporin Nup37 ;'
            'Nuclear pore complex protein Nup133;'
            'Nucleoporin SEH1;'
            'Nuclear pore complex protein Nup85'))), 9)
        self.assertEqual(len(corum.parse_list(corum.correct_protein_name_list(
            'Importin subunit alpha-4;'
            ' Prelamin-A/C[Cleaved into: Lamin-A/C;'
            ' Nucleophosmin;'
            'Poly [ADP-ribose] polymerase 1;'
            'Histone H2A.Z ;'
            'DNA topoisomerase 2-alpha;'
            'Transcriptional repressor CTCF ;'
            'Importin subunit alpha-5;'
            'Histone H2A type 2-A'))), 9)
        self.assertEqual(len(corum.parse_list(corum.correct_protein_name_list(
            'Prelamin-A/C [Cleaved into: Lamin-A/C ;'
            'Nuclear transcription factor Y subunit alpha'))), 2)

        self.assertEqual(len(corum.parse_list(corum.correct_protein_name_list(
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
