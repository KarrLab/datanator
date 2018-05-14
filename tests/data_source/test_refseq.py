import unittest
from kinetic_datanator.data_source import refseq
import tempfile
import shutil
from Bio import SeqIO
from os import path

class QuickTest(unittest.TestCase):

    def setUp(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.src = refseq.Refseq(cache_dirname=self.cache_dirname, download_backups=False, load_content=False)

    def tearDown(self):
        shutil.rmtree(self.cache_dirname)

    #@unittest.skip('skip until I figure out file path')
    def test_load_data(self):





        src = self.src
        file = "{}/test_mpn_sequence.gb".format(path.dirname(__file__))
        bio_seqio_object = SeqIO.parse(file, "genbank")
        list_of_bio_seqio_objects = [bio_seqio_object]
        src.load_content(list_of_bio_seqio_objects)
        session = src.session

        ref_genome = session.query(refseq.ReferenceGenome).filter_by(version="NC_000912.1").first()
        self.assertEqual(ref_genome.accessions[0].id, "NC_000912")
        self.assertEqual(ref_genome.organism, "Mycoplasma pneumoniae M129")
        self.assertEqual(len(ref_genome.genes), 690)

        self.assertEqual(len(session.query(refseq.Gene).all()), 690)
        gene = session.query(refseq.Gene).filter_by(name="lip2").first()
        self.assertEqual(gene.gene_synonyms[0].name, "P01_orf268")
        gene = session.query(refseq.Gene).filter_by(name="valS").first()
        self.assertEqual(gene.ec_numbers[0].ec_number, "6.1.1.9")
        the_identifier = None
        for identifier in gene.identifiers:
            if identifier.namespace == 'GeneID':
                the_identifier = identifier
        self.assertEqual(the_identifier.name, '877380')
        self.assertEqual(gene.locus_tag, "MPN480")
        #self.assertEqual(gene.essentiality, "E")

    def test_upload_data_from_kegg_org_symbol(self):
        src = self.src
        src.upload_data_from_kegg_org_symbol("ell")
        session = src.session
        ref_genome = session.query(refseq.ReferenceGenome).filter_by(version="CP002967.1").first()
        self.assertEqual(ref_genome.accessions[0].id, "CP002967")
        self.assertEqual(ref_genome.organism, "Escherichia coli W")

    @unittest.skip('too long')
    def test_new_method(self):
        src = self.src
        src.upload_ref_seq_for_all_prokaryotic_kegg_org()

