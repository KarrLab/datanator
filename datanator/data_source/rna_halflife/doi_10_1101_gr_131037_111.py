import pandas as pd
from datanator.util import rna_halflife_util, file_util
import datetime
import datanator.config.core
import datetime
from pymongo.collation import Collation, CollationStrength
import tempfile
import shutil


class Halflife(rna_halflife_util.RnaHLUtil):

    def __init__(self, cache_dir=None, server=None, src_db=None, protein_col=None,
                authDB=None, readPreference=None, username=None, password=None,
                verbose=None, max_entries=None, des_db=None, rna_col=None):
        """Init
        
        Args:
            cache_dir (:obj:`str`, optional): Cache directory for logs. Defaults to None.
            server (:obj:`str`, optional): MongoDB server address. Defaults to None.
            db (:obj:`str`, optional): Database where initial uniprot collection resides. Defaults to None.
            collection_str (:obj:`str`, optional): name of collection. Defaults to None.
            authDB (:obj:`str`, optional): MongoDB authentication database. Defaults to None.
            readPreference (:obj:`str`, optional): MongoDB read preference. Defaults to None.
            username (:obj:`str`, optional): MongoDB username. Defaults to None.
            password (:obj:`str`, optional): MongoDB password. Defaults to None.
            verbose (:obj:`bool`, optional): Wheter to display verbose messages. Defaults to None.
            max_entries (:obj:`int`, optional): Number of records to be processed. Defaults to None.
            uniprot_col_db (:obj:`int`, optional): Database to which new uniprot records will be inserted. Defaults to None.
        """
        super().__init__(server=server, username=username, password=password, src_db=src_db,
        des_db=des_db, protein_col=protein_col, rna_col=rna_col, authDB=authDB, readPreference=readPreference,
        max_entries=max_entries, verbose=verbose, cache_dir=cache_dir)
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)
        self.max_entries = max_entries
        self.verbose = verbose

    def load_uniprot(self):
        """Load new loci into uniprot collection
        """
        url = """https://genome.cshlp.org/content/suppl/2012/02/06/gr.131037.111.DC1/Supp_Table_2.xlsx"""
        names = ['Accession ID']
        df_10987 = self.make_df(url, 'V1ncodemouse_probe_annotations_', names=names, usecols='B', skiprows=list(range(0,27240)), nrows=34509)
        self.fill_uniprot_with_df(df_10987, 'Accession ID', identifier_type='sequence_embl', species=[10090])

    def fill_rna_half_life(self, df, species):
        """Load df into rna_halflife collection
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be loaded
            species (:obj:`list`): species name and ncbi_id.
        """
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            halflives = {}
            oln = row['Accession ID']
            protein_annotation = row['UniGene Name']

            halflives['transcript_size'] = row['Transcript Size']
            halflives['cds_size'] = row['CDS Size']
            halflives['intron_size'] = row['Intron Size']
            halflives['genomic_size'] = row['Genomic Size']
            halflives['intron_count'] = row['Intron Count']
            halflives['halflife'] = row['HalfLife (min)'] * 60
            halflives['r_sqaured'] = row['RSquared']
            halflives['standard_error'] = row['StandardErrorK']
            halflives['unit'] = 's'

            halflives['reference'] = [{'doi': '10.1101/gr.131037.111', 'pubmed_id': '22406755'}]
            halflives['accession_id'] = oln.split(',')
            halflives['ncbi_taxonomy_id'] = species[1]

            gene_name, protein_name = self.uniprot_query_manager.get_gene_protein_name_by_embl(oln.split(','), species=[10090])
            
            if gene_name is not None: # record exists in uniprot collection with gene_name
                self.rna_hl_collection.update_one({'gene_name': gene_name},
                                                  {'$set': {'modified': datetime.datetime.utcnow()},
                                                   '$addToSet': {'halflives': halflives,
                                                                'protein_synonyms': {'$each': [protein_annotation, protein_name]}}}, 
                                                   collation=self.collation, upsert=True)
            elif (gene_name is None and protein_name is not None and 
                  protein_name != 'Uncharacterized protein'): # record exists in uniprot collection with non-filler protein_name
                self.rna_hl_collection.update_one({'protein_name': protein_name},
                                                  {'$set': {'modified': datetime.datetime.utcnow(),
                                                            'gene_name': gene_name},
                                                   '$addToSet': {'halflives': halflives,
                                                                'protein_synonyms': {'$each': [protein_annotation, protein_name]}}}, 
                                                   collation=self.collation, upsert=True)
            else:
                query = {'halflives.accession_id': {'$in': oln}}
                doc = self.rna_hl_collection.find_one(filter=query, collation=self.collation)
                if doc is not None:
                    self.rna_hl_collection.update_one({'halflives.accession_id': halflives['accession_id']},
                                                    {'$set': {'modified': datetime.datetime.utcnow(),
                                                              'gene_name': gene_name},
                                                    '$addToSet': {'halflives': halflives,
                                                                  'protein_synonyms': protein_annotation}}, 
                                                    collation=self.collation, upsert=True)
                else:
                    doc = {'halflives': [halflives], 'modified': datetime.datetime.utcnow(),
                            'gene_name': gene_name, 'protein_name': protein_name}
                    self.rna_hl_collection.insert_one(doc)

def main():
    src_db = 'datanator'
    des_db = 'datanator'
    rna_col = 'rna_halflife'
    protein_col = 'uniprot'
    cache_dir = tempfile.mkdtemp()
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    src = Halflife(server=server, src_db=src_db,
        protein_col=protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=float('inf'),
        des_db=des_db, rna_col=rna_col, cache_dir=cache_dir)
    url = 'https://genome.cshlp.org/content/suppl/2012/02/06/gr.131037.111.DC1/Supp_Table_2.xlsx'
    # src.load_uniprot()
    # names = ['ordered_locus_name', 'half_life_ga_2', 'reads_per_kb_per_mb',
    #         'transcriptional_start_sites', 'transcriptional_end_sites', 'operon',
    #         'gene_start', 'gene_end', 'strand', 'gene_name', 'protein_annotation',
    #         'cog', 'kegg', 'half_life_qpcr', 'half_life_454']
    # df = src.make_df(url, 'V1ncodemouse_probe_annotations_', names=names, usecols='A:O', skiprows=[0,1], nrows=34509)
    # src.fill_rna_half_life(df, names, ['Mus musculus', 10090])
    # shutil.rmtree(cache_dir)

if __name__ == '__main__':
    main()