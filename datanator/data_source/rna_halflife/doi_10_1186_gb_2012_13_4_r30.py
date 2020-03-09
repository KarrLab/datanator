import pandas as pd
from datanator.util import rna_halflife_util
import datetime
import datanator.config.core
import datetime
from pymongo.collation import Collation, CollationStrength


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
        max_entries=max_entries, verbose=verbose)
        self.collation = Collation('en', strength=CollationStrength.SECONDARY)
        self.max_entries = max_entries
        self.verbose = verbose

    def load_uniprot(self):
        """Load new loci into uniprot collection
        """
        url = """https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2012-13-4-r30/MediaObjects/13059_2011_2880_MOESM3_ESM.XLSX"""
        names = ['ordered_locus_name']
        df_10987 = self.make_df(url, 'Bc10987', names=names, usecols='A', skiprows=[0,1], nrows=6014)
        df_14579 = self.make_df(url, 'Bc14579', names=names, usecols='A', skiprows=[0,1], nrows=5497)
        self.fill_uniprot_with_df(df_10987, 'ordered_locus_name')
        self.fill_uniprot_with_df(df_14579, 'ordered_locus_name')

    def fill_rna_half_life(self, df, column, species, quantification_method='Illumina GA-II'):
        """Load df into rna_halflife collection
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be loaded
            column (:obj:`str`): name of column to be used for half-life; this
            dataset has three half-life values under different methods.
            species (:obj:`list`): species name and ncbi_id.
            quantification_method (:obj:`str`): quantification method.
            'Illumina GA-II', 'RT-qPCR', or 'Roche 454'
        """
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            oln = row['ordered_locus_name']
            halflives = {}
            oln = row['ordered_locus_name']
            protein_annotation = row['protein_annotation']
            if quantification_method == 'Illumina GA-II':
                if str(row['half_life_ga_2']) == 'nan':
                    continue
                else:
                    halflives['halflife'] = row['half_life_ga_2'] * 60
                    halflives['expression_reads_per_kb_per_mb'] = row['reads_per_kb_per_mb']
            elif quantification_method == 'RT-qPCR':
                if str(row['half_life_qpcr']) == 'nan':
                    continue
                else:
                    halflives['halflife'] = row['half_life_qpcr'] * 60
            else:
                if str(row['half_life_454']) == 'nan':
                    continue
                else:
                    halflives['halflife'] = row['half_life_454'] * 60
            halflives['quantification_method'] = quantification_method
            halflives['transcriptional_start_sites'] = row['transcriptional_start_sites']
            halflives['transcriptional_end_sites'] = row['transcriptional_end_sites']
            halflives['unit'] = 's'
            if isinstance(row['operon'], float):
                halflives['operon'] = None
            else:    
                halflives['operon'] = row['operon'].split(',')
            halflives['reference'] = [{'doi': '10.1186/gb-2012-13-4-r30', 'pubmed_id': '22537947'}]
            halflives['growth_medium'] = 'Luria-Bertani (LB) broth (500 ml) at 30 degree celcius, 250 rpm.'
            halflives['ordered_locus_name'] = oln
            halflives['gene_start'] = row['gene_start']
            halflives['gene_end'] = row['gene_end']
            halflives['strand'] = row['strand']
            halflives['cog'] = row['cog']
            halflives['species'] = species[0]
            halflives['ncbi_taxonomy_id'] = species[1]

            gene_name, protein_name = self.uniprot_query_manager.get_gene_protein_name_by_oln(oln)
            
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
                query = {'halflives.ordered_locus_name': oln}
                doc = self.rna_hl_collection.find_one(filter=query, collation=self.collation)
                if doc is not None:
                    self.rna_hl_collection.update_one({'halflives.ordered_locus_name': oln},
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
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    src = Halflife(server=server, src_db=src_db,
        protein_col=protein_col, authDB='admin', readPreference='nearest',
        username=username, password=password, verbose=True, max_entries=float('inf'),
        des_db=des_db, rna_col=rna_col)
    # src.load_uniprot()

    url = 'https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2012-13-4-r30/MediaObjects/13059_2011_2880_MOESM3_ESM.XLSX'
    names = ['ordered_locus_name', 'half_life_ga_2', 'reads_per_kb_per_mb',
            'transcriptional_start_sites', 'transcriptional_end_sites', 'operon',
            'gene_start', 'gene_end', 'strand', 'gene_name', 'protein_annotation',
            'cog', 'kegg', 'half_life_qpcr', 'half_life_454']
    df = src.make_df(url, 'Bc10987', names=names, usecols='A:O', skiprows=[0,1], nrows=6014)    
    src.fill_rna_half_life(df, names, ['Bacillus cereus ATCC 10987', 222523])
    src.fill_rna_half_life(df, names, ['Bacillus cereus ATCC 10987', 222523], quantification_method='RT-qPCR')
    src.fill_rna_half_life(df, names, ['Bacillus cereus ATCC 10987', 222523], quantification_method='Roche 454')

    df = src.make_df(url, 'Bc14579', names=names, usecols='A:O', skiprows=[0,1], nrows=5497)
    src.fill_rna_half_life(df, names, ['Bacillus cereus ATCC 14579', 226900])
    src.fill_rna_half_life(df, names, ['Bacillus cereus ATCC 14579', 226900], quantification_method='RT-qPCR')
    src.fill_rna_half_life(df, names, ['Bacillus cereus ATCC 14579', 226900], quantification_method='Roche 454')  

if __name__ == '__main__':
    main()