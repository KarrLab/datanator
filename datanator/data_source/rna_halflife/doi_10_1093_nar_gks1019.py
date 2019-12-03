import pandas as pd
from datanator.util import mongo_util, rna_halflife_util
import json
import datetime
import datanator.config.core
import datetime
from pymongo.collation import Collation, CollationStrength
from pymongo import errors as mongo_error


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

    def fill_uniprot(self, url, sheet_name):
        """Fill uniprot colleciton with ordered_locus_name
        from excel sheet
        
        Args:
            url (:obj:`str`): URL for Excel sheet.
            sheet_name (:obj:`str`): sheet name within Excel.

        Return:
            (:obj:`pandas.DataFrame`): Dataframe
        """
        df = self.make_df(url, sheet_name, usecols='B:D', skiprows=[0,1,2],
        names=['ordered_locus_name', 'half_life', 'r_squared'])
        row_count = len(df.index)
        for index, row in df.iterrows():
            if index == self.max_entries:
                break
            if index % 10 == 0 and self.verbose:
                print("Inserting locus {} out of {} into uniprot collection.".format(index, row_count))
            oln = row['ordered_locus_name']
            self.fill_uniprot_by_oln(oln)
        return df

    def fill_rna_halflife(self, df, species):
        """load data into rna_halflife collection
        
        Args:
            df (:obj:`pandas.DataFrame`): dataframe to be loaded into the database
            species (:obj:`list`): species name and ncbi_id
        """
        row_count = len(df.index)
        for i, row in df.iterrows():
            if i == self.max_entries:
                break
            if i % 10 == 0 and self.verbose:
                print("Processing locus {} out {}".format(i, row_count))
            halflives = {}
            oln = row['ordered_locus_name']
            halflives['halflife'] = row['half_life'] * 60
            halflives['r_sqaured'] = row['r_sqaured']
            halflives['unit'] = 's'
            halflives['reference'] = [{'doi': '10.1093/nar/gks1019', 'pubmed_id': '23125364'}]
            halflives['growth_medium'] = 'Middlebrook 7H9 with the ADC supplement (Difco) and 0.05% Tween80, at 37 degree celcius.'
            halflives['ordered_locus_name'] = oln
            halflives['species'] = species[0]
            halflives['ncbi_taxonomy_id'] = species[1]
            gene_name, protein_name = self.uniprot_query_manager.get_gene_protein_name_by_oln(oln)
            if gene_name is not None: # record exists in uniprot collection with gene_name
                self.rna_hl_collection.update_one({'gene_name': gene_name},
                                                  {'$set': {'modified': datetime.datetime.utcnow()},
                                                   '$addToSet': {'halflives': halflives,
                                                                'protein_synonyms': protein_name}}, 
                                                   collation=self.collation, upsert=True)
            elif (gene_name is None and protein_name is not None and 
                  protein_name != 'Uncharacterized protein'): # record exists in uniprot collection with non-filler protein_name
                self.rna_hl_collection.update_one({'protein_name': protein_name},
                                                  {'$set': {'modified': datetime.datetime.utcnow(),
                                                            'gene_name': gene_name},
                                                   '$addToSet': {'halflives': halflives,
                                                                'protein_synonyms': protein_name}}, 
                                                   collation=self.collation, upsert=True)
            else:
                query = {'halflives.ordered_locus_name': oln}
                doc = self.rna_hl_collection.find_one(filter=query, collation=self.collation)
                if doc is not None:
                    self.rna_hl_collection.update_one({'halflives.ordered_locus_name': oln},
                                                    {'$set': {'modified': datetime.datetime.utcnow(),
                                                              'gene_name': gene_name},
                                                    '$addToSet': {'halflives': halflives,
                                                                  'protein_synonyms': protein_name}}, 
                                                    collation=self.collation, upsert=True)
                else:
                    doc = {'halflives': [halflives], 'modified': datetime.datetime.utcnow(),
                            'gene_name': gene_name, 'protein_name': protein_name}
                    self.rna_hl_collection.insert_one(doc)