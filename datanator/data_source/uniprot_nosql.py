'''
    Generates uniprot_swiss (reviewed) NoSQL documents
    load documents into MongoDB collections

:Author: Zhouyang Lian <zhouyang.lian@familian.life>
:Author: Jonathan <jonrkarr@gmail.com>
:Date: 2019-04-02
:Copyright: 2019, Karr Lab
:License: MIT
'''

import io
import math
import json
import pandas
import requests
import pymongo.errors
from datanator.util import mongo_util
from pymongo.collation import Collation, CollationStrength
from datanator_query_python.query import query_taxon_tree, query_kegg_orthology, query_sabiork
import datanator.config.core


class UniprotNoSQL(mongo_util.MongoUtil):
    def __init__(self, MongoDB=None, db=None, max_entries=float('inf'), verbose=False,
         username=None, password=None, authSource='admin', replicaSet=None, collection_str='uniprot'):
        self.url = 'http://www.uniprot.org/uniprot/?fil=reviewed:yes'
        self.query_url = 'https://www.uniprot.org/uniprot/?query='
        self.MongoDB = MongoDB
        self.db = db
        self.max_entries = max_entries
        self.collection_str = collection_str
        super(UniprotNoSQL, self).__init__(MongoDB=MongoDB, db=db, username=username,
                                 password=password, authSource=authSource, replicaSet=replicaSet,
                                 verbose=verbose, max_entries=max_entries)
        self.taxon_manager = query_taxon_tree.QueryTaxonTree(username=username, MongoDB=MongoDB, password=password,
                                                        authSource=authSource)
        self.ko_manager = query_kegg_orthology.QueryKO(username=username, password=password, server=MongoDB,
        authSource=authSource, verbose=verbose, max_entries=max_entries)
        self.sabio_manager = query_sabiork.QuerySabio(MongoDB=MongoDB, username=username, password=password,
        authSource=authSource)
        self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)
        self.client, self.db, self.collection = self.con_db(collection_str)

    # build dataframe for uniprot_swiss for loading into mongodb
    def load_uniprot(self, query=False, msg='', species=None):
        """Build dataframe
        
        Args:
            query (:obj:`bool`, optional): Whether download all reviewed entries of perform individual queries. Defaults to False.
            msg (:obj:`str`, optional): Query message. Defaults to ''.
            species (:obj:`list`, optional): species information to extract from df and loaded into uniprot. Defaults to None.
        """
        fields = '&columns=id,entry name,genes(PREFERRED),protein names,sequence,length,mass,ec,database(GeneID),reviewed,organism-id,database(KO),genes(ALTERNATIVE),genes(ORF),genes(OLN),database(EMBL),database(RefSeq),database(KEGG)'
        if not query:
            url = self.url + fields
        else:
            query_msg = msg
            if isinstance(species, list):
                for specie in species:
                    query_msg += '+'+str(specie)
            url = self.query_url + query_msg + '&sort=score' + fields
        url += '&format=tab'
        url += '&compress=no'
        if not math.isnan(self.max_entries):
           url += '&limit={}'.format(self.max_entries)
        
        try:
            response = requests.get(url, stream=False)
            response.raise_for_status()
        except requests.exceptions.ConnectionError:
            pass 
        

        try:
            data = pandas.read_csv(io.BytesIO(response.content), delimiter='\t', encoding='utf-8', low_memory=False)
        except pandas.errors.EmptyDataError:
            return
        except UnboundLocalError:
            return
        data.columns = [
            'uniprot_id', 'entry_name', 'gene_name', 'protein_name', 'canonical_sequence', 'length', 'mass',
            'ec_number', 'entrez_id', 'status', 'ncbi_taxonomy_id', 'ko_number', 'gene_name_alt',
            'gene_name_orf', 'gene_name_oln', 'sequence_embl', 'sequence_refseq', 'kegg_org_gene'
        ]
        data['entrez_id'] = data['entrez_id'].astype(str).str.replace(';', '')
        data['kegg_org_gene'] = data['kegg_org_gene'].astype(str).str.replace(';', '')

        try:
            data['mass'] = data['mass'].str.replace(',', '')
        except AttributeError:
            pass

        data['ko_number'] = data['ko_number'].astype(str).str.replace(';', '')
        data['gene_name_oln'] = data['gene_name_oln'].astype(str).str.split(' ')
        data['gene_name_orf'] = data['gene_name_orf'].astype(str).str.split(' ')
        data['gene_name_alt'] = data['gene_name_alt'].astype(str).str.split(' ')
        data['sequence_embl'] = self.embl_helper(data['sequence_embl'].astype(str).str)
        data['sequence_refseq'] = self.embl_helper(data['sequence_refseq'].astype(str).str)
        if species is None:
            self.load_df(data)
        else:
            self.load_df(data.loc[data['ncbi_taxonomy_id'].isin(species)])

    def embl_helper(self, s):
        """Processing emble or refseq strings into a list of standard
        format. "NP_796298.2 [E9PXF8-1];XP_006507989.1 [E9PXF8-2];" ->
        ['NP_796298.2', 'XP_006507989.1'].
        
        Args:
            s (:obj:`pandas.Dataframe`): object to be processed.

        Return:
            (:obj:`list` of :obj:`str`): list of processed strings
        """
        return s.replace('[', '').str.replace(';',' ').str.replace(']', '').str.split(' ').str[:-1]

    # load pandas.DataFrame into MongoDB
    def load_df(self, df):
        df_json = json.loads(df.to_json(orient='records'))
        try:        
            self.collection.insert(df_json)
        except pymongo.errors.InvalidOperation as e:
            return(str(e))

    def fill_species_name(self):
        ncbi_taxon_ids = self.collection.distinct('ncbi_taxonomy_id')
        start = 16650
        for i, _id in enumerate(ncbi_taxon_ids[start:]):
            if i == self.max_entries:
                break
            if i% 50 == 0 and self.verbose:
                print('Adding taxon name to {} out of {} records.'.format(i + start, len(ncbi_taxon_ids)))
            names = self.taxon_manager.get_name_by_id([_id])
            self.collection.update_many({'ncbi_taxonomy_id': _id},
                                        {'$set': {'species_name': names.get(_id)}},
                                        upsert=False)

    def fill_ko_name(self):
        ko_numbers = self.collection.distinct('ko_number')
        count = len(ko_numbers)
        start = 0
        for i, number in enumerate(ko_numbers[start:]):
            if i == self.max_entries:
                break
            if i % 50 == 0 and self.verbose:
                print('Adding ko name to {} out of {} records.'.format(i + start, count))
            kegg_name = self.ko_manager.get_def_by_kegg_id(number)
            self.collection.update_many({'ko_number': number},
                                        {'$set': {'ko_name': kegg_name}},upsert=False)

    def fill_species_info(self):
        """Fill ancestor information.
        """
        ncbi_taxon_ids = self.collection.distinct('ncbi_taxonomy_id')
        start = 0
        count = len(ncbi_taxon_ids)
        for i, _id in enumerate(ncbi_taxon_ids[start:]):
            if i == self.max_entries:
                break
            if i% 50 == 0 and self.verbose:
                print('Adding ancestor information to {} out of {} records.'.format(i + start, count))
            ancestor_taxon_id, ancestor_name = self.taxon_manager.get_anc_by_id([_id])
            self.collection.update_many({'ncbi_taxonomy_id': _id},
                                        {'$set': {'ancestor_taxon_id': ancestor_taxon_id[0],
                                                  'ancestor_name': ancestor_name[0]}},
                                        upsert=False)

    def load_abundance_from_pax(self):
        '''Load protein abundance data but interating from pax collection.
        '''
        _, _, col_pax = self.con_db('pax')
        query = {}
        projection = {'ncbi_id': 1, 'species_name': 1,
                    'observation': 1, 'organ': 1}
        docs = col_pax.find(filter=query, projection=projection, batch_size=5)
        count = col_pax.count_documents(query)
        progress = 285
        for i, doc in enumerate(docs[progress:]):            
            organ = doc['organ']
            if self.verbose and i % 1 == 0:
                print('Loading abundance info {} of {} ...'.format(
                    i + progress, min(count, self.max_entries)))
            for j, obs in enumerate(doc['observation']):
                if j == self.max_entries:
                    break
                if self.verbose and j % 100 == 0 and i % 1 == 0:
                    print('  Loading observation info {} of {} ...'.format(
                        j, len(doc['observation'])))
                try:
                    uniprot_id = obs['protein_id']['uniprot_id']
                    abundance = obs['abundance']
                    dic = {'organ': organ, 'abundance': abundance}
                    self.collection.update_one({'uniprot_id': uniprot_id},
                                        {'$push': {'abundances': dic} }, upsert=False,
                                        collation=self.collation)
                except TypeError:
                    continue

    def remove_redudant_id(self):
        """Remove redundant entries in uniprot collection.
        Priority:
            1. Has 'abundances' field
            1. Has 'ko_name' field
            2. Has 'kegg_org_gene' field
            3. Has 'orthologs' field
        """
        ids = self.collection.distinct('uniprot_id', collation=self.collation)
        projection = {'abundances': 1, 'ko_name': 1, 'kegg_org_gene': 1, 'orthologs': 1}
        for _id in ids:
            query = {'uniprot_id': _id}
            count = self.collection.count_documents(query, projection=projection, collation=self.collation)
            if count == 1:
                continue
            else:
                docs = self.collection.find(filter=query, projection=projection, collation=self.collation)
                to_be_removed = []
                to_be_saved = []
                for doc in docs:
                    pass

    def fill_reactions(self, start=0):
        """Fill reactions in which the protein acts as a catalyst.
        
        Args:
            start (:obj:`int`, optional): Starting document in sabio_rk. Defaults to 0.
        """
        docs = self.sabio_manager.collection.find({}, projection={'enzyme': 1, 'kinlaw_id': 1})
        count = self.sabio_manager.collection.count_documents({})
        for i, doc in enumerate(docs[start:]):
            if i == self.max_entries:
                break
            if i % 100 == 0 and self.verbose:
                print("Processing reaction {} out of {}...".format(i+start, count))
            enzyme = doc['enzyme']
            kinlaw_id = doc['kinlaw_id']
            if len(enzyme) == 0 or enzyme is None or enzyme[0].get('subunits') is None:
                continue
            else:
                enzyme_subunits = enzyme[0]['subunits']
                for subunit in enzyme_subunits:
                    if subunit is None:
                        continue
                    else:
                        uniprot_id = subunit['uniprot']
                        if i % 100 == 0 and self.verbose:
                            print("     Adding kinlaw_id {} into uniprot collection.".format(kinlaw_id))
                        self.collection.update_one({'uniprot_id': uniprot_id},
                                                {'$addToSet': {'sabio_kinlaw_id': kinlaw_id}},
                                                collation=self.collation, upsert=False)                
            

from multiprocessing import Pool, Process

def main():
    db = 'datanator'
    collection_str = 'uniprot'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    manager=UniprotNoSQL(MongoDB = server, db = db, 
    username = username, password = password, collection_str=collection_str,
    verbose=True)

    # manager.load_uniprot()

    # p = Process(target=manager.fill_species_name())
    # p.start()
    # p.join()

    # p = Process(target=manager.fill_ko_name())
    # p.start()
    # p.join()

    # p = Process(target=manager.fill_species_info())
    # p.start()
    # p.join()

    # p = Process(target=manager.load_abundance_from_pax())
    # p.start()
    # p.join()

    p = Process(target=manager.fill_reactions())
    p.start()
    p.join()

if __name__ == '__main__':
    main()