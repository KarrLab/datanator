import json
import requests
import re
import os
import requests
from bs4 import BeautifulSoup
from datanator_query_python.util import mongo_util
# from datanator_query_python.query import query_taxon_tree
# from pymongo.collation import Collation, CollationStrength
import datanator.config.core


class KeggOrgCode(mongo_util.MongoUtil):

    def __init__(self, MongoDB, db, cache_dirname=None, replicaSet=None, verbose=False, max_entries=float('inf'),
                username=None, password=None, readPreference=None, authSource='admin', collection_str='kegg_organism_code'):
        super().__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, verbose=verbose, max_entries=max_entries,
                        db=db, username=username, password=password, authSource=authSource, readPreference=readPreference)
        self.ENDPOINT_DOMAINS = {
            'root': 'https://www.genome.jp/kegg/catalog/org_list.html',
            'species': 'https://www.genome.jp/kegg/catalog/org_list4.html',
            'ncbi_lookup': 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name='
        }
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.collection_str = collection_str
        r = requests.get(self.ENDPOINT_DOMAINS['root'])
        self.soups = BeautifulSoup(r.content, 'html.parser')
        r = requests.get(self.ENDPOINT_DOMAINS['species'])
        self.species_soups = BeautifulSoup(r.content, 'html.parser')
        self.client, self.db, self.collection = self.con_db(self.collection_str)
        # self.collation = Collation(locale='en', strength=CollationStrength.SECONDARY)
        # self.taxon_manager = query_taxon_tree.QueryTaxonTree(collection_str='taxon_tree', 
        #         verbose=verbose, max_entries=max_entries, username=username, MongoDB=MongoDB, 
        #         password=password, db='datanator', authSource=authSource, readPreference=readPreference)

    def get_ncbi_id(self, name):
        """Given name of species, look up ncbi_taxonomy_id
        from official ncbi database by parsing html webpage.
        
        Args:
            name (:obj:`str`): name of the organism.

        Return:
            (:obj:`int`): NCBI Taxonomy ID.
        """
        endpoint = self.ENDPOINT_DOMAINS['ncbi_lookup'] + name
        r = requests.get(endpoint)
        soup = BeautifulSoup(r.content, 'html.parser')
        result = soup.find(string=re.compile('Taxonomy ID: '))
        if result is None:
            return result
        else:
            return int(str(result).split(': ')[1])
    
    def has_href_and_id(self, tag):
        return tag.has_attr('href') and tag.has_attr('id')

    def has_href_but_no_id(self, tag):
        return tag.has_attr('href') and not tag.has_attr('id')

    def parse_ids(self):
        """Parse HTML to get kegg organism codes.
        """
        for soup in self.soups.find_all(self.has_href_and_id):
            yield soup.get('id')

    def parse_names(self):
        """Parse HTML to get kegg organism names.
        """
        not_in = ['prag']
        for soup in self.soups.find_all(self.has_href_but_no_id):
            if soup.get('href') == '/dbget-bin/www_bfind?T00544': # special case
                yield 'Haemophilus influenzae PittGG (nontypeable)'
            result = re.search('.>(.*)<\/a>', str(soup))
            if result is not None and soup.get('href').startswith('/dbget-bin'):
                if result.group(1) not in not_in:
                    yield result.group(1)
                else:
                    continue
            else:
                continue

    def make_bulk(self, offset=0, bulk_size=100):
        """Make bulk objects to be inserted into MongoDB.

        Args:
            offset(:obj:`int`): Position of beginning (zero-indexed). Defaults to 0.
            bulk_size(:obj:`int`): number of objects. Defaults to 100.

        Return:
            (:obj:`list` of :obj:`dict`): list of objects to be inserted.
        """
        ids = self.parse_ids()
        names = self.parse_names()
        ncbi_ids = self.parse_species_id()
        species_names = self.parse_species_name()
        count = 0
        result = []
        for i, (_id, name, ncbi_id, species_name) in enumerate(zip(ids, names, ncbi_ids, species_names)):
            if count == bulk_size:
                break
            if i < offset:
                continue
            if count < bulk_size:
                if species_name == 'Lokiarchaeum sp. GC14 75':
                    result.append({"kegg_organism_id": _id, "org_name": [name, species_name],
                                    'ncbi_taxonomy_id': 1538547})
                elif species_name == 'archaeon GW2011 AR10':
                    result.append({"kegg_organism_id": _id, "org_name": [name, species_name],
                                    'ncbi_taxonomy_id': 1579370})
                elif species_name == 'archaeon GW2011 AR10':
                    result.append({"kegg_organism_id": _id, "org_name": [name, species_name],
                                    'ncbi_taxonomy_id': 1579378})                     
                else:   
                    result.append({"kegg_organism_id": _id, "org_name": [name, species_name],
                                    'ncbi_taxonomy_id': ncbi_id})
                count += 1
        return result

    def parse_species_name(self):
        """Parse species name from html.
        """
        not_in = ['Pragia sp. CF-458', 'Citrobacter sp. FDAARGOS 156', 'Psychrobacter sp. DAB AL43B',
                'Shewanella sp. FDAARGOS 354', 'Acidiferrobacter sp. SPIII 3', 'Janthinobacterium sp. 1 2014MBL_MicDiv',
                'Sulfurimonas sp. GYSZ 1', 'Halobacteriovorax sp. BALOs 7', 'Caulobacteraceae bacterium OTSz A_272',
                'Roseomonas sp. FDAARGOS 362', 'Streptomyces sp. CCM MD2014', 'Curtobacterium sp. MR MD2014',
                'Actinomyces sp. VUL4 3', 'Thermus sp. CCB US3_UF1', 'Capnocytophaga sp. FDAARGOS 737',
                'Maribacter sp. 1 2014MBL_MicDiv', 'Formosa sp. Hel1 33_131', 'Formosa sp. Hel3 A1_48',
                'Oceanihabitans sp. IOP 32', 'Candidatus Saccharibacteria bacterium RAAC3 TM7_1',
                'Candidatus Saccharibacteria bacterium GW2011 GWC2_44_17', 'candidate division TM6 bacterium GW2011 GWF2_28_16',
                'candidate division SR1 bacterium RAAC1 SR1_1', 'candidate division SR1 bacterium Aalborg AAW-1',
                'candidate division WWE3 bacterium RAAC2 WWE3_1', 'candidate division Kazan bacterium GW2011 GWA1_50_15',
                'Berkelbacteria bacterium GW2011 GWE1_39_12', 'Candidatus Beckwithbacteria bacterium GW2011 GWC1_49_16',
                'Candidatus Woesebacteria bacterium GW2011 GWF1_31_35', 'Candidatus Wolfebacteria bacterium GW2011 GWB1_47_1',
                'Candidatus Campbellbacteria bacterium GW2011 OD1_34_28', 'Lokiarchaeum sp. GC14 75',
                'archaeon GW2011 AR10', 'archaeon GW2011 AR20'] # no ncbi id given in html
        soups = self.species_soups.find_all(href=re.compile('/kegg-bin/show_organism'))
        for soup in soups:
            result = re.search('.>(.*)<\/a>', str(soup))
            name = result.group(1)
            if name not in not_in:
                yield name
            else:
                continue

    def parse_species_id(self):
        """Parse NCBI Taxonomy ID from html.
        """
        not_in = [2498113]
        soups = self.species_soups.find_all(href=re.compile('https://www.ncbi.nlm.nih.gov/Taxonomy/'))
        for soup in soups:
            result = re.search('.>(.*)<\/a>', str(soup))
            _id = int(result.group(1))
            if _id not in not_in and _id is not None:
                yield _id
            else:
                continue 
    
    def parse_species_html(self, obj):
        """Parse html (https://www.genome.jp/kegg/catalog/org_list4.html)
        to get NCBI Taxonomy ID. This method works
        because organisms in the source page (https://www.genome.jp/kegg/catalog/org_list4.html)
        is in the same order as (https://www.genome.jp/kegg/catalog/org_list.html)

        Args:
            obj(:obj:`list` of :obj)
        """
        pass

    def bulk_load(self, bulk_size=100):
        """Loading bulk data into MongoDB.
        
        Args:
            bulk_size(:obj:`int`): number of entries per insertion. Defaults to 100.
        """
        length = bulk_size
        count = 0
        while length != 0:
            if count == self.max_entries:
                break
            if count % 10 == 0 and self.verbose:
                print('Inserting bulk {} of {}'.format(count, bulk_size))            
            docs = self.make_bulk(offset=count * bulk_size, bulk_size=bulk_size)
            length = len(docs)
            if length != 0:
                self.collection.insert_many(docs)            
            count += 1

    # def fill_ncbi_id(self):
    #     """Fill collection with ncbi_taxonomy_id.
    #     """
    #     query = {}
    #     docs = self.collection.find(query)
    #     count = self.collection.count_documents(query)
    #     for i, doc in enumerate(docs):
    #         if i == self.max_entries:
    #             break
    #         if i % 50 == 0 and self.verbose:
    #             print('Processing doc {} out of {}.'.format(i, count))
    #         name = doc['org_name']
    #         ids = self.taxon_manager.get_ids_by_name(name)
    #         if len(ids) > 1:
    #             self.collection.update_one({'org_name': name},
    #                                         {'$set': {'ncbi_taxonomy_id': ids,
    #                                                   'ambiguous': True}}, upsert=False)
    #         else:
    #             self.collection.update_one({'org_name': name},
    #                                         {'$set': {'ncbi_taxonomy_id': ids,
    #                                                   'ambiguous': False}}, upsert=False)


def main():
    db = 'datanator'
    username = datanator.config.core.get_config()['datanator']['mongodb']['user']
    password = datanator.config.core.get_config()['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
    src = KeggOrgCode(MongoDB, db, username=username, password=password,
                    readPreference='nearest', authSource='admin', verbose=True)
    src.bulk_load()

if __name__ == '__main__':
    main()