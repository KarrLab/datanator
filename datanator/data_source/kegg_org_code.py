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
            suggestion = soup.find_all(attrs={"title": "species"})
            if suggestion == []:
                return None
            else:
                href = suggestion[0].get('href')
                id_list = re.search('.id=(.*)&lvl=', str(href))
                return int(id_list.group(1))
        else:
            return int(str(result).split(': ')[1])

    def has_align_but_no_rowspan(self, tag):
        return tag.has_attr('align') and not tag.has_attr('rowspan')

    def parse_html_iter(self):
        """Parse org code HTML iteratively.

        Yield:
            (:obj:`obj`): {'kegg_organism_id':  , 'org_name':   , 'org_synonym':  }
        """
        results = self.soups.find_all('tr', {'align': 'center'})
        for result in results:
            obj = {}
            code_row = result.find(self.has_align_but_no_rowspan, {'align': 'center'})
            if code_row is not None:
                obj['kegg_organism_id'] = str(code_row.get_text())
            name_row = result.find(self.has_align_but_no_rowspan, {'align': 'left'})
            if name_row is not None:
                name_str = str(name_row.get_text())
                name = re.search('(.*?)\(.*?\)', name_str)
                if name is not None:
                    obj['org_name'] = name.group(1)[:-1]
                    obj['org_synonym'] = re.search('.*?\((.*?)\)', name_str).group(1)
                else:
                    obj['org_name'] = name_str
                    obj['org_synonym'] = None
            yield obj

    def make_bulk(self, offset=0, bulk_size=100):
        """Make bulk objects to be inserted into MongoDB.

        Args:
            offset(:obj:`int`): Position of beginning (zero-indexed). Defaults to 0.
            bulk_size(:obj:`int`): number of objects. Defaults to 100.

        Return:
            (:obj:`list` of :obj:`dict`): list of objects to be inserted.
        """
        objs = self.parse_html_iter()
        count = 0
        skipped = 0
        result = []
        for i, obj in enumerate(objs):
            if i + skipped < offset:
                continue
            if obj == {}:
                skipped += 1
                continue
            if count == bulk_size:
                break
            if count < bulk_size:
                name = obj['org_name']
                ncbi_id = self.get_ncbi_id(name)
                if ncbi_id is None:
                    ncbi_id = self.get_ncbi_id_rest(name)
                obj['ncbi_taxonomy_id'] = ncbi_id   
                result.append(obj)
                count += 1
        return result

    def get_ncbi_id_rest(self, name):
        """Get ncbi taxonomy id of an organism using
        api.datanator.info
        
        Args:
            name (:obj:`str`): Name of the organism.

        Return:
            (:obj:`int`): NCBI Taxonomy ID.
        """
        endpoint = "https://api.datanator.info/ftx/text_search/num_of_index/?query_message={}&index=taxon_tree&from_=0&size=5&fields=tax_name&fields=name_txt".format(name)
        r = requests.get(endpoint)
        data = json.loads(r.text)
        if data.get('taxon_tree', []) !=[]:
            return data['taxon_tree'][0]['tax_id']
        else:
            return None

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
            if count % 1 == 0 and self.verbose:
                print('Inserting bulk {} ...'.format(count))            
            docs = self.make_bulk(offset=count * bulk_size, bulk_size=bulk_size)
            length = len(docs)
            if length != 0:
                self.collection.insert_many(docs)            
            count += 1


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