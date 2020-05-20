from io import BytesIO
import json
import os
import pymongo
import requests
from collections import deque
import zipfile
from datanator.util import mongo_util
import datanator.config.core


class TaxonTree(mongo_util.MongoUtil):

    def __init__(self, cache_dirname, MongoDB, db, replicaSet=None, 
        verbose=False, max_entries=float('inf'), username = None,
        password = None, authSource = 'admin'):
        self.ENDPOINT_DOMAINS = {
            'root': 'https://ftp.ncbi.nlm.nih.gov',
        }
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.collection_str = 'taxon_tree'
        self.path = os.path.join(self.cache_dirname, self.collection_str)
        super(TaxonTree, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, db=db,
                                            verbose=verbose, max_entries=max_entries, username = username,
                                            password = password, authSource = authSource)
        self.client, self.db, self.collection = self.con_db(self.collection_str)
        self.repetition = 1 # how often verbose messages show

    def load_content(self):
        '''Load contents of several .dmp files into MongoDB
        '''
        self.download_dump()
        self.parse_fullname_taxid() # taxidlineage.dmp fullnamelineage.dmp
        if self.verbose:
            print('Indexing tax_id ... \n')
        self.collection.create_index( [("tax_id", pymongo.ASCENDING)] , background=False, sparse=True)
        self.parse_nodes() # nodes.dmp
        if self.verbose:
            print('Indexing division_id and gene_code ... \n')
        index1 = pymongo.IndexModel( [("division_id", pymongo.ASCENDING)] , background=False, sparse=True)
        index2 = pymongo.IndexModel([("gene_code", pymongo.ASCENDING)] , background=False, sparse=True)
        self.collection.create_indexes([index1, index2])
        self.parse_division() # division.dmp
        self.parse_names() # names.dmp
        self.parse_gencode() # gencode.dmp


    def download_dump(self):

        os.makedirs(self.path, exist_ok=True)
        cwd = '/pub/taxonomy/new_taxdump/'
        noi = 'new_taxdump.zip'
        database_url = self.ENDPOINT_DOMAINS['root'] + cwd + noi
        local_filename = os.path.join(self.path, noi)
        if self.verbose:
            print ('\n Downloading taxdump zip file ...')
        response = requests.get(database_url)
        response.raise_for_status()

        if self.verbose:
            print (' ... Done!')
            print ('Unzipping ...')
        z = zipfile.ZipFile(BytesIO(response.content))
        z.extractall(self.path)

        if self.verbose:
            print('... Done unzipping')

    def parse_fullname_line(self, line):
        '''Parses lines in file fullnamelineage.dmp and return elements in a list
        '''
        a =  [item.replace('\t', '') for item in line.split('|')[:-1]]
        tax_id = a[0].strip()
        tax_name = a[1].strip()
        something =  [item.split(';') for item in a] 
        ancestor_name = [elem.lstrip() for elem in something[2][:-1]]
        return [tax_id, tax_name, ancestor_name]
    
    def parse_taxid_line(self, line):
        '''Parses lines in file taxidlineage.dmp and return elements in a list
        delimited by "\t|\n"
        (tab, vertical bar, and newline) characters. Each record consists of one 
        or more fields delimited by "\t|\t" (tab, vertical bar, and tab) characters.
        '''
        a =  [item.replace('\t', '') for item in line.split('|')[:-1]]
        return a[1].split(' ')[:-1]      

    def count_line(self, file):
        '''Efficiently count total number of lines in a given file
        '''
        with open(file) as f:
            for i, l in enumerate(f):
                pass
        return i + 1


    def parse_fullname_taxid(self):
        '''Parse fullnamelineage.dmp and taxidlineage.dmp store in MongoDB
           Always run first before loading anything else
           (insert_one)
        '''
        full_name = os.path.join(self.path, 'fullnamelineage.dmp')
        tax_id = os.path.join(self.path, 'taxidlineage.dmp')
        i = 0
        with open(full_name, 'r') as f1, open(tax_id, 'r') as f2:
            count = min(self.max_entries, self.count_line(full_name))
            for line_name, line_id in zip(f1, f2):
                if i == self.max_entries:
                    break
                if self.verbose and i % self.repetition == 0:
                    print ('Parsing lineage line {} of {}...'.format(i+1, count))

                lineage_dict = {}
                elem_name = self.parse_fullname_line(line_name)
                elem_id = self.parse_taxid_line(line_id)
                lineage_dict['tax_id'] = int(elem_name[0])
                lineage_dict['tax_name'] = elem_name[1]
                lineage_dict['anc_name'] = elem_name[2]
                lineage_dict['anc_id'] = [int(item) for item in elem_id]

                self.collection.insert_one( lineage_dict
                                            )

                i += 1

    def parse_nodes_line(self, line):
        '''Parse lines in nodes.dmp
        '''
        return [item.replace('\t', '') for item in line.split('|')[:-1]]

    def parse_nodes(self):
        '''nodes.dmp
        '''
        file_name = os.path.join(self.path, 'nodes.dmp')
        i = 0
        count = min(self.max_entries, self.count_line(file_name))
        with open(file_name, 'r') as f:
            for line in f:
                if i == self.max_entries:
                    break
                if self.verbose and i%self.repetition==0:
                    print ('Parsing nodes line {} of {} ...'.format(i+1, count))
                node_dict = {}
                elem = self.parse_nodes_line(line)
                tax_id = int(elem[0])
                node_dict['rank'] = elem[2]
                node_dict['locus_name_prefix'] = elem[3]
                node_dict['division_id'] = int(elem[4])
                node_dict['gene_code'] = int(elem[6])
                node_dict['comments'] = elem[-6]
                node_dict['plastid_gene_code'] = elem[-5]
                node_dict['hydrogenosome_gene_id'] = elem[-2]

                self.collection.update_one( {'tax_id': tax_id},
                                            {'$set': node_dict},
                                            upsert = True
                                            )
                i += 1

    def parse_division(self):
        '''division.dmp
        '''
        file_name = os.path.join(self.path, 'division.dmp')
        i = 0
        count = min(self.max_entries, self.count_line(file_name))
        with open(file_name, 'r') as f:
            for line in f:
                if i == self.max_entries:
                    break
                if self.verbose and i % self.repetition == 0:
                    print ('Parsing division line {} of {} ...'.format(i+1, count))
                name_dict = {}
                elem = self.parse_nodes_line(line)
                division_id = int(elem[0])
                name_dict['division_cde'] = elem[1]
                name_dict['division_name'] = elem[2]
                name_dict['division_comments'] = elem[3]

                self.collection.update_many( {'division_id': division_id},
                                            {'$set': name_dict},
                                            upsert = True
                                            )
                i += 1

    def parse_names(self):
        '''names.dmp
        1   |   all |       |   synonym |
        1   |   root    |       |   scientific name |
        2   |   bacteria    |   bacteria <blast2>   |   blast name  |
        2   |   Bacteria    |   Bacteria <prokaryotes>  |   scientific name |
        2   |   eubacteria  |       |   genbank common name 
        '''
        file_name = os.path.join(self.path, 'names.dmp')
        i = 0
        count = min(self.max_entries, self.count_line(file_name))
        with open(file_name, 'r') as f:
            for line in f:
                if i == self.max_entries:
                    break
                if self.verbose and i % self.repetition == 0:
                    print ('Parsing names line {} of {} ...'.format(i+1, count))
                name_dict = {}
                elem = self.parse_nodes_line(line)
                tax_id = int(elem[0])
                name_dict['name_txt'] = elem[1]
                name_dict['unique_variant_name'] = elem[2]
                name_dict['name_class'] = elem[3]

                self.collection.update_one( {'tax_id': tax_id},
                                            {'$addToSet': name_dict},
                                            upsert = True
                                            )
                i += 1

    def parse_gencode(self):
        '''gencode.dmp
        '''
        file_name = os.path.join(self.path, 'gencode.dmp')
        i = 0
        count = min(self.max_entries, self.count_line(file_name))
        with open(file_name, 'r') as f:
            for line in f:
                if i == self.max_entries:
                    break
                if self.verbose and i % self.repetition == 0:
                    print ('Parsing gencode line {} of {} ...'.format(i+1, count))
                gencode_dict = {}
                elem = self.parse_nodes_line(line)
                gene_code = int(elem[0])
                gencode_dict['abbreviation'] = elem[1]
                gencode_dict['gene_code_name'] = elem[2]
                gencode_dict['gene_code_cde'] = elem[3]
                gencode_dict['gene_code_starts'] = elem[4]

                self.collection.update_many( {'gene_code': gene_code},
                                            {'$set': gencode_dict},
                                            upsert = True
                                            )
                i += 1

    def insert_canon_anc(self, start=0):
        """Insert two arrays to each document, one is
        canon_anc_id, the other is canon_anc_name
        """
        query = {"canon_anc_ids": {"$exists": False}}
        canon_info = {}  # store info of canon to avoid multiple db queries
        ids_ranked = deque()
        canons = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']
        count = self.collection.count_documents(query)
        projection = {'anc_name': 1, 'anc_id': 1, 'tax_id': 1}
        docs = self.collection.find(filter=query, skip=start, projection=projection,
                                    no_cursor_timeout=True, batch_size=1000).sort('tax_id', 1).hint('tax_id_1')
        for i, doc in enumerate(docs):
            if i == self.max_entries:
                break
            if self.verbose and i % 100 == 0:
                print('Processing doc {} out of {}'.format(i + start, count))
            canon_anc_names = deque()
            canon_anc_ids = deque()
            anc_ids = doc['anc_id']
            anc_names = doc['anc_name']
            if len(anc_ids) == 0:
                continue
            for anc_id, anc_name in zip(anc_ids, anc_names):
                if anc_id not in ids_ranked: # no need to look in canon_info dictionary
                    rank = self.collection.find_one({'tax_id': anc_id})['rank']
                    if rank in canons:
                        canon_anc_ids.append(anc_id)
                        canon_anc_names.append(anc_name)
                        canon_info[anc_id] = True
                        ids_ranked.append(anc_id)
                    else:
                        ids_ranked.append(anc_id)
                else: # no need to perform db lookups
                    c = canon_info.get(anc_id, False)
                    if c:
                        canon_anc_ids.append(anc_id)
                        canon_anc_names.append(anc_name)
                    else:
                        continue                                               
            # update doc
            self.collection.update_one({'_id': doc['_id']},
                                       {'$set': {'canon_anc_ids': list(canon_anc_ids),
                                                 'canon_anc_names': list(canon_anc_names)}},
                                       upsert=False)

def main():
    db = 'datanator'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    cache_dirname = '/root/karr_lab/datanator/datanator/data_source/cache/taxon_tree'
    manager = TaxonTree(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=None, 
                    db=db, verbose=True, username=username, password=password)
    # manager.load_content()
    # manager.insert_canon_anc(start=250600)
    manager.insert_canon_anc(start=0)


if __name__ == '__main__':
    main()

