import json
import requests
import os
from six import BytesIO
from datanator.util import mongo_util
import zipfile

class TaxonTree(mongo_util.MongoUtil):

    def __init__(self, cache_dirname, MongoDB, db, replicaSet=None, verbose=False, max_entries=float('inf')):
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
                                            verbose=verbose, max_entries=max_entries)
        self.client, self.db, self.collection = self.con_db(self.collection_str)

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
        tax_id = a[0]
        tax_name = a[1]
        something =  [item.split(';') for item in a] 
        ancestor_name = [elem.lstrip() for elem in something[2][:-1]]
        return [tax_id, tax_name, ancestor_name]
    
    def parse_taxid_line(self, line):
        '''Parses lines in file taxidlineage.dmp and return elements in a list
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

        '''Parse fullnamelineage.dmp, store in MongoDB
            always runs ahead of other file parsing modules
        '''
        full_name = os.path.join(self.path, 'fullnamelineage.dmp')
        tax_id = os.path.join(self.path, 'taxidlineage.dmp')
        i = 0
        with open(full_name, 'r') as f1, open(tax_id, 'r') as f2:
            count = min(self.max_entries, self.count_line(full_name))
            for line_name, line_id in zip(f1, f2):
                if i == self.max_entries:
                    break
                if self.verbose and i % 10 == 0:
                    print ('Parsing lineage line {} of {}...'.format(i+1, count))

                lineage_dict = {}
                elem_name = self.parse_fullname_line(line_name)
                elem_id = self.parse_taxid_line(line_id)
                lineage_dict['tax_id'] = elem_name[0]
                lineage_dict['tax_name'] = elem_name[1]
                lineage_dict['anc_name'] = elem_name[2]
                lineage_dict['anc_id'] = elem_id

                self.collection.update_one( {'tax_id': lineage_dict['tax_id']},
                                            {'$set': lineage_dict},
                                            upsert = True
                                            )

                i += 1

    def parse_nodes_line(self, line):
        '''Parse lines in nodes.dmp
        '''
        return [item.replace('\t', '') for item in line.split('|')[:-1]]

    def parse_nodes(self):
        file_name = os.path.join(self.path, 'nodes.dmp')
        i = 0
        count = min(self.max_entries, self.count_line(file_name))
        with open(file_name, 'r') as f:
            for line in f:
                if i == self.max_entries:
                    break
                if self.verbose and i%10==0:
                    print ('Parsing nodes line {} of {} ...'.format(i+1, count))
                node_dict = {}
                elem = self.parse_nodes_line(line)
                node_dict['tax_id'] = elem[0]
                node_dict['rank'] = elem[2]
                node_dict['locus_name_prefix'] = elem[3]
                node_dict['division_id'] = elem[4]
                node_dict['gene_code'] = elem[6]
                node_dict['comments'] = elem[-6]
                node_dict['plastid_gene_code'] = elem[-5]
                node_dict['hydrogenosome_gene_id'] = elem[-2]

                self.collection.update_one( {'tax_id': node_dict['tax_id']},
                                            {'$set': node_dict},
                                            upsert = True
                                            )
                i += 1

    def parse_division(self):
        file_name = os.path.join(self.path, 'division.dmp')
        i = 0
        count = min(self.max_entries, self.count_line(file_name))
        with open(file_name, 'r') as f:
            for line in f:
                if i == self.max_entries:
                    break
                if self.verbose and i % 10 == 0:
                    print ('Parsing division line {} of {} ...'.format(i+1, count))
                name_dict = {}
                elem = self.parse_nodes_line(line)
                name_dict['division_id'] = elem[0]
                name_dict['division_cde'] = elem[1]
                name_dict['division_name'] = elem[2]
                name_dict['division_comments'] = elem[3]

                self.collection.update_one( {'division_id': name_dict['division_id']},
                                            {'$set': name_dict},
                                            upsert = True
                                            )
                i += 1

    def parse_names(self):
        file_name = os.path.join(self.path, 'names.dmp')
        i = 0
        count = min(self.max_entries, self.count_line(file_name))
        with open(file_name, 'r') as f:
            for line in f:
                if i == self.max_entries:
                    break
                if self.verbose and i % 10 == 0:
                    print ('Parsing names line {} of {} ...'.format(i+1, count))
                name_dict = {}
                elem = self.parse_nodes_line(line)
                name_dict['tax_id'] = elem[0]
                name_dict['name_txt'] = elem[1]
                name_dict['unique_variant_name'] = elem[2]
                name_dict['name_class'] = elem[3]

                self.collection.update_one( {'tax_id': name_dict['tax_id']},
                                            {'$set': name_dict},
                                            upsert = True
                                            )
                i += 1

    def parse_gencode(self):
        file_name = os.path.join(self.path, 'gencode.dmp')
        i = 0
        count = min(self.max_entries, self.count_line(file_name))
        with open(file_name, 'r') as f:
            for line in f:
                if i == self.max_entries:
                    break
                if self.verbose and i % 10 == 0:
                    print ('Parsing gencode line {} of {} ...'.format(i+1, count))
                gencode_dict = {}
                elem = self.parse_nodes_line(line)
                gencode_dict['gene_code'] = elem[0]
                gencode_dict['abbreviation'] = elem[1]
                gencode_dict['gene_code_name'] = elem[2]
                gencode_dict['gene_code_cde'] = elem[3]
                gencode_dict['gene_code_starts'] = elem[4]

                self.collection.update_one( {'gene_code': gencode_dict['gene_code']},
                                            {'$set': gencode_dict},
                                            upsert = True
                                            )
                i += 1



