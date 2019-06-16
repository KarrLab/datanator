import json
import requests
import os
from datanator.util import mongo_util
from datanator.util import file_util


class KeggOrthology(mongo_util.MongoUtil):

    def __init__(self, cache_dirname, MongoDB, db, replicaSet=None, 
        verbose=False, max_entries=float('inf'), username = None,
        password = None, authSource = 'admin'):
        self.ENDPOINT_DOMAINS = {
            'root': 'https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=',
        }
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.collection = 'kegg_orthology'
        self.path = os.path.join(self.cache_dirname, self.collection)
        super(KeggOrthology, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, db=db,
                                            verbose=verbose, max_entries=max_entries, username = username,
                                            password = password, authSource = authSource)
        self.file_manager = file_util.FileUtil()

    def load_content(self):
        '''Load kegg_orthologs into MongoDB
        '''
        _, _, collection = self.con_db(self.collection)
        root_url = self.ENDPOINT_DOMAINS['root']
        if self.verbose:
            print('\n Downloading root kegg orthology file ...')
        manager = requests.get(root_url)
        manager.raise_for_status()
        os.makedirs(self.path, exist_ok=True)
        file_name = manager.json()['name']
        store_path = os.path.join(self.path, file_name)
        data = manager.json()
        with open(store_path, 'w') as f:
            json.dump(data, f, indent=4)

        names = self.file_manager.extract_values(data, 'name')
        names = [name.split()[0] for name in names if name[0] == 'K']
        names = list(set(names)) # remove duplicate name
        names.sort()

        iterations = min(len(names), self.max_entries)

        file_format = '.txt'

        i = 0
        for name in names:
            if i == self.max_entries:
                break
            if self.verbose and i % 100 == 0:
                print('Downloading {} of {} kegg orthology file {}...'.format(
                    i, iterations, name))
            self.download_ko(name+file_format)
            doc = self.parse_ko_txt(name+file_format)
            if doc is not None:
                collection.replace_one(
                    {'kegg_orthology_id': doc['kegg_orthology_id']}, doc, upsert=True)

            i += 1

    def parse_definition(self, line):
        '''Definition line could be something as follows:
         " fructose-bisphosphate aldolase / 6-deoxy-5-ketofructose 1-phosphate synthase [NADP...] [EC:4.1.2.13 2.2.1.11]\n"
         EC code can be optional
        '''
        if len(line) != 0:
            head_tail_clean = line.strip().replace('\n', '')
            head_tail_clean = head_tail_clean[11:]
            sep_name = head_tail_clean.split('/')  # list

            name_sans_last = sep_name[:-1]  # list
            last_name = sep_name[-1].split('[')[0].strip()  # string
            name_list = [item.strip() for item in name_sans_last] + [last_name]

            if len(sep_name[-1].split('[EC')) > 1: # EC code exists
                ec_codes = sep_name[-1].split('[EC')[1]  # string: 'EC:1.1.1.1 2.2.2.2]'
                strip_ec = ec_codes.split(':')[-1]
                ec_list = strip_ec.replace(']', '').split(' ')
            else:
                ec_list = []

            return (name_list, ec_list)
        else:
            pass

    def parse_ko_txt(self, filename):
        '''Parse kegg_ortho txt file into dictionary object
        '''
        file_path = os.path.join(self.path, filename)
        try: 
            with open(file_path, 'r') as f:
                doc = {}
                lines = f.readlines()

                # get entry ID
                doc['kegg_orthology_id'] = lines[0].split()[1]
                # get list of gene name
                doc['gene_name'] = [name.replace(
                    ',', '') for name in lines[1].split()[1:]]
                # get definition
                enzyme_name, ec = self.parse_definition(lines[2])    
                doc['definition'] = {'name': enzyme_name, 'ec_code': ec}

                # get first word of all the lines
                first_word = [line.split()[0] for line in lines]

                # find line number of certain first words of interest
                WOI_nonrepeating = ['MODULE', 'BRITE', 'GENES']
                reference = 'REFERENCE'

                index_nonrepeating = [first_word.index(
                    word) if word in first_word else -1 for word in WOI_nonrepeating]
                index_reference = [i for i, v in enumerate(
                    first_word) if v == reference]

                module_begin = index_nonrepeating[0]
                module_end = index_nonrepeating[1]
                gene_begin = index_nonrepeating[2]
                if len(index_reference) == 0:
                    gene_end = len(lines) - 1
                else:
                    gene_end = index_reference[0]

                # get module's key,value pairs
                if module_begin == -1:
                    doc['kegg_module_id'] = None
                elif (module_end - module_begin) == 1:  # only 1 module
                    doc['kegg_module_id'] = lines[module_begin].split()[1]
                else:
                    module_list = []
                    module_list.append(lines[module_begin].split()[1])
                    module_list += [line.split()[0]
                                    for line in lines[module_begin+1:module_end]]
                    doc['kegg_module_id'] = module_list

                # get gene_id's key,value pairs
                if (gene_end - gene_begin) == 1:  # only one ortholog
                    doc['gene_ortholog'] = {'organism': lines[gene_begin].split()[1].replace(':', ''),
                                            'gene_id': lines[gene_begin].split()[2]}
                else:
                    ortholog_list = []
                    ortholog_list.append({'organism': lines[gene_begin].split()[1].replace(':', ''),
                                          'gene_id': lines[gene_begin].split()[2]})
                    organisms = [line.split()[0].replace(':', '')
                                 for line in lines[gene_begin+1:gene_end]]
                    gene_ids = [line.split()[1:]
                                for line in lines[gene_begin+1:gene_end]]
                    for organism, gene_id in zip(organisms, gene_ids):
                        ortholog_list += [{
                            'organism': organism, 'gene_id': gene_id}]
                    doc['gene_ortholog'] = ortholog_list

                # get reference's namespace:value pairs
                ref_list = []
                reference_line = [lines[i] for i in index_reference]
                try:
                    reference_info = [line.split()[1] for line in reference_line]
                    for info in reference_info:
                        ref_list.append({'namespace': info.split(':')[0],
                                         'id': info.split(':')[1]})
                except IndexError:
                    pass    

                doc['reference'] = ref_list
                return doc

        except FileNotFoundError as e:
            log_file = os.path.join(self.path, 'kegg_orthology_log.txt')
            with open(log_file, 'a') as f:
                f.write(str(e)+'\n')
            pass

            

    def download_ko(self, name):
        address = name.split('.')[0]
        try:
            info = requests.get("http://rest.kegg.jp/get/ko:{}".format(address))
            info.raise_for_status()
            file_name = os.path.join(self.path, name)
            with open(file_name, 'w') as f:
                f.write(info.text)
        except requests.exceptions.HTTPError as e:
            log_file = os.path.join(self.path, 'kegg_orthology_log.txt')
            with open(log_file, 'a') as f:
                f.write(str(e) + '\n')
            pass
