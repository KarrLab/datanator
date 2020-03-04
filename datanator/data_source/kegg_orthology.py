import json
import requests
import os
import re
from datanator.util import mongo_util
from datanator.util import file_util
import datanator.config.core
import tempfile


class KeggOrthology(mongo_util.MongoUtil):

    def __init__(self, cache_dirname=None, MongoDB=None, db=None, replicaSet='', 
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

        for i, name in enumerate(names):
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

    def parse_pathway_disease(self, lines, category='pathway'):
        """Parse parthway or disease or module information
        
        Args:
            line (:obj:`readlines()`): pathway lines.
            category (:obj:`str`): which category to parse. Defaults to pathway.

        Return:
            (:obj:`list` of :obj:`dict`): list of pathways [{"kegg_pathway_code": ..., "pathway_description": ...}]
        """
        result = []
        if category == 'pathway':
            key_0 = "kegg_pathway_code"
            key_1 = "pathway_description"
        elif category == 'disease':
            key_0 = "kegg_disease_code"
            key_1 = "kegg_disease_name"
        elif category == 'module':
            key_0 = "kegg_module_code"
            key_1 = "kegg_module_name"            

        for i, line in enumerate(lines):
            if i == 0:
                _list = line.split()
                kegg_pathway = _list[1]
                pathway_description = ' '.join(_list[2:])
                result.append({key_0: kegg_pathway, key_1: pathway_description})
            else:
                _list = line.split()
                kegg_pathway = _list[0]
                pathway_description = ' '.join(_list[1:])
                result.append({key_0: kegg_pathway, key_1: pathway_description})
        return result

    def parse_gene(self, lines):
        """Parse GENES category
        (http://rest.kegg.jp/get/ko:K00023)
        
        Args:
            lines (:obj:`readlines()`): Lines for genes.

        Return:
            (:obj:`list` of :obj:`dict`): list of parsed genes.
        """
        def parse_annotation(gene_annotation):
            """Parse 'XCC2294(phbB)' into structured object
            
            Args:
                gene_annotation (:obj:`str`): annotation string.

            Return:
                (:obj:`dict`): {'locus_id': locus, 'gene_id': gene_id}
            """
            locus = gene_annotation.split('(')[0]
            gene_id = re.search(r'\((.*?)\)', gene_annotation)
            if gene_id is not None:
                gene_id = gene_id.group(1)
            return {'locus_id': locus, 'gene_id': gene_id}

        result = []
        for i, line in enumerate(lines):
            if i == 0:
                _list = line.split()
                org_code = _list[1][:-1]
                genetic_info = []
                for gene_annotation in _list[2:]:
                    genetic_info.append(parse_annotation(gene_annotation))
                result.append({'organism': org_code, 'genetic_info': genetic_info})
            else:
                _list = line.split()
                org_code = _list[0][:-1]
                genetic_info = []
                for gene_annotation in _list[1:]:
                    genetic_info.append(parse_annotation(gene_annotation))
                result.append({'organism': org_code, 'genetic_info': genetic_info})
        return result                

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
                WOI_nonrepeating = ['PATHWAY', 'MODULE', 'DISEASE', 'BRITE', 'GENES']
                reference = 'REFERENCE'

                index_nonrepeating = [first_word.index(
                    word) if word in first_word else -1 for word in WOI_nonrepeating]
                index_reference = [i for i, v in enumerate(
                    first_word) if v == reference]

                pathway_begin = index_nonrepeating[0]
                module_begin = index_nonrepeating[1]
                disease_begin = index_nonrepeating[2]
                brite_begin = index_nonrepeating[3]
                gene_begin = index_nonrepeating[4]
                
                if pathway_begin == -1:
                    doc['kegg_pathway'] = None
                else:
                    pathway_end = list(filter(lambda i: i != -1, index_nonrepeating))[1]
                    pathway_lines = lines[pathway_begin:pathway_end] 
                    doc['kegg_pathway'] = self.parse_pathway_disease(pathway_lines)
                
                if len(index_reference) == 0:
                    gene_end = len(lines) - 1
                else:
                    gene_end = index_reference[0]

                if disease_begin == -1:
                    doc['kegg_disease'] = None
                    module_end = brite_begin
                else:
                    module_end = disease_begin
                    disease_lines = lines[disease_begin:brite_begin]
                    doc['kegg_disease'] = self.parse_pathway_disease(disease_lines, category='disease')

                if module_begin == -1:
                    doc['kegg_module'] = None
                else:
                    module_lines = lines[module_begin:module_end]
                    doc['kegg_module'] = self.parse_pathway_disease(module_lines, category='module')

                # get gene_id's key,value pairs
                doc['gene_ortholog'] = self.parse_gene(lines[gene_begin:gene_end])

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
            with open(file_name, 'w+') as f:
                f.write(info.text)
        except requests.exceptions.HTTPError as e:
            log_file = os.path.join(self.path, 'kegg_orthology_log.txt')
            with open(log_file, 'a') as f:
                f.write(str(e) + '\n')
            pass

def main():
    db = 'datanator'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    cache_dirname = tempfile.mkdtemp()
    manager = KeggOrthology(MongoDB=MongoDB,  db=db, cache_dirname=cache_dirname,
                            verbose=True, username=username,
                            password=password)
    manager.load_content()

if __name__ == '__main__':
    main()