import csv
import ete3
import os
import zipfile
from six import BytesIO
import json
import pymongo
import requests

class CorumNoSQL():
    def __init__(self,cache_dirname, MongoDB, db, verbose=False, max_entries=float('inf')):
        self.ENDPOINT_DOMAINS = {
            'corum': 'https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip',
        }
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.collection = 'corum'

    def con_db(self):
        try:
            client = pymongo.MongoClient(self.MongoDB, 400)  # 400ms max timeout
            client.server_info()
            db = client[self.db]
            collection = db[self.collection]
            return collection
        except pymongo.errors.ConnectionFailure:
            return ('Server not available')

    def load_content(self):
        """ Collect and parse all data from CORUM website into JSON files and add to NoSQL database """
        database_url = self.ENDPOINT_DOMAINS['corum']
        os.makedirs(os.path.join(
            self.cache_dirname, 'corum'), exist_ok=True)

        if self.verbose:
            print('Download list of all compounds: ...')

        response = requests.get(database_url)
        response.raise_for_status()
        # Extract All Files and save to current directory

        if self.verbose:
            print('... Done!')
        if self.verbose:
            print('Unzipping and parsing compound list ...')

        z = zipfile.ZipFile(BytesIO(response.content))
        z.extractall(self.cache_dirname)
        cwd = os.path.join(self.cache_dirname, 'allComplexes.txt')

        # create object to find NCBI taxonomy IDs
        ncbi_taxa = ete3.NCBITaxa()

        with open(cwd, 'r') as file:
            i_entry = 0
            for entry in csv.DictReader(file, delimiter='\t'):
                # entry/line number in file
                i_entry += 1

                # stop if the maximum desired number of entries has been reached
                if i_entry > self.max_entries:
                    break

                # replace 'None' strings with None
                for key, val in entry.items():
                    if val == 'None':
                        entry[key] = None

                # extract attributes
                complex_id = int(entry['ComplexID'])
                entry['ComplexID'] = complex_id #replace string value with int value
                complex_name = entry['ComplexName']
                cell_line = entry['Cell line']
                pur_method = entry['Protein complex purification method']
                # SETS OF INT IDS SEPARATED BY ; eg. GO:0005634
                go_id = entry['GO ID']
                go_dsc = entry['GO description']
                funcat_id = entry['FunCat ID']
                funcat_dsc = entry['FunCat description']
                pubmed_id = int(entry['PubMed ID'])
                entry['PubMed ID'] = pubmed_id
                gene_name = entry['subunits(Gene name)']
                gene_syn = entry['subunits(Gene name syn)']
                complex_syn = entry['Synonyms']
                disease_cmt = entry['Disease comment']
                su_cmt = entry['Subunits comment']
                complex_cmt = entry['Complex comment']

                su_uniprot = entry['subunits(UniProt IDs)']  # SETS OF STRING IDS SEPARATED BY ;\
                su_entrez = entry['subunits(Entrez IDs)']  # SETS OF INT IDS SEPARATED BY ;
                protein_name = entry['subunits(Protein name)']
                swissprot_id = entry['SWISSPROT organism']

                """ ----------------- Apply field level corrections-----------------"""
                # Split the semicolon-separated lists of subunits into protein components,
                # ignoring semicolons inside square brackets
                su_uniprot_list = parse_list(su_uniprot)
                entry['subunits(UniProt IDs)'] = su_uniprot_list
                
                su_entrez_list = parse_list(su_entrez)
                entry['subunits(Entrez IDs)'] = su_entrez_list
                
                go_id_list = parse_list(go_id)
                entry['GO ID'] = go_id_list
                
                go_dsc_list = parse_list(go_dsc)
                entry['GO description'] = go_dsc_list

                funcat_id_list = parse_list(funcat_id)
                entry['FunCat ID'] = funcat_id_list

                funcat_dsc_list = parse_list(funcat_dsc)
                entry['FunCat description'] = funcat_dsc_list

                gene_name_list = parse_list(gene_name)
                entry['subunits(Gene name)'] = gene_name_list

                gene_syn_list = parse_list(gene_syn)
                entry['subunits(Gene name syn)'] = gene_syn_list

                protein_name_list = parse_list(
                    correct_protein_name_list(protein_name))
               	entry['subunits(Protein name)'] = protein_name_list

                # check list lengths match
                if len(protein_name_list) != len(su_entrez_list):
                    msg = 'Unequal number of uniprot/entrez subunits at line {}\n  {}\n  {}'.format(
                        i_entry, '; '.join(protein_name_list), '; '.join(su_entrez_list))
                    raise Exception(msg)

                if len(su_uniprot_list) != len(su_entrez_list):
                    msg = 'Unequal number of uniprot/entrezs subunits at line {}\n  {}\n  {}'.format(
                        i_entry, '; '.join(su_uniprot_list), '; '.join(su_entrez_list))
                    raise Exception(msg)

                # Fix the redundancy issue with swissprot_id field
                if swissprot_id:
                    swissprot_id, _, _ = swissprot_id.partition(';')
                    ncbi_name, _, _ = swissprot_id.partition(' (')
                    result = ncbi_taxa.get_name_translator([ncbi_name])
                    ncbi_id = result[ncbi_name][0]
                else:
                    ncbi_id = None
                entry['SWISSPROT organism (NCBI IDs)'] = ncbi_id
                del entry['SWISSPROT organism']

                file_name = 'corum_' + str(entry['ComplexID']) + '.json'
                full_path = os.path.join(
                    self.cache_dirname, 'corum', file_name)

                with open(full_path, 'w') as f:
                    f.write(json.dumps(entry, indent=4))

                collection = self.con_db()
                collection.insert(entry)

        return collection



'''Helper functions
'''

def parse_list(str_lst):
    """ Parse a semicolon-separated list of strings into a list, ignoring semicolons that are inside square brackets

    Args:
        str_lst (:obj:`str`): semicolon-separated encoding of a list

    Returns:
        :obj:`list` of :obj:`str`: list
    """
    if str_lst:
        lst = []
        depth = 0
        phrase = ''
        for char in str_lst:
            if char == ';' and depth == 0:
                lst.append(phrase)
                phrase = ''
            else:
                if char == '[':
                    depth += 1
                elif char == ']':
                    depth -= 1
                phrase += char
        lst.append(phrase)
        return lst
    else:
        return [None]


def correct_protein_name_list(lst):
    """ Correct a list of protein names with incorrect separators involving '[Cleaved into: ...]'

    Args:
        lst (:obj:`str`): list of protein names with incorrect separators

    Returns:
        :obj:`str`: corrected list of protein names
    """
    if lst:
        lst = lst.replace('[Cleaved into: Nuclear pore complex protein Nup98;',
                          '[Cleaved into: Nuclear pore complex protein Nup98];')
        lst = lst.replace('[Cleaved into: Lamin-A/C;',
                          '[Cleaved into: Lamin-A/C];')
        lst = lst.replace('[Cleaved into: Lamin-A/C ;',
                          '[Cleaved into: Lamin-A/C ];')
        lst = lst.replace('[Includes: Maltase ;', '[Includes: Maltase ];')
    return lst
