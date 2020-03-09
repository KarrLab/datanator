from io import BytesIO
import csv
import ete3
import json
import os
import pymongo
import requests
import zipfile
from datanator.util import mongo_util
import datanator.config.core

class CorumNoSQL(mongo_util.MongoUtil):
    def __init__(self, MongoDB, db, replicaSet=None, 
        verbose=False, max_entries=float('inf'), username = None, password = None,
        authSource = 'admin', cache_dirname = None):
        self.ENDPOINT_DOMAINS = {
            'corum': 'https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip',
            'splice': 'https://mips.helmholtz-muenchen.de/corum/download/spliceComplexes.txt.zip'
        }
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.collection = 'corum'
        super(CorumNoSQL, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, db=db,
                    verbose=verbose, max_entries=max_entries, username = username, password = password,
                    authSource = authSource)

    def load_content(self, endpoint='corum'):
        """ Collect and parse all data from CORUM website into JSON files and add to NoSQL database """
        database_url = self.ENDPOINT_DOMAINS[endpoint]
        _, _, collection = self.con_db(self.collection)
        os.makedirs(os.path.join(
            self.cache_dirname, self.collection), exist_ok=True)

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
        if endpoint == 'corum':
            cwd = os.path.join(self.cache_dirname, 'allComplexes.txt')
        else:
            cwd = os.path.join(self.cache_dirname, 'spliceComplexes.txt')

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
                entry['complex_id'] = complex_id #replace string value with int value
                complex_name = entry['ComplexName']
                cell_line = entry['Cell line']
                pur_method = entry['Protein complex purification method']
                # SETS OF INT IDS SEPARATED BY ; eg. GO:0005634
                go_id = entry['GO ID']
                go_dsc = entry['GO description']
                funcat_id = entry['FunCat ID']
                funcat_dsc = entry['FunCat description']
                pubmed_id = int(entry['PubMed ID'])
                entry['pubmed_id'] = pubmed_id
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
                entry['subunits_isoform_id'] = su_uniprot_list
                parsed_su_uniprot_list = parse_subunits(su_uniprot_list)
                entry['subunits_uniprot_id'] = parsed_su_uniprot_list
                del entry['subunits(UniProt IDs)']
                
                su_entrez_list = parse_list(su_entrez)
                entry['subunits_entrez_id'] = su_entrez_list
                del entry['subunits(Entrez IDs)']
                
                go_id_list = parse_list(go_id)
                entry['go_id'] = go_id_list
                del entry['GO ID']
                
                go_dsc_list = parse_list(go_dsc)
                entry['go_description'] = go_dsc_list
                del entry['GO description']

                funcat_id_list = parse_list(funcat_id)
                entry['funcat_id'] = funcat_id_list
                del entry['FunCat ID']

                funcat_dsc_list = parse_list(funcat_dsc)
                entry['funcat_description'] = funcat_dsc_list
                del entry['FunCat description']


                gene_name_list = parse_list(gene_name)
                entry['subunits_gene_name'] = gene_name_list
                del entry['subunits(Gene name)']

                gene_syn_list = parse_list(gene_syn)
                entry['subunits_gene_name_synonym'] = gene_syn_list
                del entry['subunits(Gene name syn)']

                protein_name_list = parse_list(
                    correct_protein_name_list(protein_name))
               	entry['subunits_protein_name'] = protein_name_list
                del entry['subunits(Protein name)']


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
                entry['SWISSPROT_organism_NCBI_ID'] = ncbi_id
                del entry['SWISSPROT organism']

                file_name = 'corum_' + str(entry['complex_id']) + '.json'
                full_path = os.path.join(
                    self.cache_dirname, self.collection, file_name)

                with open(full_path, 'w') as f:
                    f.write(json.dumps(entry, indent=4))

                
                collection.update_one({'ComplexID': entry['ComplexID']},
                    {'$set': entry},
                    upsert=True
                    )
                    

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

def parse_subunits(subunits):
    """Given enzyme subunits list, separate uniprot_id
    and variant. e.g. ["P78381-2"] -> ["P78381-2", "P78381"]
    
    Args:
        subunits (list): corum enzyme subunit string representation
    """
    if subunits == [None]:
        return [None]
    result = []
    for unit in subunits:
        if '-' in unit:
            root = unit.split('-')[0]
            result.append(root)
        else:
            result.append(unit)
    return result


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

def main():
    db = 'datanator'
    cache_dirname = '../../datanator/data_source/cache'
    username = datanator.config.core.get_config()['datanator']['mongodb']['user']
    password = datanator.config.core.get_config()['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
    port = datanator.config.core.get_config()['datanator']['mongodb']['port']
    replSet = datanator.config.core.get_config()['datanator']['mongodb']['replSet']
    manager = CorumNoSQL(MongoDB, db, replicaSet=replSet, cache_dirname = cache_dirname,
            verbose = True, username = username, password = password)
    manager.load_content()

if __name__ == '__main__':
    main()
