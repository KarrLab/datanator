import csv
import ete3
import os
import warnings
import zipfile
from six import BytesIO


class CorumNoSQL():
    def __init__():
        ENDPOINT_DOMAINS = {
	        'corum': 'https://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip',
	    }

    def load_content(self):
        """ Collect and parse all data from CORUM website and add to SQLite database """
        database_url = self.ENDPOINT_DOMAINS['corum']
        req = self.requests_session
        session = self.session

        # Extract All Files and save to current directory
        response = req.get(database_url)
        z = zipfile.ZipFile(BytesIO(response.content))
        z.extractall(self.cache_dirname)
        self.cwd = os.path.join(self.cache_dirname, 'allComplexes.txt')

        # create object to find NCBI taxonomy IDs
        ncbi_taxa = ete3.NCBITaxa()

        with open(self.cwd, 'r') as file:
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
                complex_name = entry['ComplexName']
                cell_line = entry['Cell line']
                su_uniprot = entry['subunits(UniProt IDs)']  # SETS OF STRING IDS SEPARATED BY ;\
                su_entrez = entry['subunits(Entrez IDs)']  # SETS OF INT IDS SEPARATED BY ;
                pur_method = entry['Protein complex purification method']
                go_id = entry['GO ID']  # SETS OF INT IDS SEPARATED BY ; eg. GO:0005634
                go_dsc = entry['GO description']
                funcat_id = entry['FunCat ID']
                funcat_dsc = entry['FunCat description']
                pubmed_id = int(entry['PubMed ID'])
                protein_name = entry['subunits(Protein name)']
                gene_name = entry['subunits(Gene name)']
                gene_syn = entry['Synonyms']
                disease_cmt = entry['Disease comment']
                su_cmt = entry['Subunits comment']
                complex_cmt = entry['Complex comment']
                swissprot_id = entry['SWISSPROT organism']

                """ ----------------- Apply field level corrections-----------------"""
                # Split the semicolon-separated lists of subunits into protein components,
                # ignoring semicolons inside square brackets
                su_uniprot_list = parse_list(su_uniprot)
                su_entrez_list = parse_list(su_entrez)
                protein_name_list = parse_list(correct_protein_name_list(protein_name))

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
