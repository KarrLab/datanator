""" Machinery for parsing k_cats and K_ms from `BRENDA <https://www.brenda-enzymes.org/>`_.

#. Download BRENDA text file from https://www.brenda-enzymes.org/download_brenda_without_registration.php.
#. Extract brenda_download.txt from zip archive and save to `./brenda_download.txt`
#. Run the parser in this module::

   import brenda
   data = brenda.Brenda().run()

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2020-04-22
:Copyright: 2020, Karr Lab
:License: MIT
"""

import ete3
import pickle
import json
from pathlib import Path
from datanator_query_python.util import mongo_util
import datanator.config.core
import re
import warnings


class Brenda(mongo_util.MongoUtil):
    RAW_FILENAME = str(Path('~/karr_lab/datanator/docs/brenda/brenda_download.txt').expanduser())
    PROCESSED_FILENAME = str(Path('~/karr_lab/datanator/docs/brenda/brenda.pkl').expanduser())
    MAX_ENTRIES = float('inf')

    LINE_CODES = {
        'AC': 'ACTIVATING_COMPOUND',
        'AP': 'APPLICATION',
        'CF': 'COFACTOR',
        'CL': 'CLONED',
        'CR': 'CRYSTALLIZATION',
        'EN': 'ENGINEERING',
        'EXP': 'EXPRESSION',
        'GI': 'GENERAL_INFORMATION',
        'GS': 'GENERAL_STABILITY',
        'IC50': 'IC50_VALUE',
        'ID': None,
        'IN': 'INHIBITORS',
        'KI': 'KI_VALUE',
        'KKM': 'KCAT_KM_VALUE',
        'KM': 'KM_VALUE',
        'LO': 'LOCALIZATION',
        'ME': 'METALS_IONS',
        'MW': 'MOLECULAR_WEIGHT',
        'NSP': 'NATURAL_SUBSTRATE_PRODUCT',
        'OS': 'OXIDATION_STABILITY',
        'OSS': 'ORGANIC_SOLVENT_STABILITY',
        'PHO': 'PH_OPTIMUM',
        'PHR': 'PH_RANGE',
        'PHS': 'PH_STABILITY',
        'PI': 'PI_VALUE',
        'PM': 'POSTTRANSLATIONAL_MODIFICATION',
        'PR': 'PROTEIN',
        'PU': 'PURIFICATION',
        'RE': 'REACTION',
        'REN': 'RENATURED',
        'RF': 'REFERENCE',
        'RN': 'RECOMMENDED_NAME',
        'RT': 'REACTION_TYPE',
        'SA': 'SPECIFIC_ACTIVITY',
        'SN': 'SYSTEMATIC_NAME',
        'SP': 'SUBSTRATE_PRODUCT',
        'SS': 'STORAGE_STABILITY',
        'ST': 'SOURCE_TISSUE',
        'SU': 'SUBUNITS',
        'SY': 'SYNONYMS',
        'TN': 'TURNOVER_NUMBER',
        'TO': 'TEMPERATURE_OPTIMUM',
        'TR': 'TEMPERATURE_RANGE',
        'TS': 'TEMPERATURE_STABILITY'
    }

    def __init__(self, MongoDB=None, db=None, username=None, password=None,
                collection_str=None, authSource='admin', readPreference='nearest'):
        super().__init__(MongoDB=MongoDB, db=db, username=username, password=password,
                        authSource=authSource, readPreference=readPreference)
        self._ncbi_taxa = ete3.NCBITaxa()
        self.collection = self.db_obj[collection_str]

    def run(self, raw_filename=None, processed_filename=None, max_entries=None):
        raw_filename = raw_filename or self.RAW_FILENAME
        processed_filename = str(Path(processed_filename).expanduser()) or self.PROCESSED_FILENAME
        max_entries = max_entries or self.MAX_ENTRIES

        ec_data = None
        ec_code = None
        section_name = None
        val_type = None
        val_content = None
        data = {}
        with open(raw_filename, 'r') as file:
            for line in file:
                if line.endswith('\n'):
                    line = line[0:-1]

                # skip blank lines
                if not line:
                    continue

                # skip comment lines
                if line.startswith('*'):
                    continue

                # EC section separator
                if line.startswith('///'):
                    self.parse_content(ec_code, ec_data, val_type, val_content)

                    ec_data = None
                    ec_code = None
                    section_name = None
                    val_type = None
                    val_content = None
                    continue

                # process section heading
                if '\t' not in line:
                    if val_type == 'ID':
                        ec_code, ec_data = self.parse_ec_code(data, val_content)
                    else:
                        self.parse_content(ec_code, ec_data, val_type, val_content)

                    section_name = line
                    val_type = None
                    val_content = None

                # process content line
                else:
                    line_type, _, line_content = line.partition('\t')
                    if line_type:
                        assert self.LINE_CODES[line_type] == section_name

                    if val_type is None:
                        val_type = line_type
                        val_content = line_content

                    elif line_type:
                        if val_type == 'ID':
                            ec_data = self.parse_ec_code(data, val_content)
                        else:
                            self.parse_content(ec_code, ec_data, val_type, val_content)

                        val_type = line_type
                        val_content = line_content
                    else:
                        val_content += '\n' + line_content

        # link refs
        for ec_data in data.values():
            # for enz in ec_data['enzymes'].values():
            #     enz['refs'] = [ec_data['refs'][ref_id] for ref_id in enz['ref_ids']]

                # for tissue in enz['tissues']:
                #     tissue['refs'] = [ec_data['refs'][ref_id] for ref_id in tissue['ref_ids']]

                # for loc in enz['subcellular_localizations']:
                #     loc['refs'] = [ec_data['refs'][ref_id] for ref_id in loc['ref_ids']]

            # for rxn in ec_data['natural_reactions']:
            #     rxn['refs'] = [ec_data['refs'][ref_id] for ref_id in rxn['ref_ids']]

            # for rxn in ec_data['reactions']:
            #     rxn['refs'] = [ec_data['refs'][ref_id] for ref_id in rxn['ref_ids']]

            for k_cat in ec_data['k_cats']:
                k_cat['refs'] = [ec_data['refs'][ref_id] for ref_id in k_cat['ref_ids']]

            for k_m in ec_data['k_ms']:
                k_m['refs'] = [ec_data['refs'][ref_id] for ref_id in k_m['ref_ids']]

        # remove information no longer needed because refs have been deserialized
        for i, (_key, ec_data) in enumerate(data.items()):
            if i == max_entries:
                break
            if i % 10 == 0:
                print('Processing EC {}, the {}th of all ECs'.format(_key, i))
            for enz in ec_data['enzymes'].values():
                enz.pop('id')
            #     enz.pop('ref_ids')

            #     for tissue in enz['tissues']:
            #         tissue.pop('ref_ids')

            #     for loc in enz['subcellular_localizations']:
            #         loc.pop('ref_ids')

            # for rxn in ec_data['natural_reactions']:
            #     rxn.pop('ref_ids')

            # for rxn in ec_data['reactions']:
            #     rxn.pop('ref_ids')

            for k_cat in ec_data['k_cats']:
                k_cat.pop('ref_ids')

            for k_m in ec_data['k_ms']:
                k_m.pop('ref_ids')

            for ref in ec_data['refs'].values():
                ref.pop('id')

            ec_data.pop('enzymes')
            ec_data.pop('refs')

            # save to MongoDB
            self.collection.update_one({'ec_number': _key},
                                       {'$addToSet': {'ec_synonyms': {'$each': [ec_data['name'], ec_data['systematic_name']]}},
                                        '$set': {'k_ms': ec_data['k_ms'],
                                                'k_cats': ec_data['k_cats'],
                                                'comments': ec_data['comments']}}, upsert=True)
        # return extracted data
        return 'done!!'

    def parse_ec_code(self, data, val):
        match = re.match(r'^([0-9\.]+)([ \n]\((.*?)\))?$', val, re.DOTALL)

        ec_code = match.group(1)
        print('Parsing EC code {}'.format(ec_code))

        ec_data = data[ec_code] = {
            'name': None,
            'systematic_name': None,
            # 'natural_reactions': [],
            # 'reactions': [],
            'enzymes': {},
            'k_cats': [],
            'k_ms': [],
            'comments': None,
            'refs': {},
        }

        if match.group(3):
            ec_data['comments'] = match.group(3).replace('\n', ' ').strip()

        return ec_code, ec_data

    def parse_content(self, ec_code, ec_data, type, val):
        if type == 'PR':
            match = re.match((
                r'^#(\d+)#'
                r'[ \n](.*?)'
                r'([ \n]([A-Z][A-Z0-9\.]+)'
                r'[ \n](GenBank|UniProt|SwissProt|))?'
                r'([ \n]\((.*?)\))?'
                r'[ \n]<([0-9,\n]+)>$'
            ), val, re.DOTALL)

            id = match.group(1)
            taxon_name = match.group(2).replace('\n', ' ').strip()
            taxon_id = self._ncbi_taxa.get_name_translator([taxon_name]).get(taxon_name, [None])[0]
            xid = match.group(4) or None
            namespace = match.group(5) or None
            comments = self.parse_comments(match.group(7))
            ref_ids = match.group(8).replace('\n', ',').strip().split(',')

            if id not in ec_data['enzymes']:
                ec_data['enzymes'][id] = {
                    'id': id,
                    'identifiers': [{'namespace': namespace, 'id': xid}] if xid else [],
                    'taxon': {'name': taxon_name, 'id': taxon_id} if taxon_name else None,
                    'tissues': [],
                    'subcellular_localizations': [],
                    # 'comments': comments,
                    # 'ref_ids': ref_ids,
                    # 'refs': None,
                }
            else:
                if xid:
                    ec_data['enzymes'][id]['identifiers'].append({'namespace': namespace, 'id': xid})

                    ec_data['enzymes'][id]['identifiers'] = [{'namespace': ns, 'id': xid}
                                                             for ns, xid in set((id['namespace'], id['id'])
                                                                                for id in ec_data['enzymes'][id]['identifiers'])]

                if not ec_data['enzymes'][id]['taxon'] and taxon_name:
                    ec_data['enzymes'][id]['taxon'] = {'name': taxon_name, 'id': taxon_id}

                # ec_data['enzymes'][id]['comments'] += comments
                # ec_data['enzymes'][id]['ref_ids'] = sorted(set(ec_data['enzymes'][id]['ref_ids'] + ref_ids))

        elif type == 'RN':
            ec_data['name'] = val.replace('\n', ' ').strip()

        elif type == 'SN':
            ec_data['systematic_name'] = val.replace('\n', ' ').strip()

        elif type == 'ST':
            match = re.match(r'^#(.*?)#[ \n](.*?)([ \n]\((.*?)\))?[ \n]<([0-9,\n]+)>$', val, re.DOTALL)

            enz_ids = match.group(1).replace('\n', ',').strip().split(',')
            tissue = match.group(2).strip()
            comments = self.parse_comments(match.group(4))
            ref_ids = match.group(5).replace('\n', ',').strip().split(',')
            for enz_id in enz_ids:
                enzyme = ec_data['enzymes'].get(enz_id, None)
                if enzyme:
                    enzyme['tissues'].append({
                        'name': tissue,
                        'comments': self.filter_comments(comments, enz_id),
                        # 'ref_ids': ref_ids,
                        # 'refs': None,
                    })
                else:
                    warnings.warn('{} does not have enzyme with id {}. Error due to {}'.format(ec_code, enz_id, val), UserWarning)

        elif type == 'LO':
            match = re.match(r'^#(.*?)#[ \n](.*?)([ \n]\((.*?)\))?[ \n]<([0-9,\n]+)>$', val, re.DOTALL)

            enz_ids = match.group(1).replace('\n', ',').strip().split(',')
            localization = match.group(2).strip()
            comments = self.parse_comments(match.group(4))
            ref_ids = match.group(5).replace('\n', ',').strip().split(',')
            for enz_id in enz_ids:
                enzyme = ec_data['enzymes'].get(enz_id, None)
                if enzyme:
                    ec_data['enzymes'][enz_id]['subcellular_localizations'].append({
                        'name': localization,
                        'comments': self.filter_comments(comments, enz_id),
                        # 'ref_ids': ref_ids,
                        # 'refs': None,
                    })
                else:
                    warnings.warn('{} does not have enzyme with id {}. Error due to {}'.format(ec_code, enz_id, val), UserWarning)

        # elif type == 'NSP':
        #     match = re.match(r'^#(.*?)#[ \n](.*?)([ \n]\((.*?)\))?([ \n]\|.*?\|)?([ \n]\{(r|)\})?[ \n]<([0-9,\n]+)>$', val, re.DOTALL)
        #     comments = self.parse_comments(match.group(4))
        #     ref_ids = match.group(8).replace('\n', ',').strip().split(',')
        #     ec_data['natural_reactions'].append({
        #         'equation': match.group(2).replace('\n', ' ').strip(),
        #         'reversible': match.group(7) == 'r',
        #         'enz_ids': match.group(1).replace('\n', ',').strip().split(','),
        #         'comments': comments,
        #         # 'ref_ids': ref_ids,
        #         # 'refs': None,
        #     })

        # elif type == 'SP':
        #     match = re.match(r'^#(.*?)#[ \n](.*?)([ \n]\((.*?)\))?([ \n]\|.*?\|)?([ \n]\{(r|)\})?[ \n]<([0-9,\n]+)>$', val, re.DOTALL)
        #     comments = self.parse_comments(match.group(4))
        #     ref_ids = match.group(8).replace('\n', ',').strip().split(',')
        #     ec_data['reactions'].append({
        #         'equation': match.group(2).replace('\n', ' ').strip(),
        #         'reversible': match.group(7) == 'r',
        #         'enz_ids': match.group(1).replace('\n', ',').strip().split(','),
        #         'comments': comments,
        #         # 'ref_ids': ref_ids,
        #         # 'refs': None,
        #     })

        elif type in ['TN', 'KM']:
            match = re.match(
                '^#(.*?)#[ \n]((\d+(\.\d+)?)?-?(\d+(\.\d+)?)?)[ \n]{(.*?)}[ \n]([ \n]\((.*?)\))?[ \n]?<([0-9,\n]+)>$', val, re.DOTALL)

            enz_ids = match.group(1).replace('\n', ',').split(',')
            ref_ids = match.group(10).replace('\n', ',').strip().split(',')
            comments = self.parse_comments(match.group(9))

            datum = {
                'substrate': match.group(7),
                'value': match.group(2).replace('\n', ''),
                'units': 's^-1' if type == 'TN' else 'mM',
                'enzymes': [],
                'ref_ids': ref_ids,
                'refs': None,
            }

            for enz_id in enz_ids:
                enzyme = ec_data['enzymes'].get(enz_id, None)
                if enzyme:
                    filtered_comments = self.filter_comments(comments, enz_id)

                    datum['enzymes'].append({
                        'enzyme': enzyme,
                        'comments': filtered_comments,
                    })

                else:
                    warnings.warn('{} does not have enzyme with id {}. Error due to {}'.format(ec_code, enz_id, val), UserWarning)

            if type == 'TN':
                ec_data['k_cats'].append(datum)
            else:
                ec_data['k_ms'].append(datum)

        elif type == 'RF':
            match = re.match(
                '^<(\d+)>[ \n](.*?)[ \n]{(Pubmed):(.*?)}([ \n]\((.*?)\))?$', val, re.DOTALL)

            citation = match.group(2).replace('\n', ' ')
            citation_match = re.match(r'^(.*?):[ \n](.*?)\.[ \n](.*?)[ \n]\((\d+)\)[ \n](.*?),[ \n](.*?)\.$', citation, re.DOTALL)

            ec_data['refs'][match.group(1)] = {
                'id': match.group(1),
                'authors': [author.replace('\n', ' ') for author in re.split(r';[ \n]', citation_match.group(1))],
                'title': citation_match.group(2).replace('\n', ' '),
                'journal': citation_match.group(3).replace('\n', ' '),
                'volume': citation_match.group(5).replace('\n', ' '),
                'pages': citation_match.group(6).replace('\n', ' '),
                'year': int(citation_match.group(4)),
                'identifier': {
                    'namespace': match.group(3),
                    'id': match.group(4),
                },
                'comments': match.group(6).replace('\n', ' ') if match.group(6) else None,
            }

    def parse_comments(self, comment_strs):
        if not comment_strs:
            return []

        comment_strs = re.split(r';[ \n]', comment_strs.strip())
        comments = []
        for comment_str in comment_strs:
            comment_match = re.match(r'^#([0-9,\n]+)#[ \n](.*?)[ \n]<([0-9,\n]+)>$', comment_str.strip())
            if comment_match:
                enz_ids = comment_match.group(1).replace('\n', ',').strip().split(',')
                ref_ids = comment_match.group(3).replace('\n', ',').strip().split(',')

                text = comment_match.group(2).replace('\n', ' ')
                if 'wild-type' in text or 'native' in text:
                    wild_type = True
                else:
                    wild_type = None

                if 'mutant' in text:
                    genetic_variant = True
                else:
                    genetic_variant = None

                match = re.search(r'(^|,[ \n])(\d+(\.\d+)?)°C(,[ \n]|$)', text)
                if match:
                    temperature = float(match.group(2))
                else:
                    temperature = None

                match = re.search(r'(^|,[ \n])pH[ \n](\d+(\.\d+)?)(,[ \n]|$)', text, re.DOTALL)
                if match:
                    ph = float(match.group(2))
                else:
                    ph = None

                comments.append({
                    'enz_ref_ids': list(zip(enz_ids, ref_ids)),
                    'text': text,
                    'wild_type': wild_type,
                    'genetic_variant': genetic_variant,
                    'temperature': temperature,
                    'ph': ph,
                })
            else:
                comments.append({
                    'enz_ref_ids': [(None, None)],
                    'text': comment_str,
                    'wild_type': None,
                    'genetic_variant': None,
                    'temperature': None,
                    'ph': None,
                })

        return comments

    def filter_comments(self, comments, enz_id):
        filtered_comments = []
        for comment in comments:
            has_comment = False
            for enz_ref_id in comment['enz_ref_ids']:
                if enz_ref_id[0] == enz_id:
                    has_comment = True
                    break
            if has_comment:
                filtered_comments.append(comment)
        return filtered_comments


"""
names = set()
systematic_names = set()
reaction_eqs = set()
k_cat_vals = set()
k_m_vals = set()
comments = set()
for ec_data in data.values():
    names.add(ec_data['name'])
    systematic_names.add(ec_data['systematic_name'])
    reaction_eqs.update(set(rxn['equation'] for rxn in ec_data['reactions']))
    k_cat_vals.update(set(k_cat['value'] for k_cat in ec_data['k_cats']))
    k_m_vals.update(set(k_m['value'] for k_m in ec_data['k_ms']))
    comments.add(ec_data['comments'])
"""

def main():
    db = 'datanator'
    collection_str = 'ec'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    server = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    status = Brenda(MongoDB=server, username=username, password=password,
           db=db, collection_str=collection_str).run(processed_filename='~/karr_lab/datanator/docs/brenda/brenda-1.json')
    print(status)

if __name__ == '__main__':
    main()