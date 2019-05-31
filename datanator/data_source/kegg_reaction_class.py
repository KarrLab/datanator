import json
import requests
import os
from datanator.util import mongo_util


class KeggReaction(mongo_util.MongoUtil):

    def __init__(self, cache_dirname, MongoDB, db, replicaSet=None, verbose=False, max_entries=float('inf'),
                username = None, password = None):
        self.ENDPOINT_DOMAINS = {
            'root': 'https://www.genome.jp/kegg-bin/download_htext?htext=br08204.keg&format=json&filedir=',
        }
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.collection = 'kegg_reaction_class'
        self.path = os.path.join(self.cache_dirname, self.collection)
        super(KeggReaction, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, db=db,
                                            verbose=verbose, max_entries=max_entries,
                                            username = username, password = password)

    def parse_root_json(self):
        '''Parse root json file and return 
            reaction classes
        '''
        root_url = self.ENDPOINT_DOMAINS['root']
        if self.verbose:
            print('\n Downloading root kegg reactions file ...')
        manager = requests.get(root_url)
        manager.raise_for_status()
        os.makedirs(self.path, exist_ok=True)
        file_name = manager.json()['name']
        store_path = os.path.join(self.path, file_name)
        data = manager.json()
        with open(store_path, 'w') as f:
            json.dump(data, f, indent=4)

        names = self.extract_values(data, 'name')
        names = [name.split()[0] for name in names if name[:2] == 'RC']

        return names

    def load_content(self):
        '''Load kegg_reactions into MongoDB
        '''
        _, _, collection = self.con_db(self.collection)

        names = self.parse_root_json()

        iterations = min(len(names), self.max_entries)

        file_format = '.txt'

        i = 0
        for name in names:
            if i == self.max_entries:
                break
            if self.verbose and i % 100 == 0:
                print('Downloading {} of {} kegg reaction class file {}...'.format(
                    i, iterations, name))
            self.download_rxn_cls(name+file_format)
            doc = self.parse_rxn_cls_txt(name+file_format)
            if doc is not None:
                collection.replace_one(
                    {'rclass_id': doc['rclass_id']}, doc, upsert=True)

            i += 1

    def parse_rc_multiline(self, lines):
        ''' Input:
                   DEFINITION  C1y-C2y:*-*:C1b+C8y+N1y-C1b+C8y+N2y
                               N1y-N2y:*-*:C1a+C1x+C1y-C1a+C1x+C2y
                               ...
                               ...
                               O1a-O2x:*-C1z:C1b-C1x
                   RPAIR  CXXXXX 
                   ....
            Output:
                [C1y-C2y:*-*:C1b+C8y+N1y-C1b+C8y+N2y, N1y-N2y:*-*:C1a+C1x+C1y-C1a+C1x+C2y, ...]
        '''
        definition = lines[0].split()[1:]
        if len(lines) > 1:
            for line in lines[1:]:
                definition += line.split()
            return definition
        else: 
            return definition

    def parse_rc_orthology(self, lines):
        '''Input:
                ORTHOLOGY   K00260  glutamate dehydrogenase [EC:1.4.1.2]
                K00261  glutamate dehydrogenase (NAD(P)+) [EC:1.4.1.3]
                K00262  glutamate dehydrogenase (NADP+) [EC:1.4.1.4]
                K00263  leucine dehydrogenase [EC:1.4.1.9]
                ...
                K13547  L-glutamine:2-deoxy-scyllo-inosose/3-amino-2,3-dideoxy-scyllo-inosose aminotransferase [EC:2.6.1.100 2.6.1.101]
                ..
           Output
              [K00260, K00261, ...]
        '''
        ko_id = [lines[0].split()[1]]
        names_str = lines[0].split('  ')[-1].split('[')
        names = [name.strip() for name in names_str[0].split('/')]
        if len(lines) > 1:
            for line in lines[1:]:
                ko_id.append(line.split()[0])
                name_str = line.split('  ')[-1].split('[')
                names.append([name.strip() for name in name_str[0].split(' / ')])
            return (ko_id, names)
        else:
            return (ko_id, names)


    def parse_rxn_cls_txt(self, filename):
        '''Parse kegg_ortho txt file into dictionary object
            categories = ['ENTRY', 'DEFINITION', 'RPAIR', 'REACTION',
                'ENZYME', 'PATHWAY', 'ORTHOLOGY']
        '''
        file_path = os.path.join(self.path, filename)
        try: 
            with open(file_path, 'r') as f:
                doc = {}
                lines = f.readlines()
                # get first word of all the lines
                first_word = [line.split()[0] for line in lines]
                index_definition = first_word.index('DEFINITION') 
                index_rpair = first_word.index('RPAIR')
                index_reaction = first_word.index('REACTION')    # get line number of reaction
                index_enzyme = first_word.index('ENZYME')
                index_orthology = first_word.index('ORTHOLOGY')
                index_pathway = first_word.index('PATHWAY')

                doc['rclass_id'] = lines[0].split()[1]
                doc['definition'] = self.parse_rc_multiline(lines[index_definition:index_rpair])
                doc['reaction_id'] = self.parse_rc_multiline(lines[index_reaction:index_enzyme])
                doc['enzyme'] = self.parse_rc_multiline(lines[index_enzyme:index_pathway])
                ko_id, names = self.parse_rc_orthology(lines[index_orthology:-1])
                doc['orthology_id'] = []
                for _id, name in zip(ko_id, names):
                    doc['orthology_id'].append({'ko_id': _id, 'enzyme_name': name})

                return doc

        except FileNotFoundError as e:
            log_file = os.path.join(self.path, 'kegg_orthology_log.txt')
            with open(log_file, 'a') as f:
                f.write(str(e)+'\n')
            pass

            

    def download_rxn(self, name):
        address = name.split('.')[0]
        try:
            info = requests.get("http://rest.kegg.jp/get/reaction:{}".format(address))
            info.raise_for_status()
            file_name = os.path.join(self.path, name)
            with open(file_name, 'w') as f:
                f.write(info.text)
        except requests.exceptions.HTTPError as e:
            log_file = os.path.join(self.path, 'kegg_rxn_log.txt')
            with open(log_file, 'a') as f:
                f.write(str(e) + '\n')
            pass

    def download_rxn_cls(self, cls):
        address = cls.split('.')[0]
        try:
            info = requests.get("http://rest.kegg.jp/get/rclass:{}".format(address))
            info.raise_for_status()
            file_name = os.path.join(self.path, cls)
            with open(file_name, 'w') as f:
                f.write(info.text)
        except requests.exceptions.HTTPError as e:
            log_file = os.path.join(self.path, 'kegg_rxn_cls_log.txt')
            with open(log_file, 'a') as f:
                f.write(str(e) + '\n')
            pass
