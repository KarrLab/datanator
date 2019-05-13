import json
import requests
import os
from datanator.util import mongo_util

class KeggOrthology(mongo_util.MongoUtil):

    def __init__(self, cache_dirname, MongoDB, db, replicaSet=None, verbose=False, max_entries=float('inf')):
            self.ENDPOINT_DOMAINS = {
                'root': 'https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=json&filedir=',
            }
            self.cache_dirname = cache_dirname
            self.MongoDB = MongoDB
            self.db = db
            self.verbose = verbose
            self.max_entries = max_entries
            self.collection = 'kegg_orthology'
            super(KeggOrthology, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, db=db,
                                             verbose=verbose, max_entries=max_entries)

    def load_content(self):
        root_url = self.ENDPOINT_DOMAINS['root']
        if self.verbose:
            print('\n Downloading root kegg orthology file ...')
        manager = requests.get(root_url)
        manager.raise_for_status()
        path = os.path.join(self.cache_dirname, self.collection)
        os.makedirs(path, exist_ok=True)
        file_name = manager.json()['name']
        store_path = os.path.join(path, file_name)
        data = manager.json()
        with open(store_path, 'w') as f:
            json.dump(data, f, indent=4)

        names = self.extract_values(data, 'name')
        names = [name for name in names if name[0]=='K']
        self.download_ko(names)

    def parse_ko_txt(self, path):
        for filename in os.listdir(path):
            if filename.endswith('.txt'):
                collection = []
                with open(filename, 'r') as f:
                    doc = {}
                    for i, line in enumerate(f):
                        if i == 0: # ENTRY
                            doc['entry'] = line.split()[1]
                            doc['type'] = line.split()[2]
                        elif i ==1: # NAME
                            doc['gene_name'] = [item.replace(',','') for item in line.split()[1:]]
                        elif i == 2: # DEFINITION
                            doc['definition'] = line.split(' ',1)[1]
                        elif line[:5] == 'GENES':
                            continue

        

    def download_ko(self, names):
        path = os.path.join(self.cache_dirname, self.collection)
        os.makedirs(path, exist_ok=True)
        for i, name in enumerate(names):
            if i > self.max_entries:
                break
            if self.verbose and i%10 ==0:
                print('Downloading {} of {} kegg orthology file ...'.format(i, min(len(names), self.max_tries)))
            info = requests.get("http://rest.kegg.jp/get/ko:{}".format(name))
            info.raise_for_status()
            file_name = os.path.join(path, name)
            with open(file_name, 'w') as f:
                f.write(info.text)


    def extract_values(self, obj, key):
        """Pull all values of specified key from nested JSON."""
        arr = []

        def extract(obj, arr, key):
            """Recursively search for values of key in JSON tree."""
            if isinstance(obj, dict):
                for k, v in obj.items():
                    if isinstance(v, (dict, list)):
                        extract(v, arr, key)
                    elif k == key:
                        arr.append(v)
            elif isinstance(obj, list):
                for item in obj:
                    extract(item, arr, key)
            return arr

        results = extract(obj, arr, key)

        return results