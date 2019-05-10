import json
import requests
from datanator.util import mongo_util


class KeggOrthology(mongo_util.MongoUtil)

    def __init__(self, cache_dirname, MongoDB, db, replicaSet=None, verbose=False, max_entries=float('inf')):

            self.cache_dirname = cache_dirname
            self.MongoDB = MongoDB
            self.db = db
            self.verbose = verbose
            self.max_entries = max_entries
            self.collection = 'kegg_orthology'
            super(KeggOrthology, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, replicaSet=replicaSet, db=db,
                                             verbose=verbose, max_entries=max_entries)

    def get_json_ends(self, tree):
        ends = []
        nodes = [tree]
        while nodes:
            new_nodes = []
            for node in nodes:
                if 'children' in node:
                    new_nodes = new_nodes + node['children']
                else:
                    ends.append(node['name'])
            nodes = new_nodes
        return(ends)

    def load_content(self):
        data = json.load(open('ko00001.json'))
        ends = get_json_ends(data)
        ends = [line.split(" ")[0] for line in ends]
        ends = ends[1:]
        os.makedirs(os.path.join(
            self.cache_dirname, self.collection), exist_ok=True)
        for ko in ends:
            info = requests.get("http://rest.kegg.jp/get/ko:{}".format(ko))
            response.raise_for_status()

    def parse_ko_txt(self):
        
