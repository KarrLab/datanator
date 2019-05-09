import json
import requests


def get_json_ends(tree):
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


data = json.load(open('ko00001.json'))
ends = get_json_ends(data)
ends = [line.split(" ")[0] for line in ends]
ends = ends[1:]
for ko in ends:
    info = requests.get("http://rest.kegg.jp/get/ko:{}".format(ko))
