# from uniprot_base import load_uniprot

__metadata__ = {
    '__collection__': 'uniprot',
    'structure': {'uniprot': None},
}


# def load_genedoc(self=None):
#     genedoc_d = load_uniprot()
#     return genedoc_d


def get_mapping(self=None):
    mapping = {
        "uniprot": {
            "dynamic": False,
            #"path": "just_name",
            "properties": {
                "Swiss-Prot": {
                    "type": "string",
                    "analyzer": "string_lowercase"
                },
                "TrEMBL": {
                    "type": "string",
                    "analyzer": "string_lowercase"
                }
            }
        }
    }
    return mapping