import pymongo
from bson.binary import Binary
import pickle
from datanator_query_python.util import mongo_util
import datanator.config.core
from pathlib import Path


def main():
    db = 'test'
    collection_str = 'brenda_constants'
    username = datanator.config.core.get_config()[
        'datanator']['mongodb']['user']
    password = datanator.config.core.get_config(
    )['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config(
    )['datanator']['mongodb']['server']
    manager = mongo_util.MongoUtil(MongoDB=MongoDB, db=db, username=username,
                                   password=password, collection_str=collection_str)

    with open(str(Path('~/karr_lab/datanator/docs/brenda/brenda.pkl').expanduser()), 'rb') as f:
        data = pickle.load(f)
        coll.insert({'bin-data': Binary(thebytes)})