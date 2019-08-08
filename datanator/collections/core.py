import datanator.config.core as config
import datanator.util.mongo_util as mongo_util
import json
import jsonpickle


json.load()


class collection_entry_meta(object):
    """ Contains information regarding the creation of the datbase entry"""

    def __init__(self, Author=None, Date=None, Title=None, Version=None, About=None, type=None):


class collection_entry(object):
	pass


def get_mongo_connection(collection='test'):

    mongo_config = config.get_mongo_config()
    mongo_config["db"] = "datanator"

    return (mongo_util.MongoUtil(**mongo_config).con_db(collection))
