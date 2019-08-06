import datanator.config.core as config
import datanator.util.mongo_util as mongo_util

def get_mongo_connection (collection='test'):

    mongo_config = config.get_mongo_config()
    mongo_config["db"]= "datanator"


    return (mongo_util.MongoUtil(**mongo_config).con_db(collection))

