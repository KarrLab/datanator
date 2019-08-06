from pymongo import collection
import datanator.collections.core


def connect():
    client, database, collection = datanator.collections.core.get_mongo_connection("Metabolites")

    collection.insert_one({"hello":"world"})

if __name__ == "__main__":
    connect()