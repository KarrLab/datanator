from datanator_query_python.config import config
from datanator_query_python.util import mongo_util


class Demo(mongo_util.MongoUtil):

    def __init__(self, 
                 server_demo="someaddress",
                 db_demo="datanator-demo",
                 username_demo="username",
                 password_demo="password",
                 collection_str="demo-collection"):
        super().__init__(MongoDB=server_demo,
                         db=db_demo,
                         username=username_demo,
                         password=password_demo)
        self.collection = self.db_obj[collection_str]

    def update_collection(self):
        """Update collection in db.
        """
        # dic = {"uniprot_id": "P01234",
        #        "locale": "cell membrane",
        #        "array_obj": ["a", "b", "c"]}
        dic = {"uniprot_id": "P01234",
            "locale": "cell membrane",
            "array_obj": ["a", "c", "d"]}

        self.collection.update_one({"uniprot_id": dic["uniprot_id"]},
                                    {"$set": {"locale": dic["locale"]},
                                    "$addToSet": {"array_obj": {"$each": dic["array_obj"]}}},
                                    upsert=True)


def main():
    conf = config.SchemaMigration()
    username = conf.USERNAME
    password = conf.PASSWORD
    server = conf.SERVER
    src = Demo(server_demo=server,
               username_demo=username,
               password_demo=password,
               db_demo="test",
               collection_str="taxon-schema")
    src.update_collection()


if __name__ == "__main__":
    main()