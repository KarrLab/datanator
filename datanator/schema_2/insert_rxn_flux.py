from pathlib import Path
from datanator_query_python.util import mongo_util


class Transform(mongo_util.MongoUtil):
    def __init__(self,
                 MongoDB=None,
                 db=None,
                 des_col=None,
                 username=None,
                 password=None,
                 max_entries=float('inf'),
                 verbose=True):
        super().__init__(MongoDB=MongoDB,
                         db=db,
                         username=username,
                         password=password)
        self.max_entries = max_entries
        self.col = des_col
        self.db = db
        self.verbose = verbose

    def update_docs(self,
                   _dir="../../docs/"):
        """Update database with xls files from
        http://www.cecafdb.org/, stored in ../../docs/

        Args:
            _dir(:obj:`str`): Directory in which xls files are stored.
        """
        paths = Path(_dir).glob('**/*.xls')
        for path in paths:
            print(path)
