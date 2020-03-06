from datanator_query_python.util import mongo_util
import datanator.config.core
from ftplib import FTP
from pathlib import Path
import tempfile
import shutil


class EC(mongo_util.MongoUtil):

    def __init__(self, server=None, db=None, username=None, password=None, 
                 authSource='admin', readPreference='nearest', collection_str='ec',
                 verbose=True, max_entries=float('inf'), cache_dir=None):
        super().__init__(MongoDB=server, db=db, verbose=verbose, max_entries=max_entries,
                        username=username, password=password, authSource=authSource,
                        readPreference=readPreference)
        self.cache_dir = cache_dir
        self.max_entries = max_entries
        self.verbose = verbose
        self.collection_str = collection_str
        self.client, self.db, self.collection = self.con_db(collection_str)

    def establish_ftp(self):
        """establish ftp connection.
        (ftp://ftp.expasy.org/databases/enzyme/enzyme.dat)
        """
        address = "ftp.expasy.org"
        ftp = FTP(address)
        ftp.login()
        ftp.cwd('databases/enzyme')
        return ftp

    def retrieve_content(self):
        """Retrieve content of "enzyme.dat"
        """
        ftp = self.establish_ftp()
        data_file = self.cache_dir+'/enzyme.dat'
        with open(data_file, 'wb+') as fp:
            ftp.retrbinary('RETR enzyme.dat', fp.write)
        ftp.close()
        return data_file

    def parse_content(self, file_location):
        """Parse enzyme.dat file.

        Args:
            file_location(:obj:`str`): location of enzyme.dat file.
        """
        with open(file_location, 'r') as f:
            count = 0
            while count < self.max_entries:
                lines_of_interest = []
                for line in f:
                    if line.strip() == '//':
                        break
                    lines_of_interest.append(line.rstrip('\n'))
                doc = self.make_doc(lines_of_interest)
                if doc == {}:
                    continue
                if self.verbose and count % 10 == 0:
                    print('Updating EC record {}'.format(doc.get('ec_number')))
                self.collection.update_one({'ec_number': doc.get('ec_number')},
                                           {'$set': doc}, upsert=True)
                count += 1

    def make_doc(self, lines):
        """Turn a block of EC info into a dictionary object
        
        Args:
            lines (:obj:`list` of :obj:`str`): list consists of lines of information on one EC group.

        Return:
            (:obj:`dict`): dictionary object.
        """
        ec_synonyms = []
        ec_activity = []
        full_ca = ''
        document = {}        
        for line in lines:
            _id = line[:2]
            data = line[5:]
            if _id == 'ID':
                document['ec_number'] = data.rstrip('.')
            elif _id == 'DE':
                document['ec_name'] = data.rstrip('.')
            elif _id == 'AN':
                ec_synonyms.append(data.rstrip('.'))
                document['ec_synonyms'] = ec_synonyms
            elif _id == 'CA':
                if data[-1] != '.':
                    full_ca += data
                else:
                    full_ca += data.rstrip('.')
                    ec_activity.append(full_ca)
                    document['catalytic_activity'] = ec_activity
                    full_ca = ''
            elif _id == 'CF':
                document['cofactor'] = data.rstrip('.')
            else:
                continue
        return document


def main():
    cache_dir = tempfile.mkdtemp()
    db = 'datanator'
    username = datanator.config.core.get_config()['datanator']['mongodb']['user']
    password = datanator.config.core.get_config()['datanator']['mongodb']['password']
    MongoDB = datanator.config.core.get_config()['datanator']['mongodb']['server']
    src = EC(server=MongoDB, db=db, username=username, password=password, authSource='admin',
            readPreference='nearest', cache_dir=cache_dir)
    path = src.retrieve_content()
    src.parse_content(path)
    shutil.rmtree(path)


if __name__ == "__main__":
    main()