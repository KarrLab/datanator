import os
import zipfile
from six import BytesIO
import shutil
import pandas as pd
import requests
from datanator.util import mongo_util


class PaxNoSQL(mongo_util.MongoUtil):

    def __init__(self, cache_dirname, MongoDB, db, verbose=False, 
                max_entries=float('inf'), username = None, password = None,
                authSource = 'admin', replicaSet = None):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.ENDPOINT_DOMAINS = {
            'pax': 'https://pax-db.org/downloads/4.1/datasets/paxdb-abundance-files-v4.1.zip',
            'pax_protein': 'http://pax-db.org/downloads/latest/paxdb-uniprot-links-v4.1.zip'
        }
        self.collection = 'pax'
        super(PaxNoSQL, self).__init__(cache_dirname=cache_dirname, MongoDB=MongoDB, 
                replicaSet=replicaSet, db=db,
                 verbose=verbose, max_entries=max_entries, username = username, 
                 password = password, authSource = authSource)


    def load_content(self):
        """ Collects and Parses all data from Pax DB website and adds to MongoDB

        """
        client, db_obj, collection = self.con_db(self.collection)
        database_url = self.ENDPOINT_DOMAINS['pax']
        protein_conversion_url = self.ENDPOINT_DOMAINS['pax_protein']

        # Extract All Main Files and Save to Cache Directory
        response = requests.get(database_url)
        response.raise_for_status()
        z = zipfile.ZipFile(BytesIO(response.content))
        z.extractall(self.cache_dirname)
        self.cwd = self.cache_dirname+'/paxdb-abundance-files-v4.1'
        shutil.rmtree(self.cwd+'/paxdb-abundance-files-v4.1')
        self.data_files = find_files(self.cwd)

        # Extract All Protein Relation Files and Save to Cache Directory
        response = requests.get(protein_conversion_url)
        response.raise_for_status()
        z = zipfile.ZipFile(BytesIO(response.content))
        z.extractall(self.cache_dirname)
        self.cwd_prot = self.cache_dirname+'/paxdb-uniprot-links-v4.1'

        # Insert Error Report in Cache
        new_file = '/report.txt'
        new_path = self.cache_dirname + '/report.txt'
        self.report = open(new_path, 'w+')
        self.report.write('Errors found:\n')

        self.uniprot_pd = pd.read_csv(self.cwd_prot+'/paxdb-uniprot-links-v4.1.tsv',
                                      delimiter='\t', header=None, names=['string_id', 'uniprot_id'], index_col=0)

        # Find data and parse individual files
        file_count = 0
        for _,_,files in (os.walk(self.cache_dirname)):
        	file_count += len(files)
        total = min(file_count, self.max_entries)

        i = 0
        for self.file_id in range(total):
            if self.verbose:
                print('Processing file_id = '+str(self.file_id+1)+' (out of '+str(total) +
                      '; '+str(round(100*self.file_id/total, 2))+'%'+' already done)')
            if i > self.max_entries:
                break
            entry = self.parse_paxDB_files()
            collection.replace_one({'file_name':entry['file_name']},entry,upsert=True)
            i += 1

        client.close()

        return collection

    def parse_paxDB_files(self):
        """ This function parses pax DB files and adds them to the NoSQL database
        """
        file_path = self.data_files[self.file_id]
        entry = {}
        if self.verbose:
            print(file_path)

        # Get NCBI taxonomy ID from file name
        start = file_path.find('/', len(self.cwd)-1)+1
        finish = file_path.find('/', len(self.cwd)+4)
        ncbi_id = int(file_path[start:finish])

        entry['ncbi_id'] = ncbi_id

        # Get file_name
        start = file_path.find('/', len(self.cwd))+1
        file_name = file_path[start:]

        entry['file_name'] = file_name

        with open(file_path, 'r') as f:
            lines = f.readlines()

            # Get species name
            start = lines[1].find(':')+2
            finish = lines[1].find('-')-1
            species_name = lines[1][start:finish]

            field_name, _, _ = lines[1].partition(':')
            if field_name == '#name':
                entry['species_name'] = species_name
            else:
                print('Error found, see reports.txt')
                self.report.write(
                    'Warning: invalid #name field, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                entry['species_name'] = None
                return

            # Get score
            finish = len(lines[2])-1
            score = float(lines[2][8:finish])

            field_name, _, _ = lines[2].partition(':')
            if field_name == '#score':
                entry['score'] = score
            else:
                print('Error found, see reports.txt')
                self.report.write(
                    'Warning: invalid #score field, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                entry['score'] = None
                return

            # Get weight
            finish = lines[3].find('%')
            if finish == -1:
                weight = None
            else:
                weight = float(lines[3][9:finish])

            field_name, _, _ = lines[3].partition(':')
            if field_name == '#weight':
                entry['weight'] = weight
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid #weight field, excluding file form DB (file_id=' +
                                  str(self.file_id)+'; '+file_name+')\n')
                entry['weight'] = None
                return

            # Get publication link
            start = lines[4].find('http:')
            finish = lines[4].find('"', start)
            publication = lines[4][start:finish]

            field_name, _, _ = lines[4].partition(':')
            if field_name == '#description':
                entry['publication'] = publication
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid #description field, excluding file form DB (file_id=' +
                                  str(self.file_id)+'; '+file_name+')\n')
                entry['publication'] = None
                return

            # Get organ
            start = lines[5].find(':')+2
            finish = len(lines[5])-1
            organ = lines[5][start:finish]

            field_name, _, _ = lines[5].partition(':')
            if field_name == '#organ':
                entry['organ'] = organ
            else:
                print('Error found, see reports.txt')
                self.report.write(
                    'Warning: invalid #organ field, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                entry['organ'] = None
                return

            # Get coverage
            start = lines[7].find(':')+2
            finish = len(lines[7])-1
            coverage = float(lines[7][start:finish])

            field_name, _, _ = lines[7].partition(':')
            if field_name == '#coverage':
                entry['coverage'] = coverage
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid #coverage field, excluding file form DB (file_id=' +
                                  str(self.file_id)+'; '+file_name+')\n')
                entry['coverage'] = None
                return

            # Check column header
            column_headers = lines[11].split()
            if column_headers[0] == '#internal_id' \
                    and column_headers[1] == 'string_external_id' \
                    and column_headers[2] == 'abundance' \
                    and len(column_headers) < 5:
                pass
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid column headers, excluding file form DB (file_id=' +
                                  str(self.file_id)+'; '+file_name+')\n')
                return

            """ --- Parse individual measurements and add them to DB ----------- """
            observation = []
            for i in range(12, len(lines)):
                split_line = lines[i].split()
                protein_id = split_line[0]
                string_id = split_line[1]
                abundance = split_line[2]
                protein_info = None #default
                if string_id in self.uniprot_pd.index.values:
                    protein_info = {'string_id': string_id, 'uniprot_id':str(
                        self.uniprot_pd.loc[string_id]['uniprot_id'])}
                else:
                    pass

                observation.append(
                    {'protein_id': protein_info, 'string_id': string_id, 'abundance': abundance})
            entry['observation'] = observation
        return entry


'''Helper functions
'''


def find_files(path):
    """ Scan a directory (and its subdirectories) for files and sort by ncbi_id

    Args:
        path (:obj:`str`): Path containing the data_files

    Returns:
        :obj:`list`: list of files to add to DB

    """

    data_files = []
    temp = []
    for path, subdirs, files in os.walk(path):
        temp = subdirs
        break
    subdir = sorted([int(x) for x in temp])
    for items in subdir:
        compound = path + '/' + str(items)
        files = [f for f in os.listdir(compound)]
        for filename in files:
            f = os.path.join(compound, filename)
            data_files.append(f)
    return data_files
