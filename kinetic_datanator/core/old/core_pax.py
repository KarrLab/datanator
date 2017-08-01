"""
This codebase takes the txt files of the PaxDB protein abundance database
and inserts them into an SQL database

core.py(main) - loops through and parses PaxDB txt files

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Date: 2017 June 3
:Copyright: 2017, Karr Lab
:License: MIT
"""

import os
from sqlalchemy        import create_engine
from sqlalchemy.orm    import sessionmaker
from define_pax_tables import Taxon, Dataset, Protein, Observation
from parse_paxDB_files import parse_paxDB_files

data_folder = '/home/balazs/Desktop/aggregator/data/paxdb-abundance-files-v4.0/'

def main(db_name,fraction):

    # Create engine and session
    engine = create_engine('sqlite:///pax.db')
    Session = sessionmaker(bind=engine)
    session = Session()

    # Create report file_name
    report = open('report.txt', 'w+')
    report.write('Errors found:\n')

    # Find data and arse individual files
    data_files = find_files(data_folder)
    n_files = round(fraction*len(data_files),0)

    for file_id in range(0,int(n_files)):
        print('Processing file_id = '+str(file_id+1)+' (out of '+str(int(n_files))+'; '+str(round(100*file_id/n_files,2))+'%'+' already done)')
        parse_paxDB_files(session,report,file_id,data_files,data_folder)

    print('Finished parsing files, committing to DB.')
    session.commit()
    session.close()

""" ------------------------------------------------------------------------ """
def find_files(path):
    """ Scan a directory (and its subdirectories) for files
    Attributes:
        path (:obj:`str`): folder to be scanned
    """

    data_files = []
    for path, subdirs, files in os.walk(path):
        for filename in files:
            f = os.path.join(path, filename)
            data_files.append(f)
    return data_files

if __name__ == '__main__':
    main('test.db',1)
