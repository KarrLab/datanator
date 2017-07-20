"""
This code takes jaspar database text files and implements them into a SQLlite database

core.py(main) - loops through and parses Jaspar txt files

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017 July 12th
:Copyright: 2017, Karr Lab
:License: MIT
"""

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from parse_jaspar_real import parse_Jaspar_db

database_url = 'http://jaspar.genereg.net/html/DOWNLOAD/database/'

def main(db_name):

    #Create Engine and Session
    print 'Creating Session...'
    engine = create_engine('sqlite:///jaspar.db')
    Session = sessionmaker(bind=engine)
    session = Session()

    #Parse and Add to DB
    parse_Jaspar_db(session, database_url)

    #Commiting and Closing Session
    print 'Finished parsing files, committing to DB.'
    session.commit()
    session.close()
    print 'Successful Creation... Done...'


if __name__ == '__main__':
    main('test.db')
