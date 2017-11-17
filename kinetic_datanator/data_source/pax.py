"""
This codebase takes the txt files of the PaxDB protein abundance database
and inserts them into an SQL database

define_tables.py - defines the python classes corresponding to the tables in
the resulting SQL database

:Author: Balazs Szigeti <balazs.szigeti@mssm.edu>
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017 June 3
:Copyright: 2017, Karr Lab
:License: MIT
"""

from sqlalchemy import create_engine, ForeignKey, exists, Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref, sessionmaker
from kinetic_datanator.core import data_source
import os, zipfile
from operator import itemgetter
from six import BytesIO
import shutil
import pandas as pd

Base = declarative_base()

""" --------------------------- Table Definitions ----------------------------------"""
class Taxon(Base):
    """  Represents a species
    Attributes:
        ncbi_id (:obj:`int`): NCBI id
        species_name (:obj:`str`): name and possibly genetic variant
    """

    __tablename__ = 'taxon'
    ncbi_id      = Column(Integer,primary_key=True)
    species_name = Column(String(255))

class Dataset(Base):
    """  Represents a given dataset (typically results form a single paper)
    Attributes:
        ncbi_id     (:obj:`int`): NCBI id - linked to the 'taxon' table
        publication (:obj:`str`): URL of the corresponding publication
        file_name   (:obj:`str`): the name of text file corresponding to the dataset
        score       (:obj:`flt`): PaxDb's internal quality score
        weight      (:obj:`int`): TBA
        coverage    (:obj:`int`): what percentage of the genome is coevred by the datatset
    """

    __tablename__  = 'dataset'
    id          = Column(Integer, primary_key=True)
    publication = Column(String(255))
    file_name   = Column(String(255))
    score       = Column(Float)
    weight      = Column(Integer)
    coverage    = Column(Integer)

    taxon_ncbi_id = Column(Integer, ForeignKey('taxon.ncbi_id'))
    taxon         = relationship('Taxon', backref=backref('datasets'), foreign_keys=[taxon_ncbi_id])

class Protein(Base):
    """  Represents a protein
    Attributes:
        protein_id (:obj:`int`): PaxDB's internal numerical protein ID
        string_id  (:obj:`str`): Ensembl ID of protein
    """

    __tablename__ = 'protein'
    protein_id = Column(Integer, primary_key=True)
    string_id  = Column(String(255))
    uniprot_id = Column(String(255))

class Observation(Base):
    """  Represents a protein
    Attributes:
        protein_id (:obj:`int`): PaxDB's internal numerical protein ID
        dataset_id (:obj:`int`): ID of the database - linked to the 'dataset' table
        abundance  (:obj:`flt`): Normalized abudnance of the protein
    """

    __tablename__ = 'observation'
    id         = Column(Integer, primary_key=True)
    abundance  = Column(String(255))

    dataset_id = Column(Integer, ForeignKey('dataset.id'))
    protein_id = Column(Integer, ForeignKey('protein.protein_id'))
    dataset    = relationship('Dataset',backref=backref('observation'),foreign_keys=[dataset_id])
    protein    = relationship('Protein',backref=backref('observation'),foreign_keys=[protein_id])

    #FORMAT: foreign_table = relationship('foreign_class',backref=backref('self_table'),foreign_keys=[self_column])

""" ------------------------ Method for Finding Files ---------------------------------"""

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



""" ----------------------------- Pax DB Class  ---------------------------------"""


class Pax(data_source.HttpDataSource):
    """ A local sqlite copy of the Pax database

    """

    base_model = Base
    ENDPOINT_DOMAINS = {
        'pax': 'http://pax-db.org/downloads/latest/datasets/paxdb-abundance-files-v4.1.zip',
        'pax_protein': 'http://pax-db.org/downloads/latest/paxdb-uniprot-links-v4.1.zip'
    }

    def load_content(self):
        """ Collects and Parses all data from Pax DB website and adds to SQLlite DB

        Args:
            req (:obj:`requests object`): Requests session object
        """

        database_url = self.ENDPOINT_DOMAINS['pax']
        protein_conversion_url = self.ENDPOINT_DOMAINS['pax_protein']
        req = self.requests_session

        # Extract All Main Files and Save to Cache Directory
        response = req.get(database_url)
        z = zipfile.ZipFile(BytesIO(response.content))
        z.extractall(self.cache_dirname)
        self.cwd = self.cache_dirname+'/paxdb-abundance-files-v4.1'
        shutil.rmtree(self.cwd+'/paxdb-abundance-files-v4.1')
        self.data_files = find_files(self.cwd)

        # Extract All Protein Relation Files and Save to Cache Directory
        response = req.get(protein_conversion_url)
        z = zipfile.ZipFile(BytesIO(response.content))
        z.extractall(self.cache_dirname)
        self.cwd_prot = self.cache_dirname+'/paxdb-uniprot-links-v4.1'


        # Insert Error Report in Cache
        new_file = '/report.txt'
        new_path = self.cache_dirname + '/report.txt'
        self.report = open(new_path, 'w+')
        self.report.write('Errors found:\n')

        self.uniprot_pd = pd.read_csv(self.cwd_prot+'/paxdb-uniprot-links-v4.1.tsv',
            delimiter='\t', names = ['string_id', 'uniprot_id'], index_col=0)

        # Find data and parse individual files
        for self.file_id in range(int(self.max_entries)):
            if self.verbose:
                print('Processing file_id = '+str(self.file_id+1)+' (out of '+str(int(self.max_entries))+'; '+str(round(100*self.file_id/self.max_entries,2))+'%'+' already done)')
            self.parse_paxDB_files()

        # Commit Session
        if self.verbose:
          print('Finished parsing files, committing to DB.')
        self.session.commit()


    def parse_paxDB_files(self):
        """ This function parses pax DB files and adds them to the SQL database
        Attributes:
            session     (:obj:)     : SQLalchemy object
            file_id     (:obj:`str`): internal ID of the file
            data_files  (:obj:`str`): list of the files to be processed
            data_folder (:obj:`str`): root folder of the database
        """
        file_path = self.data_files[self.file_id]
        if self.verbose:
            print(file_path)

        # Get NCBI taxonomy ID from file name
        start  = file_path.find('/',len(self.cwd)-1)+1
        finish = file_path.find('/',len(self.cwd)+4)
        ncbi_id = int(file_path[start:finish])


        # Get file_name
        start   = file_path.find('/',len(self.cwd))+1
        file_name = file_path[start:]

        with open(file_path,'r') as f:
            lines=f.readlines()

            # Get species name
            start  = lines[1].find(':')+2
            finish = lines[1].find('-')-1
            species_name = lines[1][start:finish]

            field_name,_,_ = lines[1].partition(':')
            if field_name=='#name':
                pass
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid #name field, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                return

            # Get score
            finish = len(lines[2])-1
            score = float(lines[2][8:finish])

            field_name,_,_ = lines[2].partition(':')
            if field_name=='#score':
                pass
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid #score field, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                return

            # Get weight
            finish = lines[3].find('%')
            if finish==-1:
                weight = None
            else:
                weight = float(lines[3][9:finish])

            field_name,_,_ = lines[3].partition(':')
            if field_name=='#weight':
                pass
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid #weight field, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                return

            # Get publication link
            start  = lines[4].find('http:')
            finish = lines[4].find('"',start)
            publication = lines[4][start:finish]

            field_name,_,_ = lines[4].partition(':')
            if field_name=='#description':
                pass
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid #description field, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                return

            # Get organ
            start  = lines[5].find(':')+2
            finish = len(lines[5])-1
            organ  = lines[5][start:finish]

            field_name,_,_ = lines[5].partition(':')
            if field_name=='#organ':
                pass
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid #organ field, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                return

            # Get coverage
            start    = lines[7].find(':')+2
            finish   = len(lines[7])-1
            coverage = float(lines[7][start:finish])

            field_name,_,_ = lines[7].partition(':')
            if field_name=='#coverage':
                pass
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid #coverage field, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                return

            # Check column header
            column_headers = lines[11].split()
            if column_headers[0]=='#internal_id' and column_headers[1]=='string_external_id' and column_headers[2]=='abundance' and len(column_headers)<5:
                pass
            else:
                print('Error found, see reports.txt')
                self.report.write('Warning: invalid column headers, excluding file form DB (file_id='+str(self.file_id)+'; '+file_name+')\n')
                return

            """ --- Add taxon and database (metadata info) to session ---------- """

            q = self.session.query(Taxon).filter(Taxon.ncbi_id==ncbi_id)
            if self.session.query(q.exists()).scalar():
                taxon = q.first()
            else:
                taxon = Taxon(ncbi_id=ncbi_id, species_name=species_name)
                self.session.add(taxon)

            dataset = Dataset(publication=publication, file_name=file_name, score=score, weight=weight, coverage=coverage, taxon=taxon)
            self.session.add(dataset)
            #print(taxon.species_name)

            """ --- Parse individual measurements and add them to DB ----------- """

            for i in range(11,len(lines)):
                split_line = lines[i].split()
                protein_id = split_line[0]
                string_id  = split_line[1]
                abundance  = split_line[2]


                # Insert relevant table entries
                q = self.session.query(Protein).filter(Protein.protein_id==protein_id)
                if self.session.query(q.exists()).scalar():
                    protein = q.first()
                else:
                    if string_id in self.uniprot_pd.index.values:
                        protein = Protein(string_id=string_id, uniprot_id = str(self.uniprot_pd.loc[string_id]))
                        self.session.add(protein)
                    else:
                        protein = Protein(string_id=string_id)
                        self.session.add(protein)

                observation = Observation(dataset=dataset, abundance=abundance, protein=protein)
                self.session.add(observation)

        return
