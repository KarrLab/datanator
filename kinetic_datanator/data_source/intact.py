# -*- coding: utf-8 -*-

import pandas as pd
from sqlalchemy import Column, Integer, String
import sqlalchemy.ext.declarative
from kinetic_datanator.core import data_source
from six.moves.urllib.request import urlretrieve
import zipfile
from six import BytesIO
from ftplib import FTP
import os


Base = sqlalchemy.ext.declarative.declarative_base()

class ProteinInteractions(Base):
    """ Represents protein interactions in from the IntAct Database

    Attributes:
        Index (:obj:`int`): Index of the DB
        interactor_a (:obj:`str`): represents participant A
        interactor_b (:obj:`str`): represents participant B
        publications (:obj:`str`): resource
        interaction (:obj:`str`): interaction ID
        feature_a (:obj:`str`): binding site of participant A
        feature_b (:obj:`str`): binding site of participant B
        stoich_a (:obj:`str`): stoichiometry of participant A
        stoich_b (:obj:`str`): stoichiometry of participant B

    """

    __tablename__ = 'Protein_Interactions'

    index = Column(Integer, primary_key = True)
    interactor_a = Column(String(255))
    interactor_b = Column(String(255))
    publications = Column(String(255))
    interaction = Column(String(255))
    feature_a = Column(String(255))
    feature_b = Column(String(255))
    stoich_a = Column(String(255))
    stoich_b = Column(String(255))
    interaction_type = Column(String(255))

class ProteinComplex(Base):
    """ Represents protein complexes from the IntAct Database

    Attributes:
        identifier (:obj:`str`):
        name (:obj:`str`):
        ncbi (:obj:`str`):
        subunits (:obj:`str`):
        evidence (:obj:`str`):
        go_annot (:obj:`str`):
        desc (:obj:`str`):
        source (:obj:`str`):
    """
    __tablename__ = 'Protein_Complex'

    identifier = Column(String(255), primary_key = True)
    name = Column(String(255))
    ncbi = Column(String(255))
    subunits = Column(String(255))
    evidence = Column(String(255))
    go_annot = Column(String(255))
    desc = Column(String(255))
    source = Column(String(255))

class IntAct(data_source.HttpDataSource):
    """ A local SQLite copy of the IntAct Database"""
    base_model = Base

    ENDPOINT_DOMAINS = {'intact' : 'ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip',
                        'complex' : 'ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/'}


    def load_content(self):

        #Downloads Content from FTP Server
        self.add_complex()
        self.add_interactions()

    def add_complex(self):
        if not os.path.exists(self.cache_dirname+'/intact_complex'):
            os.makedirs(self.cache_dirname+'/intact_complex')

        ftp = FTP('ftp.ebi.ac.uk')
        ftp.login()
        ftp.cwd('/pub/databases/intact/complex/current/complextab/')
        filenames = ftp.nlst()
        if not os.path.exists(self.cache_dirname+'/intact_complex/'+filenames[0]):
            for filename in filenames:
                local_filename = self.cache_dirname+'/intact_complex/'+filename
                file = open(local_filename, 'wb')
                ftp.retrbinary('RETR '+filename, file.write)
                file.close()
        ftp.quit()

        columns = ['#Complex ac', 'Recommended name', 'Taxonomy identifier',
        'Identifiers (and stoichiometry) of molecules in complex', 'Experimental evidence' ,
        'Go Annotations', 'Description', 'Source']

        new_columns = ['identifier', 'name', 'ncbi', 'subunits', 'evidence', 'go_annot', 'desc', 'source']

        files = os.listdir(self.cache_dirname+'/intact_complex')
        for tsv in files:
            if 'README' in tsv:
                continue
            else:
                dt = pd.read_csv(self.cache_dirname+'/intact_complex/'+tsv, delimiter = '\t', encoding='utf-8')
                pand = dt.loc[:, columns]
                pand.columns = new_columns
                pand = pand.set_index('identifier')
                pand.to_sql(name = 'Protein_Complex', con=self.engine, if_exists = 'append')
                self.session.commit()


    def add_interactions(self):

        if not self.cache_dirname + '/intact.txt':
            path = urlretrieve(self.ENDPOINT_DOMAINS['intact'])
            zipped = zipfile.ZipFile(BytesIO(path[0]))
            zipped.extractall(self.cache_dirname)

        columns = ['#ID(s) interactor A', 'ID(s) interactor B', 'Publication Identifier(s)', 'Interaction identifier(s)',
                    'Feature(s) interactor A', 'Feature(s) interactor B' , 'Stoichiometry(s) interactor A', 'Stoichiometry(s) interactor B',
                    'Interaction type(s)']

        dt = pd.read_csv(self.cache_dirname + '/intact.txt', delimiter = '\t', encoding = 'utf-8')

        pand = dt.loc[:, columns]
        new_columns = ['interactor_a', 'interactor_b', 'publications', 'interaction', 'feature_a', 'feature_b', 'stoich_a', 'stoich_b', 'interaction_type']
        pand.columns = new_columns

        if not self.max_entries == float('inf'):
            pand = pand[0:self.max_entries]

        pand.to_sql(name = 'Protein_Interactions', con=self.engine, if_exists = 'replace', chunksize = 1000)
        self.session.commit()

    # column_list = ["#ID(s) interactor A") , "ID(s) interactor B"), "Alt. ID(s) interactor A"),\
    #     "Alt. ID(s) interactor B"), "Alias(es) interactor A"), "Alias(es) interactor B"),\
    #     "Interaction detection method(s)"), "Publication 1st author(s)"), "Publication Identifier(s)"),\
    #     "Taxid interactor A"), "Taxid interactor B"), "Interaction type(s)"), "Source database(s)"),\
    #     "Interaction identifier(s)"), "Confidence value(s)"), "Expansion method(s)"), \
    #     "Biological role(s) interactor A"), "Biological role(s) interactor B"), \
    #     "Experimental role(s) interactor A"), "Experimental role(s) interactor B"), \
    #     "Type(s) interactor A"), "Type(s) interactor B"), "Xref(s) interactor A"),\
    #     "Xref(s) interactor B"), "Interaction Xref(s)"), "Annotation(s) interactor A"),\
    #     "Annotation(s) interactor B"), "Interaction annotation(s)"), "Host organism(s)"), \
    #     "Interaction parameter(s)"), "Creation date"), "Update date"), \
    #     "Checksum(s) interactor A"), "Checksum(s) interactor B"), "Interaction Checksum(s)"),\
    #     "Negative"), "Feature(s) interactor A"), "Feature(s) interactor B"),\
    #     "Stoichiometry(s) interactor A"), "Stoichiometry(s) interactor B"), "Identification method participant A"), \
    #     "Identification method participant B")]
