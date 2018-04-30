""" Downloads and parses the Intact
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-08-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

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


class ProteinInteraction(Base):
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

    __tablename__ = 'Protein_Interaction'

    index = Column(Integer, primary_key = True)
    protein_a = Column(String(255))
    protein_b = Column(String(255))
    gene_a = Column(String(255))
    gene_b = Column(String(255))
    type_a = Column(String(255))
    type_b = Column(String(255))
    role_a = Column(String(255))
    role_b = Column(String(255))
    feature_a = Column(String(255))
    feature_b = Column(String(255))
    stoich_a = Column(String(255))
    stoich_b = Column(String(255))
    method = Column(String(255))
    interaction_id = Column(String(255))
    interaction_type = Column(String(255))
    publication = Column(String(255))
    publication_author = Column(String(255))
    confidence = Column(String(255))

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
                        'intact-partial': 'ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact_negative.txt' ,
                        'complex' : 'ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/'}

    test = True

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
            dt = pd.read_csv(self.cache_dirname+'/intact_complex/'+tsv, delimiter = '\t', encoding='utf-8')
            pand = dt.loc[:, columns]
            pand.columns = new_columns
            pand = pand.set_index('identifier')
            pand.to_sql(name = 'Protein_Complex', con=self.engine, if_exists = 'append')
            self.session.commit()

    def parse_protein(self, interactor, alias):
        protein = gene = None
        if 'uniprotkb' in interactor:
            protein = self.split_colon(interactor)[1]
        else:
            if 'display_short' in alias:
                protein = self.find_between(alias, 'psi-mi:', '(display_short)')
            else:
                protein = None

        for item in self.split_line(alias):
            if '(gene name)' in item:
                gene = self.find_between(item, 'uniprotkb:', '(gene name)')

        return protein, gene

    def split_colon(self, str):
        return str.split(':')

    def split_line(self,str):
        return str.split('|')

    def find_between(self, s, first, last):
        return s[s.index(first) + len(first):s.index(last, s.index(first) + len(first))]

    def find_between_parentheses(self, string):
        if 'psi-mi:' in string:
            return self.find_between(string, '(', ')')
        else:
            return None

    def parse_pubmed(self, string):
        for item in self.split_line(string):
            if 'pubmed:' in item:
                return self.split_colon(item)[1]
        return None

    def add_interactions(self):


        if self.test:
            dt = pd.read_csv(self.ENDPOINT_DOMAINS['intact-partial'], delimiter = '\t', encoding = 'utf-8')
        else:
            if os.path.exists(self.cache_dirname+'/intact.txt') and os.path.exists(self.cache_dirname+'/interact.pkl'):
                pass
            else:
                path = urlretrieve(self.ENDPOINT_DOMAINS['intact'])
                zipped = zipfile.ZipFile(file = path[0])
                zipped.extractall(self.cache_dirname)
                dt = pd.read_csv(self.cache_dirname + '/intact.txt', delimiter = '\t', encoding = 'utf-8')
                dt.to_pickle(self.cache_dirname+'/interact.pkl')

            dt = pd.read_pickle(self.cache_dirname+'/interact.pkl')

        for index, row in dt.iterrows():

            if index == self.max_entries:
                break

            interaction = ProteinInteraction()
            interaction.protein_a, interaction.gene_a = self.parse_protein(row['#ID(s) interactor A'], row['Alias(es) interactor A'])
            interaction.protein_b, interaction.gene_b = self.parse_protein(row['ID(s) interactor B'], row['Alias(es) interactor B'])
            interaction.interaction_type = self.find_between_parentheses(row['Interaction type(s)'])
            interaction.method = self.find_between_parentheses(row['Interaction detection method(s)'])
            interaction.type_a = self.find_between_parentheses(row['Type(s) interactor A'])
            interaction.type_b = self.find_between_parentheses(row['Type(s) interactor B'])
            interaction.role_a = self.find_between_parentheses(row['Biological role(s) interactor A'])
            interaction.role_b = self.find_between_parentheses(row['Biological role(s) interactor B'])
            interaction.feature_a = row['Feature(s) interactor A']
            interaction.feature_b = row['Feature(s) interactor B']
            interaction.stoich_a = row['Stoichiometry(s) interactor A']
            interaction.stoich_b = row['Stoichiometry(s) interactor B']
            interaction.interaction_id = row['Interaction identifier(s)']
            interaction.publication = self.parse_pubmed(row['Publication Identifier(s)'])
            interaction.publication_author = row['Publication 1st author(s)']
            interaction.confidence = row['Confidence value(s)']
            self.session.add(interaction)

        self.session.commit()
