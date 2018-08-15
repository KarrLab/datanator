""" Downloads and parses the IntAct database of protein-protein interactions

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-08-13
:Copyright: 2017, Karr Lab
:License: MIT
"""

from ftplib import FTP
from kinetic_datanator.core import data_source
from six import BytesIO
from sqlalchemy import Column, Integer, String
import glob
import pandas
import os
import sqlalchemy.ext.declarative
import zipfile


Base = sqlalchemy.ext.declarative.declarative_base()


class ProteinInteraction(Base):
    """ Represents protein interactions in from the IntAct database

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

    index = Column(Integer, primary_key=True)
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
    """ Represents protein complexes from the IntAct database

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

    identifier = Column(String(255), primary_key=True)
    name = Column(String(255))
    ncbi = Column(String(255))
    subunits = Column(String(255))
    evidence = Column(String(255))
    go_annot = Column(String(255))
    desc = Column(String(255))
    source = Column(String(255))


class IntAct(data_source.FtpDataSource):
    """ A local SQLite copy of the IntAct database """
    base_model = Base

    ENDPOINT_DOMAINS = {
        'psimitab': 'ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact_negative.txt',
        'complextab': 'ftp://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/',
    }

    def get_paths_to_backup(self, download=False):
        """ Get a list of the files to backup/unpack

        Args:
            download (:obj:`bool`, optional): if :obj:`True`, prepare the files for uploading

        Returns:
            :obj:`list` of :obj:`str`: list of paths to backup
        """
        paths = super(data_source.FtpDataSource, self).get_paths_to_backup(download=download)
        paths.append('intact/')
        return paths

    def load_content(self):
        """ Load the content of the local copy of the data source """

        # Download data from FTP Server
        self.download_content()

        # parse data and build SQLite database
        self.add_complexes()
        self.add_interactions()

        # commit changes to database
        self.session.commit()

    def download_content(self):
        """ Download data from FTP server """
        if not os.path.exists(os.path.join(self.cache_dirname, 'intact')):
            os.makedirs(os.path.join(self.cache_dirname, 'intact'))
        if not os.path.exists(os.path.join(self.cache_dirname, 'intact', 'complextab')):
            os.makedirs(os.path.join(self.cache_dirname, 'intact', 'complextab'))
        if not os.path.exists(os.path.join(self.cache_dirname, 'intact', 'psimitab')):
            os.makedirs(os.path.join(self.cache_dirname, 'intact', 'psimitab'))

        ftp = FTP('ftp.ebi.ac.uk')
        ftp.login()

        ftp.cwd('/pub/databases/intact/complex/current/complextab/')
        rel_filenames = ftp.nlst()
        for rel_filename in rel_filenames:
            local_filename = os.path.join(self.cache_dirname, 'intact', 'complextab', rel_filename)
            if not os.path.exists(local_filename):
                with open(local_filename, 'wb') as file:
                    ftp.retrbinary('RETR ' + rel_filename, file.write)

        ftp.cwd('/pub/databases/intact/current/psimitab/')
        local_filename = os.path.join(self.cache_dirname, 'intact', 'psimitab', 'intact_negative.txt')
        with open(local_filename, 'wb') as file:
            ftp.retrbinary('RETR ' + 'intact_negative.txt', file.write)

        ftp.quit()

    def add_complexes(self):
        """ Parse complexes from data and add complexes to SQLite database """
        raw_columns = [
            '#Complex ac', 'Recommended name', 'Taxonomy identifier',
            'Identifiers (and stoichiometry) of molecules in complex', 'Experimental evidence',
            'Go Annotations', 'Description', 'Source',
        ]
        relabeled_columns = ['identifier', 'name', 'ncbi', 'subunits', 'evidence', 'go_annot', 'desc', 'source']

        filenames = glob.glob(os.path.join(self.cache_dirname, 'intact', 'complextab', '*.tsv'))
        for filename in filenames:
            print(filename)
            raw_data = pandas.read_csv(filename, delimiter='\t', encoding='utf-8')
            relabeled_data = raw_data.loc[:, raw_columns]
            relabeled_data.columns = relabeled_columns
            relabeled_data = relabeled_data.set_index('identifier')
            relabeled_data.to_sql(name='Protein_Complex', con=self.engine, if_exists='append')

    def add_interactions(self):
        """ Parse interactions from data and add interactions to SQLite database """
        data = pandas.read_csv(os.path.join(self.cache_dirname, 'intact', 'psimitab', 'intact_negative.txt'),
                               delimiter='\t', encoding='utf-8')
        for index, row in data.iterrows():
            if index == self.max_entries:
                break

            interaction = ProteinInteraction()
            interaction.protein_a, interaction.gene_a = self.find_protein_gene(row['#ID(s) interactor A'], row['Alias(es) interactor A'])
            interaction.protein_b, interaction.gene_b = self.find_protein_gene(row['ID(s) interactor B'], row['Alias(es) interactor B'])
            interaction.interaction_type = self.find_between_psi_mi_parentheses(row['Interaction type(s)'])
            interaction.method = self.find_between_psi_mi_parentheses(row['Interaction detection method(s)'])
            interaction.type_a = self.find_between_psi_mi_parentheses(row['Type(s) interactor A'])
            interaction.type_b = self.find_between_psi_mi_parentheses(row['Type(s) interactor B'])
            interaction.role_a = self.find_between_psi_mi_parentheses(row['Biological role(s) interactor A'])
            interaction.role_b = self.find_between_psi_mi_parentheses(row['Biological role(s) interactor B'])
            interaction.feature_a = row['Feature(s) interactor A']
            interaction.feature_b = row['Feature(s) interactor B']
            interaction.stoich_a = row['Stoichiometry(s) interactor A']
            interaction.stoich_b = row['Stoichiometry(s) interactor B']
            interaction.interaction_id = row['Interaction identifier(s)']
            interaction.publication = self.find_pubmed_id(row['Publication Identifier(s)'])
            interaction.publication_author = row['Publication 1st author(s)']
            interaction.confidence = row['Confidence value(s)']
            self.session.add(interaction)

    def find_protein_gene(self, interactor, alias):
        """ Parse the protein and gene identifiers from key-value pairs of interactors and their aliases

        Args:
            interactor (:obj:`str`): key-value pairs of interactor
            alias (:obj:`str`): key-value pairs of the alias of the interactor

        Returns:
            :obj:`str`: protein identifier
            :obj:`str`: gene identifier
        """
        protein = None
        gene = None
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

    def find_pubmed_id(self, string):
        """ Parse PubMed identifier from annotated key-value pair of publication type-identifier

        Args:
            string (:obj:`str`): key-value pair of publication type-identifier

        Returns:
            :obj:`str`: PubMed identifier
        """
        for item in self.split_line(string):
            if 'pubmed:' in item:
                return self.split_colon(item)[1]
        return None

    def find_between_psi_mi_parentheses(self, string):
        """ Find the text between parentheses in values of psi-mi key-value pairs

        Args:
            string (:obj:`str`): string

        Returns:
            :obj:`str`: substring between the first occurrence of the substring :obj:`first` and the 
                last occurrence of the substring :obj:`last
        """
        if 'psi-mi:' in string:
            return self.find_between(string, '(', ')')
        else:
            return None

    def find_between(self, string, first, last):
        """ Get the substring between the first occurrence of the substring :obj:`first` and the 
        last occurrence of the substring :obj:`last`

        Args:
            string (:obj:`str`): string
            first (:obj:`str`): starting substring
            last (:obj:`str`): ending substring

        Returns:
            :obj:`str`: substring between the first occurrence of the substring :obj:`first` and the 
                last occurrence of the substring :obj:`last
        """
        return string[string.index(first) + len(first):string.index(last, string.index(first) + len(first))]

    def split_colon(self, string):
        """ Split a string into substrings separated by ':'

        Args:
            string (:obj:`str`): string

        Returns:
            :obj:`list`: substring separated by ':'
        """
        return string.split(':')

    def split_line(self, string):
        """ Split a string into substrings separated by '|'

        Args:
            string (:obj:`str`): string

        Returns:
            :obj:`list`: substring separated by '|'
        """
        return string.split('|')
