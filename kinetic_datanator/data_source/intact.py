# -*- coding: utf-8 -*-

import pandas as pd
from sqlalchemy import Column, Integer, String
import sqlalchemy.ext.declarative
from kinetic_datanator.core import data_source
from six.moves.urllib.request import urlretrieve
import zipfile
from six import BytesIO

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


class IntAct(data_source.HttpDataSource):
    """ A local SQLite copy of the IntAct Database"""
    base_model = Base

    ENDPOINT_DOMAINS = {'intact' : 'ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip'}

    def load_content(self):

        #Downloads Content from FTP Server
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
