""" Downloads and parses the UnitProt database for protein-protein interactions

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2018-08-15
:Copyright: 2017-2018, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_source
from sqlalchemy import Column, Integer, String, Float
import io
import math
import os
import pandas
import sqlalchemy.ext.declarative


Base = sqlalchemy.ext.declarative.declarative_base()


class UniprotData(Base):
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

    __tablename__ = 'uniprot'
    index = Column(Integer, primary_key=True)
    uniprot_id = Column(String(255), unique=True)
    entry_name = Column(String(255))
    gene_name = Column(String(255))
    protein_name = Column(String(255))
    canonical_sequence = Column(String(255))
    length = Column(Integer)
    mass = Column(Integer)
    ec_number = Column(String(255))
    entrez_id = Column(Integer)
    status = Column(String(255))


class Uniprot(data_source.HttpDataSource):
    """

    """
    base_model = Base
    ENDPOINT_DOMAINS = {
        'uniprot': 'http://www.uniprot.org/uniprot/?fil=reviewed:yes',
    }

    def load_content(self):
        # download data
        url = self.ENDPOINT_DOMAINS['uniprot']
        url += '&columns=id,entry name,genes(PREFERRED),protein names,sequence,length,mass,ec,database(GeneID),reviewed'
        url += '&format=tab'
        url += '&compress=no'
        if not math.isnan(self.max_entries):
            url += '&limit={}'.format(self.max_entries)
        response = self.requests_session.get(url)
        response.raise_for_status()

        # parse data
        table = pandas.read_csv(io.BytesIO(response.content), delimiter='\t', encoding='utf-8')

        # put data into SQLite database
        table.columns = [
            'uniprot_id', 'entry_name', 'gene_name', 'protein_name', 'canonical_sequence', 'length', 'mass',
            'ec_number', 'entrez_id', 'status',
        ]
        table['entrez_id'] = table['entrez_id'].str.replace(';', '')
        table['mass'] = table['mass'].str.replace(',', '')

        table.to_sql(name='uniprot', con=self.engine, if_exists='append', chunksize=1000)
        self.session.commit()
