import datetime
import dateutil.parser
import io
import json
import jxmlease
import requests.exceptions
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
import warnings
import zipfile
from kinetic_datanator.core import data_source


Base = sqlalchemy.ext.declarative.declarative_base()
# :obj:`Base`: base model for local sqlite database

sample_characteristic = sqlalchemy.Table(
    'sample_characteristic', Base.metadata,
    sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
    sqlalchemy.Column('characteristic_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('characteristic._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Sample:Characteristic one-to-many association table

sample_variable = sqlalchemy.Table(
    'sample_variable', Base.metadata,
    sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
    sqlalchemy.Column('variable_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('variable._id')),
)
# :obj:`sqlalchemy.Table`: Sample:Variable one-to-many association table


class Characteristic(Base):
    """ Represents an expiremental characteristic
    Attributes:
        name (:obj:`str`): name of the characteristic (e.g. organism)
        value (:obj:`str`): value of characteristic (e.g Mus musculus)
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    sample_id = sqlalchemy.ForeignKey('sample._id')
    name = sqlalchemy.Column(sqlalchemy.String())
    value = sqlalchemy.Column(sqlalchemy.String())


    __tablename__ = 'characteristic'

class Variable(Base):
    """ Represents an expiremental variable
    Attributes:
        name (:obj:`str`): name of the variable (e.g. genotype)
        value (:obj:`str`): value of characteristic (e.g control)
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    sample_id = sqlalchemy.ForeignKey('sample._id')
    name = sqlalchemy.Column(sqlalchemy.String())
    value = sqlalchemy.Column(sqlalchemy.String())


    __tablename__ = 'variable'


class Sample(Base):
    """ Represents an observed concentration
    Attributes:
        assay_name (:obj:`str`): name of the particular assay in the experiment
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    assay_name = sqlalchemy.Column(sqlalchemy.String())

    characteristics = sqlalchemy.orm.relationship('Characteristic', secondary=sample_characteristic, backref=sqlalchemy.orm.backref('samples'))
    variable = sqlalchemy.orm.relationship('Variable', secondary=sample_variable, backref=sqlalchemy.orm.backref('samples'))

    __tablename__ = 'samples'



class ArrayExpress(data_source.HttpDataSource):
    """ A local sqlite copy of the ECMDB database
    Attributes:
        DOWNLOAD_INDEX_URL (:obj:`str`): URL to download an index of array express experiments
        DOWNLOAD_COMPOUND_URL (:obj:`str`): URL pattern to download samples of a single experiment
    """

    base_model = Base
    ENDPOINT_DOMAIN = 'https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments'
    DOWNLOAD_SAMPLE_URL = ENDPOINT_DOMAIN +'/{}/samples'




    def load_content(self):
        """ Download the content of SABIO-RK and store it to a local sqlite database. """
        db_session = self.session
        req_session = self.requests_session


        list_experiments = ['E-MTAB-5775']

        xml_parser = jxmlease.Parser()
        for i_entry, entry in enumerate(list_experiments):
            
            sample = Sample()

            print self.DOWNLOAD_SAMPLE_URL.format(entry)
            response = req_session.get(self.DOWNLOAD_SAMPLE_URL.format(entry))
            try:
                response.raise_for_status()
            except requests.exceptions.HTTPError:
                warnings.warn('Unable to download data for compound {}'.format(entry), data_source.DataSourceWarning)
                continue

            entry_details = xml_parser(response.text)#['compound']
            print entry_details['experiment']['sample'][0]
            for thing in entry_details['experiment']['sample'][0]:
            	print thing

            print entry_details['experiment']['sample'][0]['assay']['name']
            sample.assay_name = entry_details['experiment']['sample'][0]['assay']['name']


            #for thing in entry_details['experiment']['sample']:
            #	print thing











if __name__ == '__main__':
	blue = ArrayExpress()
	blue.load_content()





