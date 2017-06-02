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
    characteristics = synonyms = sqlalchemy.orm.relationship('Characteristic', secondary=compound_synonym, backref=sqlalchemy.orm.backref('compounds'))


    __tablename__ = 'characteristics'

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


    __tablename__ = 'variables'


class Sample(Base):
    """ Represents an observed concentration
    Attributes:
        assay_name (:obj:`str`): name of the particular assay in the experiment
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)

    assay_name = sqlalchemy.Column(sqlalchemy.String())

    __tablename__ = 'concentration'



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





