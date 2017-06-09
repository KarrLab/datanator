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
# :obj:`sqlalchemy.Table`: Sample:Characteristic many-to-many association table

sample_variable = sqlalchemy.Table(
    'sample_variable', Base.metadata,
    sqlalchemy.Column('sample__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('sample._id'), index=True),
    sqlalchemy.Column('variable_id', sqlalchemy.Integer, sqlalchemy.ForeignKey('variable._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Sample:Variable many-to-many association table


class Characteristic(Base):
    """ Represents an expiremental characteristic
    Attributes:
        samples
        name (:obj:`str`): name of the characteristic (e.g. organism)
        value (:obj:`str`): value of characteristic (e.g Mus musculus)
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String())
    value = sqlalchemy.Column(sqlalchemy.String())

    sqlalchemy.schema.UniqueConstraint(name, value)

    __tablename__ = 'characteristic'


class Variable(Base):
    """ Represents an expiremental variable
    Attributes:
        samples
        name (:obj:`str`): name of the variable (e.g. genotype)
        value (:obj:`str`): value of variable (e.g control)
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String())
    value = sqlalchemy.Column(sqlalchemy.String())

    sqlalchemy.schema.UniqueConstraint(name, value)

    __tablename__ = 'variable'


class Sample(Base):
    """ Represents an observed concentration
    Attributes:
        experiment
        index
        name (:obj:`str`): name of the source of the sample
        characteristics
        variables
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    experiment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('experiment._id'), index=True)
    experiment = sqlalchemy.orm.relationship('Experiment', backref=sqlalchemy.orm.backref('samples'), foreign_keys=[experiment_id])
    index = sqlalchemy.Column(sqlalchemy.Integer())
    name = sqlalchemy.Column(sqlalchemy.String())
    characteristics = sqlalchemy.orm.relationship('Characteristic',
                                                  secondary=sample_characteristic, backref=sqlalchemy.orm.backref('samples'))
    variables = sqlalchemy.orm.relationship('Variable',
                                            secondary=sample_variable, backref=sqlalchemy.orm.backref('samples'))

    sqlalchemy.schema.UniqueConstraint(experiment_id, index)
    sqlalchemy.schema.UniqueConstraint(experiment_id, name)

    __tablename__ = 'sample'


class Experiment(Base):
    """
    Attributes:
        id
        samples
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    id = sqlalchemy.Column(sqlalchemy.String(), unique=True)

    organism = sqlalchemy.Column(sqlalchemy.String())
    name = sqlalchemy.Column(sqlalchemy.String())
    description = sqlalchemy.Column(sqlalchemy.String())
    experiment_type = sqlalchemy.Column(sqlalchemy.String())
    #ax_accession = sqlalchemy.Column(sqlalchemy.String())
    # samples = sqlalchemy.orm.relationship('Sample',
    #                                             secondary=experiment_sample, backref=sqlalchemy.orm.backref('experiment'))

    __tablename__ = 'experiment'


class ExperimentDesign(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    experiment_id = sqlalchemy.Column(sqlalchemy.Integer(), sqlalchemy.ForeignKey('Experiment._id'))
    experiment = sqlalchemy.orm.relationship('Experiment', backref=sqlalchemy.orm.backref(
        'experiment_designs'), foreign_keys=[experiment_id])
    name = sqlalchemy.Column(sqlalchemy.String())

    __tablename__ = 'experiment_design'


class ArrayExpress(data_source.HttpDataSource):
    """ A local sqlite copy of the ECMDB database
    Attributes:
        DOWNLOAD_INDEX_URL (:obj:`str`): URL to download an index of array express experiments
        DOWNLOAD_COMPOUND_URL (:obj:`str`): URL pattern to download samples of a single experiment
    """

    base_model = Base
    ENDPOINT_DOMAIN = 'https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments'
    DOWNLOAD_SAMPLE_URL = ENDPOINT_DOMAIN + '/{}/samples'
    DOWNLOAD_COMPLETE_SAMPLE_URL = 'https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/samples'

    def load_samples(self, experiments):
        """ Download the content of SABIO-RK and store it to a local sqlite database. """
        req_session = self.requests_session

        xml_parser = jxmlease.Parser()
        for experiment in experiments:

            response = req_session.get(self.DOWNLOAD_SAMPLE_URL.format(experiment.id))
            response.raise_for_status()
            entry_details = xml_parser(response.text)

            # create a sample object for each sample in the experiment
            for num, entry in enumerate(entry_details['experiment']['sample']):
                sample = Sample()
                sample.experiment = experiment
                sample.assay_name = entry_details['experiment']['sample'][num]['source']['name']

                # create a characteristic object for each characteristic and append that to the sample's characteristic field
                characteristics = entry_details['experiment']['sample'][num]['characteristic']
                for entry in characteristics:
                    new_charachteristic = Characteristic()
                    new_charachteristic.name = entry['category']
                    new_charachteristic.value = entry['value']
                    sample.characteristics.append(new_charachteristic)

                # create a variable object for each variable and append that to the sample's variable field
                variables = entry_details['experiment']['sample'][num]['variable']
                if isinstance(variables, list):
                    for entry in variables:
                        new_variable = Variable()
                        new_variable.name = entry['name']
                        new_variable.value = entry['value']
                        sample.variables.append(new_variable)
                else:
                    new_variable = Variable()
                    new_variable.name = variables['name']
                    new_variable.value = entry['value']
                    sample.variables.append(new_variable)

    def load_experiments(self, experiment_ids=None):
        db_session = self.session
        req_session = self.requests_session

        if experiment_ids is None:
            url = self.ENDPOINT_DOMAIN
        else:
            url = self.ENDPOINT_DOMAIN + '/' + ','.join([str(id) for id in experiment_ids])
        response = req_session.get(url)
        response.raise_for_status()
        xml_parser = jxmlease.Parser()
        entry_details = xml_parser(response.text)

        for single_entry in entry_details['experiments']['experiment']:
            experiment = Experiment()
            experiment.id = single_entry['accession']
            experiment.name = single_entry['name']
            experiment.experiment_type = single_entry['experimenttype']
            if 'experimentdesign' in single_entry:
                if isinstance(single_entry['experimentdesign'], list):
                    entries = single_entry['experimentdesign']
                else:
                    entries = [single_entry['experimentdesign']]

                for entry in single_entry['experimentdesign']:
                    experiment.experiment_designs.append(ExperimentDesign(name=entry))

            db_session.add(experiment)

    def load_content(self):
        db_session = self.session

        # retrieve list of all experiments
        self.load_experiments()

        # retrieve the samples for all of the experiments
        experiments = db_session.query(Experiment.id).all()
        self.load_samples(experiments)

        # save the changes to the database file
        session.commit()
