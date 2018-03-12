""" Downloads and parses the ArrayExpress database
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-08-16
:Copyright: 2017, Karr Lab
:License: MIT
"""

import datetime
import dateutil.parser
import pkg_resources
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
from kinetic_datanator.core import data_source
from Bio import SeqIO


Base = sqlalchemy.ext.declarative.declarative_base()
# :obj:`Base`: base model for local sqlite database

geneEntry_geneIdentifier = sqlalchemy.Table(
    'geneEntry_geneIdentifier', Base.metadata,
    sqlalchemy.Column('gene__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('gene._id'), index=True),
    sqlalchemy.Column('identifier__id', sqlalchemy.Integer, sqlalchemy.ForeignKey('identifier._id'), index=True),
)
# :obj:`sqlalchemy.Table`: Sample:Characteristic many-to-many association table

class GeneEntry(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String())
    organism = sqlalchemy.Column(sqlalchemy.String())
    identifiers = sqlalchemy.orm.relationship(
        'GeneIdentifier', secondary=geneEntry_geneIdentifier, backref=sqlalchemy.orm.backref('samples'))
    exp_id = sqlalchemy.Column(sqlalchemy.String())
    samp_id = sqlalchemy.Column(sqlalchemy.String())

    __tablename__ = 'gene'

class GeneIdentifier(Base):
    """ Represents a url
    Attributes:
        _id (:obj:`int`): unique id
        category (:obj:`str`): name of the characteristic (e.g. organism)
        value (:obj:`str`): value of characteristic (e.g. Mus musculus)
        samples (:obj:`list` of :obj:`Sample`): samples
    """
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String())
    sqlalchemy.schema.UniqueConstraint(name)

    __tablename__ = 'identifier'



class GetGenes(data_source.HttpDataSource):
    """ A local sqlite copy of the ArrayExpress database
    Attributes:
        EXCLUDED_DATASET_IDS (:obj:`list` of :obj:`str`): list of IDs of datasets to exclude
    """

    base_model = Base


    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backup=True, verbose=False,
                 clear_requests_cache=False, download_request_backup=False):
        """
        Args:
            name (:obj:`str`, optional): name
            cache_dirname (:obj:`str`, optional): directory to store the local copy of the data source and the HTTP requests cache
            clear_content (:obj:`bool`, optional): if :obj:`True`, clear the content of the sqlite local copy of the data source
            load_content (:obj:`bool`, optional): if :obj:`True`, load the content of the local sqlite database from the external source
            max_entries (:obj:`float`, optional): maximum number of entries to save locally
            commit_intermediate_results (:obj:`bool`, optional): if :obj:`True`, commit the changes throughout the loading
                process. This is particularly helpful for restarting this method when webservices go offline.
            download_backup (:obj:`bool`, optional): if :obj:`True`, load the local copy of the data source from the Karr Lab server
            verbose (:obj:`bool`, optional): if :obj:`True`, print status information to the standard output
            clear_requests_cache (:obj:`bool`, optional): if :obj:`True`, clear the HTTP requests cache
            download_request_backup (:obj:`bool`, optional): if :obj:`True`, download the request backup
        """
        super(GetGenes, self).__init__(name=name, cache_dirname=cache_dirname, clear_content=clear_content,
                                           load_content=load_content, max_entries=max_entries,
                                           commit_intermediate_results=commit_intermediate_results,
                                           download_backup=download_backup, verbose=verbose,
                                           clear_requests_cache=clear_requests_cache, download_request_backup=download_request_backup)



    def load_content(self):
        """
        Downloads all medatata from array exrpess on their samples and experiments. The metadata
        is saved as the text file. Within the text files, the data is stored as a JSON object.
        Args:
            start_year (:obj:`int`, optional): the first year to retrieve experiments for
            end_year (:obj:`int`, optional): the last year to retrieve experiments for
        """

        session = self.session

        # download and parse experiment ids

        path_to_file = "/home/yosef/Desktop/CDNA_FILES/Streptococcus_pneumoniae_r6.ASM704v1.cdna.all.fa"
        with open(path_to_file, mode='r') as handle:

            # Use Biopython's parse function to process individual
            # FASTA records (thus reducing memory footprint)
            for record in SeqIO.parse(handle, 'fasta'):

                # Extract individual parts of the FASTA record

                #print(record.gene)
                new_gene = GeneEntry(name=record.name)
                session.add(new_gene)
                #session.commit()
                #print(record.name)
                #print(record.id)
                #print(record.description)



        self.session.commit()




    def get_or_create_object(self, cls, **kwargs):
        """ Get the first instance of :obj:`cls` that has the property-values pairs described by kwargs, or create an instance of :obj:`cls`
        if there is no instance with the property-values pairs described by kwargs
        Args:
            cls (:obj:`class`): type of object to find or create
            **kwargs: values of the properties of the object
        Returns:
            :obj:`Base`: instance of :obj:`cls` hat has the property-values pairs described by kwargs
        """
        q = self.session.query(cls).filter_by(**kwargs)
        if self.session.query(q.exists()).scalar():
            return q.first()

        obj = cls(**kwargs)
        self.session.add(obj)
        return obj

if __name__ == '__main__':
    GetGenes().load_content()
