"""
import pprint
from Bio import SeqIO
import datetime
import dateutil.parser
import pkg_resources
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
from kinetic_datanator.core import data_source
"""
import sqlalchemy
import sqlalchemy.ext.declarative
import sqlalchemy.orm
from kinetic_datanator.core import data_source
from Bio import SeqIO
import Bio
import json
import math
Base = sqlalchemy.ext.declarative.declarative_base()

def create_orm(class_1, class_2):
    
    new_table = sqlalchemy.Table(
    "{}_{}".format(class_1,class_2), Base.metadata,
    sqlalchemy.Column('{}{}_id'.format(class_1[0].lower(), class_1[1:]), sqlalchemy.Integer, sqlalchemy.ForeignKey('{}{}._id'.format(class_1[0].lower(), class_1[1:])), index=True),
    sqlalchemy.Column('{}{}_id'.format(class_2[0].lower(), class_2[1:]), sqlalchemy.Integer, sqlalchemy.ForeignKey('{}{}._id'.format(class_2[0].lower(), class_2[1:])), index=True),
    )
    
    return sqlalchemy.orm.relationship(class_2, secondary=new_table, backref=sqlalchemy.orm.backref('{}{}'.format(class_1[0].lower(), class_1[1:])))


class Location(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    _end = sqlalchemy.Column(sqlalchemy.Integer())
    _start = sqlalchemy.Column(sqlalchemy.Integer())
    end = sqlalchemy.Column(sqlalchemy.Integer())
    nofuzzy_end = sqlalchemy.Column(sqlalchemy.Integer())
    nofuzzy_start = sqlalchemy.Column(sqlalchemy.Integer())
    ref = sqlalchemy.Column(sqlalchemy.String())
    ref_db = sqlalchemy.Column(sqlalchemy.String())
    start = sqlalchemy.Column(sqlalchemy.Integer())
    strand = sqlalchemy.Column(sqlalchemy.Integer())

    __tablename__ = 'location'

class Qualifier(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    key = sqlalchemy.Column(sqlalchemy.String())
    value = sqlalchemy.Column(sqlalchemy.String())

    __tablename__ = 'qualifier'

class GeneSynonym(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    name = sqlalchemy.Column(sqlalchemy.String())
    __tablename__ = 'geneSynonym'

class EcNumber(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    ec_number = sqlalchemy.Column(sqlalchemy.String())
    __tablename__ = 'ecNumber'

class Identifier(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    namespace = sqlalchemy.Column(sqlalchemy.String())
    name = sqlalchemy.Column(sqlalchemy.String())
    __tablename__ = 'identifier'


class Gene(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    id = sqlalchemy.Column(sqlalchemy.String())
    name = sqlalchemy.Column(sqlalchemy.String())
    locus_tag = sqlalchemy.Column(sqlalchemy.String())
    gene_synonyms = create_orm("Gene", "GeneSynonym")#sqlalchemy.orm.relationship('GeneSynonym', secondary=gene_gene_synonym, backref=sqlalchemy.orm.backref('genes'))
    location = create_orm("Gene", "Location")
    qualifiers = create_orm('Gene','Qualifier')
    ec_numbers = create_orm('Gene','EcNumber')
    identifiers = create_orm('Gene', 'Identifier')
    essentiality = sqlalchemy.Column(sqlalchemy.String())
    __tablename__ = 'gene'


class ReferenceGenome(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    accessions = create_orm("ReferenceGenome", "ReferenceGenomeAccession")
    #accession = sqlalchemy.Column(sqlalchemy.String())
    version = sqlalchemy.Column(sqlalchemy.String())
    organism = sqlalchemy.Column(sqlalchemy.String())
    genes = create_orm("ReferenceGenome", "Gene")#sqlalchemy.orm.relationship('GeneSynonym', secondary=gene_gene_synonym, backref=sqlalchemy.orm.backref('genes'))


    __tablename__ = 'referenceGenome'


class ReferenceGenomeAccession(Base):
    _id = sqlalchemy.Column(sqlalchemy.Integer(), primary_key=True)
    id = sqlalchemy.Column(sqlalchemy.String())

    __tablename__ = 'referenceGenomeAccession'

class Refseq(data_source.HttpDataSource):
    base_model = Base


    def __init__(self, name=None, cache_dirname=None, clear_content=False, load_content=False, max_entries=float('inf'),
                 commit_intermediate_results=False, download_backups=False, verbose=False,
                 clear_requests_cache=False, download_request_backup=False):

        super(Refseq, self).__init__(name=name, cache_dirname=cache_dirname, clear_content=clear_content,
                                           load_content=load_content, max_entries=max_entries,
                                           commit_intermediate_results=commit_intermediate_results,
                                           download_backups=download_backups, verbose=verbose,
                                            clear_requests_cache=clear_requests_cache, download_request_backup=download_request_backup)


    def load_content(self, list_bio_seqio_objects):
        session = self.session


        for bio_seqio_object in list_bio_seqio_objects:

            for seq_record in bio_seqio_object:

                #accession = seq_record.annotations['accessions'][0]
                ref_genome = self.get_or_create_object(ReferenceGenome,
                    version = seq_record.id
                    )

                if 'accessions' in seq_record.annotations:
                    for accession in seq_record.annotations['accessions']:
                        print(accession)
                        acc = self.get_or_create_object(ReferenceGenomeAccession,
                            id=accession)
                        ref_genome.accessions.append(acc)


                ref_genome.version = seq_record.id
                if 'organism' in seq_record.annotations:
                    ref_genome.organism = seq_record.annotations['organism']

                #session.add(ref_genome)
                for seq_feature in seq_record.features:
                    if seq_feature.type == 'CDS':

                        gene = Gene()
                        gene.id = seq_feature.id
                        
                        parts = []
                        if type(seq_feature.location)==Bio.SeqFeature.FeatureLocation:
                            parts.append(seq_feature.location)

                        elif type(seq_feature.location)==Bio.SeqFeature.CompoundLocation:
                            x = seq_feature
                            for part in seq_feature.location.parts:
                                parts.append(part)
                        for part in parts:
                            location = Location()
                            location._end = part._end
                            location._start = part._start
                            location.end = part.end
                            location.nofuzzy_end = part.nofuzzy_end
                            location.nofuzzy_start = part.nofuzzy_start
                            location.ref = part.ref
                            location.ref_db = part.ref_db
                            location.start = part.start
                            location.strand = part.strand
                            session.add(location)
                            gene.location.append(location)

                        qual = seq_feature.qualifiers
                        #print(seq_feature.qualifiers)
                        if 'gene' in qual:
                            gene.name = qual['gene'][0]
                            if len(qual['gene'])>1:
                                raise ValueError("More than one value")
                        if "gene_synonym" in qual:
                            for name in qual['gene_synonym']:
                                for embedded_name in name.split(";"):
                                    gene_synonym = self.get_or_create_object(GeneSynonym,
                                        name = embedded_name)
                                    gene.gene_synonyms.append(gene_synonym)
                        if 'EC_number' in qual:
                            for number in qual['EC_number']:
                                ec_number = self.get_or_create_object(EcNumber, 
                                    ec_number = number)
                                gene.ec_numbers.append(ec_number)
                        if 'db_xref' in qual:
                            for identifier in qual['db_xref']:
                                if len(identifier)>3: #make sure its not nan
                                    identifier = self.get_or_create_object(Identifier,
                                        namespace = identifier.split(":")[0],
                                        name = identifier.split(":")[1])
                                    gene.identifiers.append(identifier)
                        if 'locus_tag' in qual:
                            gene.locus_tag = qual['locus_tag'][0]
                            if len(qual['locus_tag'])>1:
                                raise ValueError("More than one value")
                        if 'essentiality2016_assigned' in qual:
                            gene.essentiality = qual['essentiality2016_assigned'][0]

                            
                        session.add(gene)
                        ref_genome.genes.append(gene)
                    elif seq_feature.type == 'source':
                        pass
                    else: 
                        pass

                        #to do: finish this up. its metadata on the reference genome

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

    filenames = [
            "/home/yosef/Desktop/Genome_Database/mpn_annotation.gbff",
            #"/home/yosef/Desktop/Genome_Database/sequence.gb",
            ]
    get_genes = Refseq(cache_dirname = '/home/yosef/Desktop/Genome_Database')
    get_genes.load_content(filenames)
