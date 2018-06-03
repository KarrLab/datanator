from kinetic_datanator.data_source import refseq
import os
from os import path
from Bio import SeqIO
import openpyxl
from kinetic_datanator.core import flask_common_schema, models
import json
import sqlalchemy_utils
import sqlalchemy


CACHE_DIRNAME = os.path.join(os.path.dirname(__file__), '..', 'data_source', 'cache')

class Uploader():
    def __init__(self, parent_node, dict_of_new_objects):
        self.parent_node = parent_node
        self.dict_of_new_objects = dict_of_new_objects

class UploadData():
    def __init__(self, cache_dirname=CACHE_DIRNAME):
        self.cache_dirname = cache_dirname
        self.flask = flask_common_schema.FlaskCommonSchema(cache_dirname=cache_dirname)


    def upload_reference_genome(self, path_to_annotation_file):
        #bio_seqio_object = SeqIO.parse(path_to_annotation_file, "genbank")
        refseq.Refseq(cache_dirname=self.cache_dirname).load_content([SeqIO.parse(path_to_annotation_file, "genbank")])


    def upload_rna_rseq_experiment(self, path_to_template_directory):
        print(path_to_template_directory)
        filename = "{}/RNA-SeqMetadataTemplate.xlsx".format(path_to_template_directory)
        wb = openpyxl.load_workbook(filename=filename)
        ws = wb.get_sheet_by_name('Experiments')
        for i in range(2, ws.max_row+1):
            experiment_name = wb.get_sheet_by_name('Experiments').cell(row=i,column=1).value
            experiment_description = wb.get_sheet_by_name('Experiments').cell(row=i,column=2).value
            new_experiment = self.flask.get_or_create_object(models.RNASeqExperiment,
                exp_name=experiment_name,
                accession_number=experiment_name)
            print(experiment_name)
            print(experiment_description)
            self.flask.session.add(new_experiment)
        self.flask.session.commit()

    def get_or_create_object_upload(self, session, cls, kwargs):
        """ Get the first instance of :obj:`cls` that has the property-values pairs described by kwargs, or create an instance of :obj:`cls`
        if there is no instance with the property-values pairs described by kwargs
        Args:
            cls (:obj:`class`): type of object to find or create
            **kwargs: values of the properties of the object
        Returns:
            :obj:`Base`: instance of :obj:`cls` hat has the property-values pairs described by kwargs
        """
        q = session.query(cls).filter_by(**kwargs)
        if session.query(q.exists()).scalar():
            return q.first()

        obj = cls(**kwargs)
        session.add(obj)
        return obj

    def get_primitive_fields(self, object, object_type):
        print(object_type)
        list_of_primitives = {}
        #new_to upload = []
        for field in object:
            if (type(object[field]) is not list) and (type(object[field]) is not dict):
                list_of_primitives[field] = object[field]
            if type(object[field]) is dict:
                #print(object[field])
                #new_type = getattr(object_type, field)
                new_type = sqlalchemy_utils.get_type(object_type.__dict__[field])
                #print(new_type)
                #print("above")
                new_object = self.upload_from_json(object[field], new_type)
                list_of_primitives[field] = new_object

        return list_of_primitives

    def get_dict_of_new_objects(self, object):
        dict_of_new_objects = {}
        for field in object:
            if type(object[field]) is list:
                dict_of_new_objects[field] = object[field]
        return dict_of_new_objects

    def fill_out_total_fields(self, new_object, a_dict):
        return new_object
        for key in self.get_dict_of_new_objects(a_dict):
            new_object.__dict__[key] = []
        return new_object




    def upload_from_json(self, objects, first_data_type):
        session = self.flask.session
        #objects = json.loads(a_json)
        #print(objects)
        #print(type(first_data_type))
        if type(first_data_type) is str:
            print("blueberry muffins")
            models_object = models.__dict__[first_data_type]
        else:
            models_object = first_data_type


        #print(models_object)
        #print(models.RNASeqExperiment)

        data_type = models_object

        #while True:
        #kwargs = {"accession_number":"blue"}
        to_upload = []
        #list_of_fields = {}
        dict_of_new_objects = {}
        print(data_type)
        list_of_fields = self.get_primitive_fields(objects, data_type)
        print(list_of_fields)
        dict_of_new_objects = self.get_dict_of_new_objects(objects)
        """
        for field in objects:
            if type(objects[field]) is not list:
                list_of_fields[field] = objects[field]
            if type(objects[field]) is list:
                dict_of_new_objects[field] = objects[field]
        """

        
        first_parent = self.get_or_create_object_upload(session, models_object, list_of_fields)
        #print(first_parent)
        #new_experiment = self.flask.get_or_create_object(models.RNASeqExperiment, accession_number="E-MTAB-3")
        #new_experiment.samples.append(self.flask.get_or_create_object(models.RNASeqDataSet, name="blue"))
        #new_experiment.__dict__["samples"].append(self.flask.get_or_create_object(models.RNASeqDataSet, name="blue"))
        #new_sample = self.flask.get_or_create_object(models.RNASeqDataSet, name="blue")
        #getattr(new_experiment, "samples").append(new_sample)
        #print(new_experiment.samples[0].name)

        #print(new_experiment.samples.append("blue"))
        #print(new_experiment.__dict__['samples'])
        """
        for key in self.get_dict_of_new_objects(objects):
            new_experiment.__dict__[key] = []
        """
        #print(new_experiment.__dict__)
        #exp = session.query(models.RNASeqExperiment).filter_by(accession_number='E-MTAB-3').first()



        
        first_upload = Uploader(parent_node=first_parent, dict_of_new_objects=dict_of_new_objects)
        to_upload.append(first_upload)
        #print(to_upload[0].parent_node.__dict__)
        #session.commit()
        #print(to_upload[0].parent_node.__dict__)


        while to_upload:
            new_to_upload = []
            for uploader in to_upload:
                for key in uploader.dict_of_new_objects:
                    #print(uploader.parent_node)
                    #print(uploader.parent_node.__dict__)
                    #print(type(uploader.parent_node).__dict__)
                    new_obj_type = sqlalchemy_utils.get_type(type(uploader.parent_node).__dict__[key])
                    list_of_new_objects = uploader.dict_of_new_objects[key]
                    for new_object in list_of_new_objects:
                        the_new_object =self.get_or_create_object_upload(session, new_obj_type, self.get_primitive_fields(new_object, new_obj_type))
                        #print(uploader.parent_node.__dict__)
                        #print(uploader.parent_node._experimentmetadata)
                        #the_new_object = self.fill_out_total_fields(the_new_object, new_object)
                        #uploader.parent_node.__dict__[key].append(the_new_object)
                        #metadata = self.flask.get_or_create_object(models.ExperimentMetadata, description="green")
                        #print(type(uploader.parent_node._experimentmetadata) is sqlalchemy.orm.collections.InstrumentedList)
                        #print(type(getattr(uploader.parent_node, "samples")) is sqlalchemy.orm.collections.InstrumentedList)
                        #setattr(uploader.parent_node, "_experimentmetadata", metadata)
                        #print(uploader.parent_node)
                        #print(key)
                        getattr(uploader.parent_node, key).append(the_new_object)
                        dict_of_new_objects = self.get_dict_of_new_objects(new_object)
                        if dict_of_new_objects:
                            new_uploader = Uploader(parent_node = the_new_object, dict_of_new_objects = dict_of_new_objects)
                            new_to_upload.append(new_uploader)
            to_upload = new_to_upload
        session.commit()
        return first_parent

        
        





        """
        for 






        session.commit()
        print(new_experiment)
        exp = session.query(models.RNASeqExperiment).filter_by(accession_number='E-MTAB-3').first()
        print(exp.accession_number)
            #for field in 

        #print(models_object.__dict__)
        
        #new_object = models_object(**objects)
        new_type = sqlalchemy_utils.get_type(models_object.__dict__['samples'])
        new_sample = new_type()
        """












        #print(new_sample.__dict__)
        #print(models_object._sa_class_manager).__dict__
        #print(type(models_object.samples))
        


































    def upload_processed_data(sample_name, csv_file, top_directory=path.dirname(__file__)):
        #this takes in a csv file, and puts the file in the proper place so that it can be queried. 
        new_pandas = pd.read_csv(csv_file, sep=',').set_index("gene_locus")
        new_pandas.to_pickle("{}/{}".format(top_directory, sample_name))




