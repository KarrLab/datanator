
from capturer import CaptureOutput
from cement.utils import test
from kinetic_datanator.core import upload_data
from kinetic_datanator.util import warning_util
import os
import re
import shutil
import tempfile
import unittest
from kinetic_datanator.core import flask_common_schema, models
import json



warning_util.disable_warnings()


@unittest.skip("Can't run together with json test for some reason")
class TestUploadData(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        self.flk = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_upload_rna_rseq_experiment(self):
    	dir_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "fixtures/RNA-Seq_Experiment_Test_Template"))
    	upload_data.UploadData(self.cache_dirname).upload_rna_rseq_experiment(dir_path)
    	experiment = self.flk.session.query(models.RNASeqExperiment).filter_by(accession_number = 'test_experiment_1').first()
    	self.assertEqual(experiment.exp_name, 'test_experiment_1')
    	sample = self.flk.session.query(models.RNASeqDataSet).filter_by(experiment_accession_number = 'test_experiment_1', sample_name='test_sample_1').first()
        #self.assertEqual(sam

class TestUploadDataFromJson(unittest.TestCase):
    
    
    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        #self.cache_dirname = "/home/yosef/Desktop/angular"
        self.flk = flask_common_schema.FlaskCommonSchema(cache_dirname=self.cache_dirname)


    @classmethod
    def tearDownClass(self):
        #pass
        shutil.rmtree(self.cache_dirname)
    
    
    



    def test_upload_from_json(self):
        
        flk = self.flk
        session = flk.session
        

        test_json = """
             {
                "_experimentmetadata": [
                    {
                        "description": "a test of uploading",
                        "experiment_design": [
                            {
                                "name": "test_design"
                            }
                        ],
                        "experiment_type": [],
                        "method": [],
                        "resource": [],
                        "taxon": [
                            {
                                "_metadata": [],
                                "name": "test_organism"
                            }
                        ]
                    }
                ],
                "accession_number": "E-MTAB-3",
                "exp_name": "test_name",
                "has_fastq_files": true,
                "samples": [
                    {
                        "_metadata": [],
                        "assay": "test_assay",
                        "ensembl_organism_strain": "test_strain",
                        "experiment_accession_number": "E-MTAB-3",
                        "full_strain_specificity": false,
                        "read_type": "single",
                        "reference_genome": [],
                        "sample_name": "test_sample_name"
                    }
                ]
            }
        """
        test_json_simple = """
             {
             "_experimentmetadata": [
                    {
                        "description": "a test of uploading",
                        "experiment_design": [
                            {
                                "name": "test_design"
                            }
                        ],
                        "experiment_type": [],
                        "method": [],
                        "resource": [],
                        "taxon": [
                            {
                                "_metadata": [],
                                "name": "test_organism"
                            }
                        ]
                    }
                ],
                "accession_number": "E-MTAB-3",
                "exp_name": "test_name",
                "has_fastq_files": true,
                "samples": [
                    {
                        "assay": "test_assay",
                        "ensembl_organism_strain": "test_strain",
                        "experiment_accession_number": "E-MTAB-3",
                        "full_strain_specificity": false,
                        "read_type": "single",
                        "reference_genome": [],
                        "sample_name": "test_sample_name"
                    }
                ]
            }
            """

        test_updated_json = """
                {
            "_experimentmetadata": {
                "description": "a test of experiment metadata",
                "experiment_design": [
                    {
                        "name": "a new experiment design"
                    }
                ],
                "experiment_type": [],
                "method": [],
                "resource": [],
                "taxon": []
            },
            "accession_number": "E-MTAB-3",
            "exp_name": "test experiment name",
            "has_fastq_files": true,
            "samples": [
                {
                    "_metadata": {
                        "cell_compartment": [],
                        "cell_line": [],
                        "characteristic": [],
                        "conditions": [],
                        "method": [],
                        "synonym": [],
                        "taxon": [],
                        "variable": []
                    },
                    "assay": "test assay",
                    "ensembl_organism_strain": "test strain",
                    "experiment_accession_number": "E-MTAB-3",
                    "full_strain_specificity": false,
                    "read_type": "single",
                    "reference_genome": [],
                    "sample_name": "test sample name"
                }
            ]
        }
        """
        
        upload_data.UploadData(self.cache_dirname).upload_from_json(json.loads(test_updated_json), "RNASeqExperiment")
        #upload_data.UploadData().upload_from_json(json.loads(test_updated_json), "RNASeqExperiment")
        
        exp = session.query(models.RNASeqExperiment).filter_by(accession_number = "E-MTAB-3").first()
        self.assertEqual(exp.accession_number, "E-MTAB-3")
        self.assertEqual(exp.exp_name, "test experiment name")
        sample = exp.samples[0]
        self.assertEqual(sample.ensembl_organism_strain, "test strain")
        self.assertEqual(sample.full_strain_specificity, False)
        self.assertEqual(sample.sample_name, "test sample name")
        experiment_metadata = exp._experimentmetadata
        self.assertEqual(experiment_metadata.description, "a test of experiment metadata")
        self.assertEqual(experiment_metadata.experiment_design[0].name, "a new experiment design")

        