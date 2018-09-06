
from kinetic_datanator.api.serializer import *
from kinetic_datanator.api.query import metabolite_concentrations
from kinetic_datanator.core import models, common_schema
import tempfile
import shutil
import unittest

class TestSerializers(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        flk = common_schema.CommonSchema(cache_dirname=self.cache_dirname)
        self.proline = flk.session.query(models.Metabolite).filter_by(metabolite_name = 'L-Proline').first()

        q = metabolite_concentrations.MetaboliteConcentrationQuery(cache_dirname=self.cache_dirname, include_variants=True)
        self.obs = q.run(self.proline)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)


    def test_observed_value(self):
        serialized = ObservedValueSerializer().dump(self.obs.observed_results, many=True)
        print(serialized)
