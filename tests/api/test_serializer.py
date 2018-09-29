
from kinetic_datanator.api.serializer import *
from kinetic_datanator.api.query import metabolite_concentrations
from kinetic_datanator.core import models, common_schema
import tempfile
import shutil
import unittest

class TestSerializers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.cache_dirname = tempfile.mkdtemp()
        flk = common_schema.CommonSchema(cache_dirname=cls.cache_dirname)
        cls.proline = flk.session.query(models.Metabolite).filter_by(metabolite_name = 'L-Proline').first()

        q = metabolite_concentrations.MetaboliteConcentrationQuery(cache_dirname=cls.cache_dirname, include_variants=True)
        cls.obs = q.run(cls.proline)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.cache_dirname)

    def test_observed_value(self):
        serialized = ObservedValueSerializer().dump(self.obs.observed_results, many=True)
