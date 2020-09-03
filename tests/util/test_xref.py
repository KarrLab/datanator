from datanator.util import x_ref
from datanator_query_python.config import config
import unittest


class TestTransform(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        conf = config.SchemaMigration()
        cls.des_col = "orthodb"
        cls.src = x_ref.XRef(MongoDB=conf.SERVER,
                            db="datanator",
                            des_col=cls.des_col,
                            username=conf.USERNAME,
                            password=conf.PASSWORD,
                            max_entries=10,
                            verbose=True)

    @classmethod
    def tearDownClass(cls):
        pass

    @unittest.skip("passed")
    def test_get_kegg_rxn(self):
        sub, s_n, pro, p_n = self.src.get_kegg_rxn("rn:R08549")
        self.assertEqual(sub, ["KPGXRSRHYNQIFN-UHFFFAOYSA-L", 
                               "RGJOEKWQDUBAIZ-IBOSZNHHSA-N", 
                               "BAWFJGJZGIEFAR-NNYOXOHSSA-O"])
        self.assertEqual(pro, ["VNOYUJKHFWYWIR-ITIYDSSPSA-N",
                               "CURLTUGMZLYLDI-UHFFFAOYSA-N",
                               "BOPGDPNILDQYTO-NNYOXOHSSA-N",
                               "GPRLSGONYQIRFK-UHFFFAOYSA-N"])
        self.assertEqual(p_n, ['Succinyl-CoA', 'CO2', 'NADH', 'H+'])
        sub, s_n, pro, p_n = self.src.get_kegg_rxn("rn:R00627")
        self.assertEqual(sub, ["ZSLZBFCDCINBPY-ZSJPKINUSA-N", "O[*]"])

    def test_uniprot_id_to_orthodb(self):
        cache = {"mock": "something"}
        r, _ = self.src.uniprot_id_to_orthodb("P12345")
        self.assertEqual(r, {'orthodb_id': '1104596at2759', 'orthodb_name': 'Pyridoxal phosphate-dependent transferase'})
        r, _ = self.src.uniprot_id_to_orthodb("mock", cache)
        self.assertEqual(r, "something")