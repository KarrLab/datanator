import unittest
from datanator.schema_2 import transform_metabolite_concentrations
from datanator_query_python.config import config
from datanator_query_python.util import mongo_util


class TestTMC(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        conf = config.DatanatorTest()
        cls.src = transform_metabolite_concentrations.TMC(MongoDB=conf.SERVER,
                                                          db="test",
                                                          username=conf.USERNAME,
                                                          password=conf.PASSWORD,
                                                          max_entries=20,
                                                          des_col="observation")
        cls.obj = {
            "inchikey": "XJLXINKUBYWONI-NNYOXOHSSA-O",
            "chebi_id": "18009",
            "concentrations": [
                {
                    "affinities": [
                        {
                            "ec_number": "1.1.1.10",
                            "k_m": 0.000003,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.105",
                            "k_m": 0.000027,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.145",
                            "k_m": 3e-7,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.146",
                            "k_m": 0.0000010809606516578,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.284",
                            "k_m": 0.00232701525564402,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.300",
                            "k_m": 0.00000397287770711798,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.315",
                            "k_m": 6.9e-7,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.40",
                            "k_m": 0.0000042386962268565,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.42",
                            "k_m": 0.00000455192267069642,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.44",
                            "k_m": 0.00003,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.49",
                            "k_m": 0.0000238522852342474,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.50",
                            "k_m": 2.3e-7,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.51",
                            "k_m": 0.0096,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.71",
                            "k_m": 0.00065,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.87",
                            "k_m": 0.000027,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.48",
                            "k_m": 0.0087,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.5",
                            "k_m": 0.0014,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.15",
                            "k_m": 0.000352,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.2",
                            "k_m": 0.00196356818063443,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.5",
                            "k_m": 0.000008,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.2.2.6",
                            "k_m": 0.000065,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.16.1.8",
                            "k_m": 0.000036,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.6.5.10",
                            "k_m": 0.000057,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.8.1.7",
                            "k_m": 0.0000307064187926607,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.10",
                            "k_m": 0.000003,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.105",
                            "k_m": 0.000027,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.145",
                            "k_m": 3e-7,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.146",
                            "k_m": 0.0000010809606516578,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.284",
                            "k_m": 0.00232701525564402,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.300",
                            "k_m": 0.00000397287770711798,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.315",
                            "k_m": 6.9e-7,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.40",
                            "k_m": 0.0000042386962268565,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.42",
                            "k_m": 0.00000455192267069642,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.44",
                            "k_m": 0.00003,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.49",
                            "k_m": 0.0000238522852342474,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.50",
                            "k_m": 2.3e-7,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.51",
                            "k_m": 0.0096,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.71",
                            "k_m": 0.00065,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.87",
                            "k_m": 0.000027,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.48",
                            "k_m": 0.0087,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.5",
                            "k_m": 0.0014,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.15",
                            "k_m": 0.000352,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.2",
                            "k_m": 0.00196356818063443,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.5",
                            "k_m": 0.000008,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.2.2.6",
                            "k_m": 0.000065,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.16.1.8",
                            "k_m": 0.000036,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.6.5.10",
                            "k_m": 0.000057,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.8.1.7",
                            "k_m": 0.0000307064187926607,
                            "k_i": None,
                            "unit": "M"
                        }
                    ],
                    "concentration": 0.000028417783438096,
                    "species_name": "Mus musculus",
                    "ncbi_taxonomy_id": 10090,
                    "reference": {
                        "namespace": "doi",
                        "id": "10.1038/nchembio.2077"
                    },
                    "unit": "M",
                    "canon_anc_ids": [
                        131567,
                        2759,
                        33208,
                        7711,
                        40674,
                        9989,
                        10066,
                        10088
                    ],
                    "canon_anc_names": [
                        "cellular organisms",
                        "Eukaryota",
                        "Metazoa",
                        "Chordata",
                        "Mammalia",
                        "Rodentia",
                        "Muridae",
                        "Mus"
                    ]
                },
                {
                    "affinities": [
                        {
                            "ec_number": "1.1.1.2",
                            "k_m": 0.0000281602556806574,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.276",
                            "k_m": 0.0005,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.42",
                            "k_m": 0.0000204767251107922,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.44",
                            "k_m": 0.0000307589086065726,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.65",
                            "k_m": 0.00002,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.71",
                            "k_m": 0.0000202,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.11",
                            "k_m": 0.000036,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.4",
                            "k_m": 0.000014,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.5",
                            "k_m": 0.000420184504240796,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.4.1.4",
                            "k_m": 0.0000260933527794003,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.10",
                            "k_m": 0.00022,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.5",
                            "k_m": 0.000275499546279118,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.3",
                            "k_m": 0.0003,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.7",
                            "k_m": 0.00204939015319192,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.1.3.84",
                            "k_m": 0.002,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.2",
                            "k_m": 0.0000281602556806574,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.276",
                            "k_m": 0.0005,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.42",
                            "k_m": 0.0000204767251107922,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.44",
                            "k_m": 0.0000307589086065726,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.65",
                            "k_m": 0.00002,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.71",
                            "k_m": 0.0000202,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.11",
                            "k_m": 0.000036,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.4",
                            "k_m": 0.000014,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.5",
                            "k_m": 0.000420184504240796,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.4.1.4",
                            "k_m": 0.0000260933527794003,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.10",
                            "k_m": 0.00022,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.5",
                            "k_m": 0.000275499546279118,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.3",
                            "k_m": 0.0003,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.7",
                            "k_m": 0.00204939015319192,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.1.3.84",
                            "k_m": 0.002,
                            "k_i": None,
                            "unit": "M"
                        }
                    ],
                    "concentration": 0.000182814523452369,
                    "species_name": "Saccharomyces cerevisiae",
                    "ncbi_taxonomy_id": 4932,
                    "reference": {
                        "namespace": "doi",
                        "id": "10.1038/nchembio.2077"
                    },
                    "unit": "M",
                    "canon_anc_ids": [
                        131567,
                        2759,
                        4751,
                        4890,
                        4891,
                        4892,
                        4893,
                        4930
                    ],
                    "canon_anc_names": [
                        "cellular organisms",
                        "Eukaryota",
                        "Fungi",
                        "Ascomycota",
                        "Saccharomycetes",
                        "Saccharomycetales",
                        "Saccharomycetaceae",
                        "Saccharomyces"
                    ]
                },
                {
                    "affinities": [
                        {
                            "ec_number": "1.1.1.169",
                            "k_m": 0.000007,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.25",
                            "k_m": 0.000140945974641298,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.276",
                            "k_m": 0.00054,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.282",
                            "k_m": 0.000223606797749979,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.3",
                            "k_m": 0.000073,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.40",
                            "k_m": 0.0000415,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.41",
                            "k_m": 0.000017,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.42",
                            "k_m": 0.0000392,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.44",
                            "k_m": 0.000049,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.49",
                            "k_m": 0.000012818610191887,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.54",
                            "k_m": 0.00004,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.17.1.7",
                            "k_m": 0.000056,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.11",
                            "k_m": 0.00009,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.21",
                            "k_m": 0.0000985833593944261,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.38",
                            "k_m": 0.00017,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.39",
                            "k_m": 0.00022,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.4",
                            "k_m": 0.0001,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.41",
                            "k_m": 0.0000948683298050514,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.60",
                            "k_m": 0.0036,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.79",
                            "k_m": 0.0000446,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.8",
                            "k_m": 0.000428832719244705,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.4.1.4",
                            "k_m": 0.000149631169487067,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.1.3.7",
                            "k_m": 0.0027,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.1.4.16",
                            "k_m": 0.00015,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.1.4.37",
                            "k_m": 0.00015,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.271",
                            "k_m": 0.000069,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.94",
                            "k_m": 0.00025,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.17.1.8",
                            "k_m": 0.024,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.30",
                            "k_m": 0.0000139,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.6.1.2",
                            "k_m": 0.00017,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.8.1.9",
                            "k_m": 0.000015,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.169",
                            "k_m": 0.000007,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.25",
                            "k_m": 0.000140945974641298,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.276",
                            "k_m": 0.00054,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.282",
                            "k_m": 0.000223606797749979,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.3",
                            "k_m": 0.000073,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.40",
                            "k_m": 0.0000415,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.41",
                            "k_m": 0.000017,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.42",
                            "k_m": 0.0000392,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.44",
                            "k_m": 0.000049,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.49",
                            "k_m": 0.000012818610191887,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.54",
                            "k_m": 0.00004,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.17.1.7",
                            "k_m": 0.000056,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.11",
                            "k_m": 0.00009,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.21",
                            "k_m": 0.0000985833593944261,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.38",
                            "k_m": 0.00017,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.39",
                            "k_m": 0.00022,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.4",
                            "k_m": 0.0001,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.41",
                            "k_m": 0.0000948683298050514,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.60",
                            "k_m": 0.0036,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.79",
                            "k_m": 0.0000446,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.2.1.8",
                            "k_m": 0.000428832719244705,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.4.1.4",
                            "k_m": 0.000149631169487067,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.1.3.7",
                            "k_m": 0.0027,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.1.4.16",
                            "k_m": 0.00015,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "3.1.4.37",
                            "k_m": 0.00015,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.271",
                            "k_m": 0.000069,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.1.1.94",
                            "k_m": 0.00025,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.17.1.8",
                            "k_m": 0.024,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.5.1.30",
                            "k_m": 0.0000139,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.6.1.2",
                            "k_m": 0.00017,
                            "k_i": None,
                            "unit": "M"
                        },
                        {
                            "ec_number": "1.8.1.9",
                            "k_m": 0.000015,
                            "k_i": None,
                            "unit": "M"
                        }
                    ],
                    "concentration": 0.00000208,
                    "species_name": "Escherichia coli",
                    "ncbi_taxonomy_id": 562,
                    "reference": {
                        "namespace": "doi",
                        "id": "10.1038/nchembio.2077"
                    },
                    "unit": "M",
                    "canon_anc_ids": [
                        131567,
                        2,
                        1224,
                        1236,
                        91347,
                        543,
                        561
                    ],
                    "canon_anc_names": [
                        "cellular organisms",
                        "Bacteria",
                        "Proteobacteria",
                        "Gammaproteobacteria",
                        "Enterobacterales",
                        "Enterobacteriaceae",
                        "Escherichia"
                    ]
                },
                {
                    "ncbi_taxonomy_id": 562,
                    "species_name": "Escherichia coli",
                    "growth_media": "0.2 g/L NH4Cl, 2.0 g/L (NH4)2SO4, 3.25 g/L KH2PO4, 2.5 g/L K2HPO4, 1.5 g/L NaH2PO4, 0.5 g/L MgSO4; trace substances: 10 mg/L CaCl2, 0.5 mg/L ZnSO4, 0.25 mg/L CuCl2, 0.25 mg/L  MnSO4, 0.175 mg/L CoCl2, 0.125 mg/L H3BO3, 2.5 mg/L AlCl3, 0.5 mg/L Na2MoO4, 10",
                    "growth_system": "Bioreactor, pH controlled, aerated, dilution rate=0.125 L/h",
                    "concentration": "80.0",
                    "concentration_units": "uM",
                    "internal": None,
                    "error": "3.0",
                    "temperature": "37 oC",
                    "strain": "K12",
                    "growth_status": "Stationary Phase, glucose limited",
                    "molecules": "320000",
                    "molecules_error": "12000",
                    "reference": {
                        "namespace": "pubmed",
                        "id": "11488613",
                        "text": "Buchholz, A., Takors, R., Wandrey, C. (2001). \"Quantification of intracellular metabolites in Escherichia coli K12 using liquid chromatographic-electrospray ionization tandem mass spectrometric techniques.\" Anal Biochem 295:129-137."
                    },
                    "canon_anc_ids": [
                        131567,
                        2,
                        1224,
                        1236,
                        91347,
                        543,
                        561
                    ],
                    "canon_anc_names": [
                        "cellular organisms",
                        "Bacteria",
                        "Proteobacteria",
                        "Gammaproteobacteria",
                        "Enterobacterales",
                        "Enterobacteriaceae",
                        "Escherichia"
                    ]
                },
                {
                    "ncbi_taxonomy_id": 562,
                    "species_name": "Escherichia coli",
                    "growth_media": "Gutnick minimal complete medium (4.7 g/L KH2PO4; 13.5 g/L K2HPO4; 1 g/L K2SO4; 0.1 g/L MgSO4-7H2O; 10 mM NH4Cl) with 4 g/L glucose",
                    "growth_system": "Shake flask and filter culture",
                    "concentration": "2.08",
                    "concentration_units": "uM",
                    "internal": None,
                    "error": "0.0",
                    "temperature": "37 oC",
                    "strain": "K12 NCM3722",
                    "growth_status": "Mid-Log Phase",
                    "molecules": "8320",
                    "molecules_error": "0",
                    "reference": {
                        "namespace": "pubmed",
                        "id": "19561621",
                        "text": "Bennett, B. D., Kimball, E. H., Gao, M., Osterhout, R., Van Dien, S. J., Rabinowitz, J. D. (2009). \"Absolute metabolite concentrations and implied enzyme active site occupancy in Escherichia coli.\" Nat Chem Biol 5:593-599."
                    },
                    "canon_anc_ids": [
                        131567,
                        2,
                        1224,
                        1236,
                        91347,
                        543,
                        561
                    ],
                    "canon_anc_names": [
                        "cellular organisms",
                        "Bacteria",
                        "Proteobacteria",
                        "Gammaproteobacteria",
                        "Enterobacterales",
                        "Enterobacteriaceae",
                        "Escherichia"
                    ]
                },
                {
                    "ncbi_taxonomy_id": 562,
                    "species_name": "Escherichia coli",
                    "growth_media": "48 mM Na2HPO4, 22 mM KH2PO4, 10 mM NaCl, 45 mM (NH4)2SO4, supplemented with 1 mM MgSO4, 1 mg/l thiamine·HCl, 5.6 mg/l CaCl2, 8 mg/l FeCl3, 1 mg/l MnCl2·4H2O, 1.7 mg/l ZnCl2, 0.43 mg/l CuCl2·2H2O, 0.6 mg/l CoCl2·2H2O and 0.6 mg/l Na2MoO4·2H2O.  4 g/L Gluco",
                    "growth_system": "Bioreactor, pH controlled, O2 and CO2 controlled, dilution rate: 0.2/h",
                    "concentration": "110.0",
                    "concentration_units": "uM",
                    "internal": None,
                    "error": "0.0",
                    "temperature": "37 oC",
                    "strain": "BW25113",
                    "growth_status": "Stationary Phase, glucose limited",
                    "molecules": "440000",
                    "molecules_error": "0",
                    "reference": {
                        "namespace": "pubmed",
                        "id": "17379776",
                        "text": "Ishii, N., Nakahigashi, K., Baba, T., Robert, M., Soga, T., Kanai, A., Hirasawa, T., Naba, M., Hirai, K., Hoque, A., Ho, P. Y., Kakazu, Y., Sugawara, K., Igarashi, S., Harada, S., Masuda, T., Sugiyama, N., Togashi, T., Hasegawa, M., Takai, Y., Yugi, K., Arakawa, K., Iwata, N., Toya, Y., Nakayama, Y., Nishioka, T., Shimizu, K., Mori, H., Tomita, M. (2007). \"Multiple high-throughput analyses monitor the response of E. coli to perturbations.\" Science 316:593-597."
                    },
                    "canon_anc_ids": [
                        131567,
                        2,
                        1224,
                        1236,
                        91347,
                        543,
                        561
                    ],
                    "canon_anc_names": [
                        "cellular organisms",
                        "Bacteria",
                        "Proteobacteria",
                        "Gammaproteobacteria",
                        "Enterobacterales",
                        "Enterobacteriaceae",
                        "Escherichia"
                    ]
                },
                {
                    "ncbi_taxonomy_id": 562,
                    "species_name": "Escherichia coli",
                    "growth_media": None,
                    "growth_system": None,
                    "concentration": "630.0",
                    "concentration_units": "uM",
                    "internal": None,
                    "error": "0.0",
                    "temperature": None,
                    "strain": "K-12",
                    "growth_status": None,
                    "molecules": "2520000",
                    "molecules_error": "0",
                    "reference": {
                        "namespace": "pubmed",
                        "id": None,
                        "text": "1. Cybercell Database: <a href='http://ccdb.wishartlab.com/CCDB/cgi-bin/STAT_NEW.cgi'>http://ccdb.wishartlab.com/CCDB/cgi-bin/STAT_NEW.cgi</a> <br>\n\t2. Phillips R., Kondev, J., Theriot, J. (2008) “Physical Biology of the Cell” Garland Science, New York, NY."
                    },
                    "canon_anc_ids": [
                        131567,
                        2,
                        1224,
                        1236,
                        91347,
                        543,
                        561
                    ],
                    "canon_anc_names": [
                        "cellular organisms",
                        "Bacteria",
                        "Proteobacteria",
                        "Gammaproteobacteria",
                        "Enterobacterales",
                        "Enterobacteriaceae",
                        "Escherichia"
                    ]
                },
                {
                    "ncbi_taxonomy_id": 4932,
                    "species_name": "Saccharomyces cerevisiae",
                    "growth_media": "Minimal medium  supplemented with ammonia salts and glucose",
                    "concentration": "325.0",
                    "concentration_units": "&#181;M",
                    "internal": "false",
                    "error": "305.0",
                    "strain": None,
                    "reference": {
                        "namespace": "pubmed",
                        "id": "4578278",
                        "text": "Gancedo, J. M., Gancedo, C. (1973). \"Concentrations of intermediary metabolites in yeast.\" Biochimie 55:205-211."
                    },
                    "canon_anc_ids": [
                        131567,
                        2759,
                        4751,
                        4890,
                        4891,
                        4892,
                        4893,
                        4930
                    ],
                    "canon_anc_names": [
                        "cellular organisms",
                        "Eukaryota",
                        "Fungi",
                        "Ascomycota",
                        "Saccharomycetes",
                        "Saccharomycetales",
                        "Saccharomycetaceae",
                        "Saccharomyces"
                    ]
                },
                {
                    "ncbi_taxonomy_id": 4932,
                    "species_name": "Saccharomyces cerevisiae",
                    "growth_media": "Minimal medium supplemented with ammonia salts and (glucose or galactose)",
                    "concentration": "85.0",
                    "concentration_units": "&#181;M",
                    "internal": "false",
                    "error": "65.0",
                    "strain": None,
                    "reference": {
                        "namespace": "pubmed",
                        "id": "4578278",
                        "text": "Gancedo, J. M., Gancedo, C. (1973). \"Concentrations of intermediary metabolites in yeast.\" Biochimie 55:205-211."
                    },
                    "canon_anc_ids": [
                        131567,
                        2759,
                        4751,
                        4890,
                        4891,
                        4892,
                        4893,
                        4930
                    ],
                    "canon_anc_names": [
                        "cellular organisms",
                        "Eukaryota",
                        "Fungi",
                        "Ascomycota",
                        "Saccharomycetes",
                        "Saccharomycetales",
                        "Saccharomycetaceae",
                        "Saccharomyces"
                    ]
                }
            ],
            "kegg_id": "C00006",
            "metabolite": "NADP",
            "schema_version": "2",
            "synonyms": [
                "Adenine-nicotinamide dinucleotide phosphate",
                "b-NADP",
                "b-Nicotinamide adenine dinucleotide phosphate",
                "b-TPN",
                "beta-NADP",
                "beta-nicotinamide adenine dinucleotide phosphate",
                "beta-TPN",
                "Codehydrase II",
                "Codehydrogenase II",
                "Coenzyme II",
                "Cozymase II",
                "NAD phosphate",
                "nadp",
                "nadp+",
                "NAP",
                "Nicotinamide adenine dinucleotide phosphate",
                "Nicotinamide-adenine dinucleotide phosphate",
                "oxidized nicotinamide-adenine dinucleotide phosphate",
                "TPN",
                "Triphosphopyridine nucleotide"
            ]
        }

    @classmethod
    def tearDownClass(cls):
        pass

    @unittest.skip('for now')
    def test_build_conc_obs(self):
        a = self.src.build_conc_observation(self.obj)
        for x  in a:
            print(x[0])

    def test_process_docs(self):
        self.src.process_docs()

