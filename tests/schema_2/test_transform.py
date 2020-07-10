from datanator.schema_2 import transform
from datanator_query_python.config import config
import unittest
import numpy as np


class TestTransform(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        conf = config.SchemaMigration()
        cls.des_col = "transformation-test"
        cls.src = transform.Transform(MongoDB=conf.SERVER,
                                      db="test",
                                      des_col=cls.des_col,
                                      username=conf.USERNAME,
                                      password=conf.PASSWORD,
                                      max_entries=20,
                                      verbose=True)

    @classmethod
    def tearDownClass(cls):
        pass
        # cls.src.db_obj.drop_collection(cls.des_col)

    @unittest.skip("for now")
    def test_parse_docs(self):
        self.src.process_docs("uniprot")
    
    @unittest.skip("passed")
    def test_build_uniprot_entity(self):
        obj = {
                "uniprot_id": "Q75IW1",
                "add_id": [
                    {
                        "name_space": "gene_name_alt",
                        "value": None
                    },
                    {
                        "name_space": "gene_name_orf",
                        "value": "OsJ_11271 OSJNBb0059G13.19"
                    },
                    {
                        "name_space": "gene_name_oln",
                        "value": "Os03g0416300 LOC_Os03g30260"
                    }
                ],
                "ancestor_name": [
                    "cellular organisms",
                    "Eukaryota",
                    "Viridiplantae",
                    "Streptophyta",
                    "Streptophytina",
                    "Embryophyta",
                    "Tracheophyta",
                    "Euphyllophyta",
                    "Spermatophyta",
                    "Magnoliophyta",
                    "Mesangiospermae",
                    "Liliopsida",
                    "Petrosaviidae",
                    "commelinids",
                    "Poales",
                    "Poaceae",
                    "BOP clade",
                    "Oryzoideae",
                    "Oryzeae",
                    "Oryzinae",
                    "Oryza",
                    "Oryza sativa"
                ],
                "ancestor_taxon_id": [
                    131567,
                    2759,
                    33090,
                    35493,
                    131221,
                    3193,
                    58023,
                    78536,
                    58024,
                    3398,
                    1437183,
                    4447,
                    1437197,
                    4734,
                    38820,
                    4479,
                    359160,
                    147367,
                    147380,
                    1648021,
                    4527,
                    4530
                ],
                "canon_anc_ids": [
                    131567,
                    2759,
                    33090,
                    35493,
                    4447,
                    38820,
                    4479,
                    4527,
                    4530
                ],
                "canon_anc_names": [
                    "cellular organisms",
                    "Eukaryota",
                    "Viridiplantae",
                    "Streptophyta",
                    "Liliopsida",
                    "Poales",
                    "Poaceae",
                    "Oryza",
                    "Oryza sativa"
                ],
                "canonical_sequence": "MARFLLGAAAIALLAGVSSLLLMVPFAEAYDPLDPNGNITIKWDITQWTPDGYVAVVTIYNFQKYRHIQAPGWSLGWAWAKKEIIWSMAGGQATEQGDCSAFKANIPHCCKRDPRVVDLVPGAPYNMQFGNCCKGGVLTSWVQDPLNAVASFQITVGHSGTSNKTVKAPKNFTLKAPGPGYSCGLAQEVKPPTRFISLDGRRTTQAHVTWNVTCTYSQFVAQRAPTCCVSLSSFYNETIVNCPKCACGCQNKKPGSCVEGNSPYLASVVNGPGKGSLTPLVQCTPHMCPIRVHWHVKLNYRDYWRVKVTITNWNYRMNYSQWNLVVQHPNFENVSTVFSFNYKSLNPYGVINDTAMMWGVKYYNDLLMVAGPDGNVQSELLFRKDRSTFTFDKGWAFPRRIYFNGESCVMPSPDLYPWLPPSSTPRFRTVFLLMSFLVCGTLAFLHNHLVLDKNCGKC",
                "ec_number": None,
                "entrez_id": "4333115",
                "entry_name": "COBL2_ORYSJ",
                "gene_name": "BC1L2",
                "ko_name": [
                    None
                ],
                "ko_number": None,
                "length": 458,
                "mass": "51107",
                "ncbi_taxonomy_id": 39947,
                "protein_name": "COBRA-like protein 2 (Protein BRITTLE CULM1-like 2)",
                "schema_version": "2",
                "species_name": "Oryza sativa Japonica Group",
                "status": "reviewed",
                "abundances": [
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1447"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "119"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1219"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2443"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1984"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2883"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2984"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2595"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "389"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2373"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "3052"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2763"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1992"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1321"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "8918"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1730"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1463"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1730"
                    }
                ],
                "sabio_kinlaw_id": [
                    4573,
                    4574,
                    4575,
                    4576
                ],
                "modifications": [
                    {
                        "pro_id": "PR:000024921",
                        "uniprot_id": "P0AFG8-1",
                        "processing": "2-887",
                        "deletions": np.NaN,
                        "processsed_sequence_iubmb": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIYKLETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "processsed_formula": "C4436H6965N1217O1216S27",
                        "processsed_molecular_weight": 97668.439,
                        "processsed_charge": 98,
                        "modifications": "K --> MOD:00064 (716)",
                        "crosslinks": np.NaN,
                        "modified_sequence_abbreviated_bpforms": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIY{AA0055}LETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "modified_sequence_bpforms": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIY{AA0055}LETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "concrete": True,
                        "modified_formula": "C4438H6966N1217O1217S27",
                        "modified_molecular_weight": 97709.46800000001,
                        "modified_charge": 97,
                        "modifications_formula": "C2HO",
                        "modifications_molecular_weight": 41.028999999994994,
                        "modifications_charge": -1,
                        "pro_issues": np.NaN,
                        "monomeric_form_issues": np.NaN,
                        "reference": {
                            "doi": "10.1093/nar/gkw1075"
                        }
                    },
                    {
                        "pro_id": "PR:000036675",
                        "uniprot_id": "P0AFG8-1",
                        "processing": "2-887",
                        "deletions": np.NaN,
                        "processsed_sequence_iubmb": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIYKLETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "processsed_formula": "C4436H6965N1217O1216S27",
                        "processsed_molecular_weight": 97668.439,
                        "processsed_charge": 98,
                        "modifications": "K --> MOD:00064 (716)",
                        "crosslinks": np.NaN,
                        "modified_sequence_abbreviated_bpforms": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIY{AA0055}LETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "modified_sequence_bpforms": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIY{AA0055}LETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "concrete": True,
                        "modified_formula": "C4438H6966N1217O1217S27",
                        "modified_molecular_weight": 97709.46800000001,
                        "modified_charge": 97,
                        "modifications_formula": "C2HO",
                        "modifications_molecular_weight": 41.028999999994994,
                        "modifications_charge": -1,
                        "pro_issues": np.NaN,
                        "monomeric_form_issues": np.NaN,
                        "reference": {
                            "doi": "10.1093/nar/gkw1075"
                        }
                    }
                ]   
            }
        result = self.src.build_uniprot_entity(obj)

    @unittest.skip("passed")
    def test_build_uniprot_obs(self):
        obj = { 
                "uniprot_id": "Q75IW1",
                "add_id": [
                    {
                        "name_space": "gene_name_alt",
                        "value": None
                    },
                    {
                        "name_space": "gene_name_orf",
                        "value": "OsJ_11271 OSJNBb0059G13.19"
                    },
                    {
                        "name_space": "gene_name_oln",
                        "value": "Os03g0416300 LOC_Os03g30260"
                    }
                ],
                "ancestor_name": [
                    "cellular organisms",
                    "Eukaryota",
                    "Viridiplantae",
                    "Streptophyta",
                    "Streptophytina",
                    "Embryophyta",
                    "Tracheophyta",
                    "Euphyllophyta",
                    "Spermatophyta",
                    "Magnoliophyta",
                    "Mesangiospermae",
                    "Liliopsida",
                    "Petrosaviidae",
                    "commelinids",
                    "Poales",
                    "Poaceae",
                    "BOP clade",
                    "Oryzoideae",
                    "Oryzeae",
                    "Oryzinae",
                    "Oryza",
                    "Oryza sativa"
                ],
                "ancestor_taxon_id": [
                    131567,
                    2759,
                    33090,
                    35493,
                    131221,
                    3193,
                    58023,
                    78536,
                    58024,
                    3398,
                    1437183,
                    4447,
                    1437197,
                    4734,
                    38820,
                    4479,
                    359160,
                    147367,
                    147380,
                    1648021,
                    4527,
                    4530
                ],
                "canon_anc_ids": [
                    131567,
                    2759,
                    33090,
                    35493,
                    4447,
                    38820,
                    4479,
                    4527,
                    4530
                ],
                "canon_anc_names": [
                    "cellular organisms",
                    "Eukaryota",
                    "Viridiplantae",
                    "Streptophyta",
                    "Liliopsida",
                    "Poales",
                    "Poaceae",
                    "Oryza",
                    "Oryza sativa"
                ],
                "canonical_sequence": "MARFLLGAAAIALLAGVSSLLLMVPFAEAYDPLDPNGNITIKWDITQWTPDGYVAVVTIYNFQKYRHIQAPGWSLGWAWAKKEIIWSMAGGQATEQGDCSAFKANIPHCCKRDPRVVDLVPGAPYNMQFGNCCKGGVLTSWVQDPLNAVASFQITVGHSGTSNKTVKAPKNFTLKAPGPGYSCGLAQEVKPPTRFISLDGRRTTQAHVTWNVTCTYSQFVAQRAPTCCVSLSSFYNETIVNCPKCACGCQNKKPGSCVEGNSPYLASVVNGPGKGSLTPLVQCTPHMCPIRVHWHVKLNYRDYWRVKVTITNWNYRMNYSQWNLVVQHPNFENVSTVFSFNYKSLNPYGVINDTAMMWGVKYYNDLLMVAGPDGNVQSELLFRKDRSTFTFDKGWAFPRRIYFNGESCVMPSPDLYPWLPPSSTPRFRTVFLLMSFLVCGTLAFLHNHLVLDKNCGKC",
                "ec_number": None,
                "entrez_id": "4333115",
                "entry_name": "COBL2_ORYSJ",
                "gene_name": "BC1L2",
                "ko_name": [
                    None
                ],
                "ko_number": None,
                "length": 458,
                "mass": "51107",
                "ncbi_taxonomy_id": 39947,
                "protein_name": "COBRA-like protein 2 (Protein BRITTLE CULM1-like 2)",
                "schema_version": "2",
                "species_name": "Oryza sativa Japonica Group",
                "status": "reviewed",
                "abundances": [
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1447"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "119"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1219"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2443"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1984"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2883"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2984"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2595"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "389"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2373"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "3052"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "2763"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1992"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1321"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "8918"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1730"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1463"
                    },
                    {
                        "organ": "WHOLE_ORGANISM",
                        "abundance": "1730"
                    }
                ],
                "sabio_kinlaw_id": [
                    4573,
                    4574,
                    4575,
                    4576
                ],
                "modifications": [
                    {
                        "pro_id": "PR:000024921",
                        "uniprot_id": "P0AFG8-1",
                        "processing": "2-887",
                        "deletions": np.NaN,
                        "processsed_sequence_iubmb": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIYKLETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "processsed_formula": "C4436H6965N1217O1216S27",
                        "processsed_molecular_weight": 97668.439,
                        "processsed_charge": 98,
                        "modifications": "K --> MOD:00064 (716)",
                        "crosslinks": np.NaN,
                        "modified_sequence_abbreviated_bpforms": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIY{AA0055}LETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "modified_sequence_bpforms": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIY{AA0055}LETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "concrete": True,
                        "modified_formula": "C4438H6966N1217O1217S27",
                        "modified_molecular_weight": 97709.46800000001,
                        "modified_charge": 97,
                        "modifications_formula": "C2HO",
                        "modifications_molecular_weight": 41.028999999994994,
                        "modifications_charge": -1,
                        "pro_issues": np.NaN,
                        "monomeric_form_issues": np.NaN,
                        "reference": {
                            "doi": "10.1093/nar/gkw1075"
                        }
                    },
                    {
                        "pro_id": "PR:000036675",
                        "uniprot_id": "P0AFG8-1",
                        "processing": "2-887",
                        "deletions": np.NaN,
                        "processsed_sequence_iubmb": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIYKLETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "processsed_formula": "C4436H6965N1217O1216S27",
                        "processsed_molecular_weight": 97668.439,
                        "processsed_charge": 98,
                        "modifications": "K --> MOD:00064 (716)",
                        "crosslinks": np.NaN,
                        "modified_sequence_abbreviated_bpforms": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIY{AA0055}LETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "modified_sequence_bpforms": "SERFPNDVDPIETRDWLQAIESVIREEGVERAQYLIDQLLAEARKGGVNVAAGTGISNYINTIPVEEQPEYPGNLELERRIRSAIRWNAIMTVLRASKKDLELGGHMASFQSSATIYDVCFNHFFRARNEQDGGDLVYFQGHISPGVYARAFLEGRLTQEQLDNFRQEVHGNGLSSYPHPKLMPEFWQFPTVSMGLGPIGAIYQAKFLKYLEHRGLKDTSKQTVYAFLGDGEMDEPESKGAITIATREKLDNLVFVINCNLQRLDGPVTGNGKIINELEGIFEGAGWNVIKVMWGSRWDELLRKDTSGKLIQLMNETVDGDYQTFKSKDGAYVREHFFGKYPETAALVADWTDEQIWALNRGGHDPKKIYAAFKKAQETKGKATVILAHTIKGYGMGDAAEGKNIAHQVKKMNMDGVRHIRDRFNVPVSDADIEKLPYITFPEGSEEHTYLHAQRQKLHGYLPSRQPNFTEKLELPSLQDFGALLEEQSKEISTTIAFVRALNVMLKNKSIKDRLVPIIADEARTFGMEGLFRQIGIYSPNGQQYTPQDREQVAYYKEDEKGQILQEGINELGAGCSWLAAATSYSTNNLPMIPFYIYYSMFGFQRIGDLCWAAGDQQARGFLIGGTSGRTTLNGEGLQHEDGHSHIQSLTIPNCISYDPAYAYEVAVIMHDGLERMYGEKQENVYYYITTLNENYHMPAMPEGAEEGIRKGIY{AA0055}LETIEGSKGKVQLLGSGSILRHVREAAEILAKDYGVGSDVYSVTSFTELARDGQDCERWNMLHPLETPRVPYIAQVMNDAPAVASTDYMKLFAEQVRTYVPADDYRVLGTDGFGRSDSRENLRHHFEVDASYVVVAALGELAKRGEIDKKVVADAIAKFNIDADKVNPRLA",
                        "concrete": True,
                        "modified_formula": "C4438H6966N1217O1217S27",
                        "modified_molecular_weight": 97709.46800000001,
                        "modified_charge": 97,
                        "modifications_formula": "C2HO",
                        "modifications_molecular_weight": 41.028999999994994,
                        "modifications_charge": -1,
                        "pro_issues": np.NaN,
                        "monomeric_form_issues": np.NaN,
                        "reference": {
                            "doi": "10.1093/nar/gkw1075"
                        }
                    }
                ]
            }
        self.assertEqual(self.src.build_uniprot_observation({}), {})
        self.assertEqual(self.src.build_uniprot_observation(obj)["entity"]["type"], "protein")

    def test_build_rna_observation(self):
        obj = {
            "uniprot_id": "Q8TUR2",
            "halflives": [
                {
                    "halflife": 4041.87006,
                    "std": 523.16592,
                    "std_over_avg": 0.1294366004,
                    "unit": "s",
                    "reference": [
                        {
                            "doi": "10.1186/s12864-016-3219-8"
                        }
                    ],
                    "growth_medium": "TMA",
                    "ordered_locus_name": "MA0001",
                    "ar_cog": "arCOG00468",
                    "cog_class": "L",
                    "cog": "COG1474",
                    "species": "Methanosarcina acetivorans",
                    "ncbi_taxonomy_id": 188937
                },
                {
                    "systematic_name": "YJL194W",
                    "halflife": 1200,
                    "unit": "s",
                    "species": "Saccharomyces cerevisiae W303",
                    "ncbi_taxonomy_id": 580240,
                    "r_squared": 0.98,
                    "reference": [
                        {
                            "doi": "10.1091/mbc.e11-01-0028"
                        }
                    ]
                },
                {
                    "accession_id": "NM_031449 ",
                    "probeset_id": 3000010,
                    "values": [
                        {
                            "gm07029": 6.200316296854718,
                            "biological_replicates": "a1",
                            "note": "independent cell cultures for the same cell line",
                            "unit": "hr"
                        },

                        {
                            "gm07029": 5.817285876322001,
                            "biological_replicates": "a3",
                            "note": "independent cell cultures for the same cell line",
                            "unit": "hr"
                        },
                        {
                            "gm10835": 4.167696688892588,
                            "biological_replicates": "a1",
                            "note": "independent cell cultures for the same cell line",
                            "unit": "hr"
                        },
                        {
                            "gm10835": 4.454436766714646,
                            "biological_replicates": "a2",
                            "note": "independent cell cultures for the same cell line",
                            "unit": "hr"
                        },
                        {
                            "gm10835": 4.0912138205438024,
                            "biological_replicates": "a3",
                            "note": "independent cell cultures for the same cell line",
                            "unit": "hr"
                        },
                        {
                            "gm12813": 7.853596564318888,
                            "biological_replicates": "a1",
                            "note": "independent cell cultures for the same cell line",
                            "unit": "hr"
                        },
                        {
                            "gm12813": 8.231318451169917,
                            "technical_replicates": "a1",
                            "note": "separate RNA aliquots from the same cell culture",
                            "unit": "hr"
                        },
                        {
                            "gm12813": 7.958703606479381,
                            "biological_replicates": "a2",
                            "note": "independent cell cultures for the same cell line",
                            "unit": "hr"
                        },
                        {
                            "gm12813": 7.798393420876806,
                            "technical_replicates": "a2",
                            "note": "separate RNA aliquots from the same cell culture",
                            "unit": "hr"
                        },
                        {
                            "gm12813": 7.167623222693315,
                            "technical_replicates": "a3",
                            "note": "separate RNA aliquots from the same cell culture",
                            "unit": "hr"
                        },
                        {
                            "gm07019": 5.640622176,
                            "unit": "hr"
                        },
                        {
                            "gm12812": 6.162088116,
                            "unit": "hr"
                        },
                        {
                            "gm12814": 6.042021467,
                            "unit": "hr"
                        },
                        {
                            "gm12815": 6.758158592,
                            "unit": "hr"
                        }
                    ],
                    "anova_3": 3.923e-7,
                    "anova_7": 0.00000515675001947275,
                    "false_discovery_rate_3": 0.0004363744,
                    "false_discovery_rate_7": 0.00154363733333333,
                    "species": "Homo sapiens",
                    "ncbi_taxonomy_id": 9606,
                    "reference": [
                        {
                            "doi": "10.1038/srep01318"
                        }
                    ],
                    "gene_symbol": "ZMIZ2 "
                },
                {
                    "chromosome": "chr10",
                    "systematic_name": "YJL194W",
                    "gene_name": "CDC6",
                    "type": "coding_dna_sequences",
                    "halflife": 1593.6956414259646,
                    "unit": "s",
                    "species": "Saccharomyces cerevisiae S288C",
                    "ncbi_taxonomy_id": 559292,
                    "reference": [
                        {
                            "doi": "10.1016/j.cell.2013.12.026"
                        }
                    ]
                },
                {
                    "halflife": 353.6181058470822,
                    "r_sqaured": 0.9974211229873777,
                    "unit": "s",
                    "reference": [
                        {
                            "doi": "10.1093/nar/gks1019",
                            "pubmed_id": "23125364"
                        }
                    ],
                    "growth_medium": "Middlebrook 7H9 with the ADC supplement (Difco) and 0.05% Tween80, at 37 degree celcius.",
                    "ordered_locus_name": "MSMEG_1867",
                    "species": "Mycolicibacterium smegmatis MC2 155",
                    "ncbi_taxonomy_id": 246196
                },
                {
                    "halflife": 489.00036341439363,
                    "variation_coefficient": 18.2337631285916,
                    "species": "Escherichia coli str. K-12 substr. MG1655",
                    "ncbi_taxonomy_id": 511145,
                    "unit": "s",
                    "reference": [
                        {
                            "doi": "10.1093/nar/gkt1150"
                        }
                    ],
                    "growth_medium": "M9 minimal medium supplemented with glucose",
                    "ordered_locus_name": "b0060",
                    "doubling_time": {
                        "value": 6.9,
                        "unit": "h"
                    }
                },
                {
                    "transcript_size": 1938,
                    "cds_size": 918,
                    "intron_size": 6906,
                    "genomic_size": 8844,
                    "intron_count": 8,
                    "halflife": 21797.96271,
                    "r_sqaured": 0.988426426,
                    "standard_error": 0.019283288,
                    "unit": "s",
                    "reference": [
                        {
                            "doi": "10.1101/gr.131037.111",
                            "pubmed_id": "22406755"
                        }
                    ],
                    "accession_id": [
                        "AK088066",
                        "AK133695",
                        "BC003426",
                        "NM_145371"
                    ],
                    "ncbi_taxonomy_id": 10090,
                    "species": "Mus musculus"
                },
                {
                    "halflife": 113.13485976,
                    "expression_reads_per_kb_per_mb": 19.37889675,
                    "quantification_method": "Illumina GA-II",
                    "transcriptional_start_sites": 160034,
                    "transcriptional_end_sites": 162771,
                    "unit": "s",
                    "operon": [
                        "TU_160034-162771_F"
                    ],
                    "reference": [
                        {
                            "doi": "10.1186/gb-2012-13-4-r30",
                            "pubmed_id": "22537947"
                        }
                    ],
                    "growth_medium": "Luria-Bertani (LB) broth (500 ml) at 30 degree celcius, 250 rpm.",
                    "ordered_locus_name": "BCE_0159",
                    "gene_start": 160105,
                    "gene_end": 162642,
                    "strand": "F",
                    "cog": "R",
                    "species": "Bacillus cereus ATCC 10987",
                    "ncbi_taxonomy_id": 222523
                },
                {
                    "halflife": 4041.87006,
                    "std": 523.16592,
                    "std_over_avg": 0.1294366004,
                    "unit": "s",
                    "reference": [
                        {
                            "doi": "10.1186/s12864-016-3219-8"
                        }
                    ],
                    "growth_medium": "TMA",
                    "ordered_locus_name": "MA0001",
                    "ar_cog": "arCOG00468",
                    "cog_class": "L",
                    "cog": "COG1474",
                    "species": "Methanosarcina acetivorans",
                    "ncbi_taxonomy_id": 188937
                },
                {
                    "halflife": 1436.4,
                    "species": "Lactococcus lactis subsp. lactis Il1403",
                    "ncbi_taxonomy_id": 272623,
                    "unit": "s",
                    "reference": [
                        {
                            "doi": "10.1371/journal.pone.0059059"
                        }
                    ],
                    "doubling_time": {
                        "value": 6.301338005090412,
                        "unit": "h"
                    }
                }
            ],
            "ko_number": "K10725",
            "protein_names": [
                "ORC1-type DNA replication protein 1",
                "ORC1-type DNA replication protein 1"
            ]
        }
        result = self.src.build_rna_observation(obj)        
        self.assertEqual(len(result), 23)
        self.assertEqual(result[2]["environment"]["replicate"], "a1")