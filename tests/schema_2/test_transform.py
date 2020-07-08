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

    def test_parse_docs(self):
        self.src.process_docs("uniprot")
    
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