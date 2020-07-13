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

    # @unittest.skip("for now")
    def test_parse_docs(self):
        self.src.process_docs("rna_modification",
                              db="datanator")
    
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

    @unittest.skip("passed")
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
        self.assertEqual(result[3]["environment"]["replicate"], "a3")

    @unittest.skip("passed")
    def test_build_rna_modification_entity(self):
        null = None
        obj = {
            "amino_acid": "Ala",
            "aa_code": "A",
            "aa_name": "Alanine",
            "kegg_orthology_id": "K14218",
            "kegg_gene_name": "tRNA-Ala",
            "definition": "tRNA Ala",
            "kegg_pathway_id": "ko00970 ",
            "kegg_pathway_name": "Aminoacyl-tRNA biosynthesis",
            "modifications": [
                {
                    "anticodon": "VGC",
                    "organism": "Escherichia coli",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCUAUAGCUCAGCDGGGAGAGCGCCUGCUUVGCACGCAGGAG7UCUGCGGTPCGAUCCCGCAUAGCUCCACCA",
                    "sequence_bpforms": "GGGGCUAUAGCUCAGC{8U}GGGAGAGCGCCUGCUU{502U}GCACGCAGGAG{7G}UCUGCGG{5U}{9U}CGAUCCCGCAUAGCUCCACCA",
                    "sequence_iupac": "GGGGCUAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCAUAGCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 5,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 4,
                    "formula": "C726H832N289O538P76",
                    "molecular_weight": 24568.13291,
                    "charge": -77,
                    "canonical_formula": "C722H822N289O535P76",
                    "canonical_molecular_weight": 24462.01191,
                    "canonical_charge": -77,
                    "extra_formula": "C4H10O3",
                    "extra_molecular_weight": 106.12100000000001,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:09.439Z"
                    },
                    "ncbi_taxonomy_id": 562
                },
                {
                    "anticodon": "GGC",
                    "organism": "Escherichia coli",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCUAUAGCUCAGCDGGGAGAGCGCUUGCAUGGCAUGCAAGAG7UCAGCGGTPCGAUCCCGCUUAGCUCCACCA",
                    "sequence_bpforms": "GGGGCUAUAGCUCAGC{8U}GGGAGAGCGCUUGCAUGGCAUGCAAGAG{7G}UCAGCGG{5U}{9U}CGAUCCCGCUUAGCUCCACCA",
                    "sequence_iupac": "GGGGCUAUAGCUCAGCUGGGAGAGCGCUUGCAUGGCAUGCAAGAGGUCAGCGGUUCGAUCCCGCUUAGCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 4,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 3,
                    "formula": "C726H831N293O533P76",
                    "molecular_weight": 24543.157909999998,
                    "charge": -76,
                    "canonical_formula": "C724H822N293O533P76",
                    "canonical_molecular_weight": 24510.06391,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9",
                    "extra_molecular_weight": 33.094,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:09.559Z"
                    },
                    "ncbi_taxonomy_id": 562
                },
                {
                    "anticodon": "AGC",
                    "organism": "Saccharomyces cerevisiae",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGCGUGUKGCGUAGDCGGDAGCGCRCUCCCUUIGCOPGGGAGAGGDCUCCGGTPCGAUUCCGGACUCGUCCACCA",
                    "sequence_bpforms": "GGGCGUGU{1G}GCGUAG{8U}CGG{8U}AGCGC{22G}CUCCCUU{9A}GC{19A}{9U}GGGAGAGG{8U}CUCCGG{5U}{9U}CGAUUCCGGACUCGUCCACCA",
                    "sequence_iupac": "GGGCGUGUGGCGUAGUCGGUAGCGCGCUCCCUUAGCAUGGGAGAGGUCUCCGGUUCGAUUCCGGACUCGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 10,
                    "number_of_modified_a": 2,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 2,
                    "number_of_modified_u": 6,
                    "formula": "C726H834N283O542P76",
                    "molecular_weight": 24550.102909999998,
                    "charge": -77,
                    "canonical_formula": "C721H820N285O540P76",
                    "canonical_molecular_weight": 24471.95191,
                    "canonical_charge": -77,
                    "extra_formula": "C5H14N-2O2",
                    "extra_molecular_weight": 78.15100000000001,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:10.684Z"
                    },
                    "ncbi_taxonomy_id": 4932
                },
                {
                    "anticodon": "UGC",
                    "organism": "Bacillus subtilis",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGAGCCUUAGCUCAGCDGGGAGAGCGCCUGCUU5GC=CGCAGGAG7UCAGCGGTPCGAUCCCGCUAGGCUCCACCA",
                    "sequence_bpforms": "GGAGCCUUAGCUCAGC{8U}GGGAGAGCGCCUGCUU{501U}GC{6A}CGCAGGAG{7G}UCAGCGG{5U}{9U}CGAUCCCGCUAGGCUCCACCA",
                    "sequence_iupac": "GGAGCCUUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCAGCGGUUCGAUCCCGCUAGGCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 4,
                    "formula": "C726H836N290O535P76",
                    "molecular_weight": 24538.174909999998,
                    "charge": -76,
                    "canonical_formula": "C722H823N290O534P76",
                    "canonical_molecular_weight": 24461.02791,
                    "canonical_charge": -77,
                    "extra_formula": "C4H13O",
                    "extra_molecular_weight": 77.14699999999999,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:11.889Z"
                    },
                    "ncbi_taxonomy_id": 1423
                },
                {
                    "anticodon": "UGC",
                    "organism": "Mycoplasma capricolum",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCCU4AGCUCAGCDGGGAGAGCACCUGCCUUGC=CGCAGGGG7UCGACGGUPCGAUCCCGUUAGGGUCCACCA",
                    "sequence_bpforms": "GGGCCCU{74U}AGCUCAGC{8U}GGGAGAGCACCUGCCUUGC{6A}CGCAGGGG{7G}UCGACGGU{9U}CGAUCCCGUUAGGGUCCACCA",
                    "sequence_iupac": "GGGCCCUUAGCUCAGCUGGGAGAGCACCUGCCUUGCACGCAGGGGGUCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 5,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 3,
                    "formula": "C724H832N290O534P76S",
                    "molecular_weight": 24526.18191,
                    "charge": -76,
                    "canonical_formula": "C722H823N290O535P76",
                    "canonical_molecular_weight": 24477.02691,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9O-1S",
                    "extra_molecular_weight": 49.155,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:11.918Z"
                    },
                    "ncbi_taxonomy_id": 2095
                },
                {
                    "anticodon": "GGC",
                    "organism": "Bacillus subtilis",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCCAUAGCUCAGCDGGGAGAGCGCUACGCUGGCAGCGUAGAG7UCAGGGGTPCGAGCCCCCUUGGCUCCACCA",
                    "sequence_bpforms": "GGGGCCAUAGCUCAGC{8U}GGGAGAGCGCUACGCUGGCAGCGUAGAG{7G}UCAGGGG{5U}{9U}CGAGCCCCCUUGGCUCCACCA",
                    "sequence_iupac": "GGGGCCAUAGCUCAGCUGGGAGAGCGCUACGCUGGCAGCGUAGAGGUCAGGGGUUCGAGCCCCCUUGGCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 4,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 3,
                    "formula": "C727H834N298O532P76",
                    "molecular_weight": 24612.228909999998,
                    "charge": -76,
                    "canonical_formula": "C725H825N298O532P76",
                    "canonical_molecular_weight": 24579.13491,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9",
                    "extra_molecular_weight": 33.094,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:11.947Z"
                    },
                    "ncbi_taxonomy_id": 1423
                },
                {
                    "anticodon": "CGC",
                    "organism": "Halobacterium salinarum",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCUCGUAGAUCAGCGGUAGAUCRCUUCCUUCGCAAGGAAGAGGCC?UGGG]PBOAAUCCCAGCGAGUCCACCA",
                    "sequence_bpforms": "GGGCUCGUAGAUCAGCGGUAGAUC{22G}CUUCCUUCGCAAGGAAGAGGCC{5C}UGGG{19U}{9U}{0C}{19A}AAUCCCAGCGAGUCCACCA",
                    "sequence_iupac": "GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGAGGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA",
                    "length": 75,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 2,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C720H823N288O524P75",
                    "molecular_weight": 24218.02815,
                    "charge": -76,
                    "canonical_formula": "C714H812N289O523P75",
                    "canonical_molecular_weight": 24132.88215,
                    "canonical_charge": -76,
                    "extra_formula": "C6H11N-1O",
                    "extra_molecular_weight": 85.146,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.398Z"
                    },
                    "ncbi_taxonomy_id": 2242
                },
                {
                    "anticodon": "CGC",
                    "organism": "Haloferax volcanii",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCUCGUAGAUCAGUGGCAGAUCRCUUCCUUCGCAAGGAAGAGGC??GGGG]PBOAAUCCCCGCGAGUCCACCA",
                    "sequence_bpforms": "GGGCUCGUAGAUCAGUGGCAGAUC{22G}CUUCCUUCGCAAGGAAGAGGC{5C}{5C}GGGG{19U}{9U}{0C}{19A}AAUCCCCGCGAGUCCACCA",
                    "sequence_iupac": "GGGCUCGUAGAUCAGUGGCAGAUCGCUUCCUUCGCAAGGAAGAGGCCCGGGGUUCAAAUCCCCGCGAGUCCACCA",
                    "length": 75,
                    "number_of_modifications": 7,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 3,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C721H826N289O524P75",
                    "molecular_weight": 24247.07015,
                    "charge": -76,
                    "canonical_formula": "C714H813N290O523P75",
                    "canonical_molecular_weight": 24147.89715,
                    "canonical_charge": -76,
                    "extra_formula": "C7H13N-1O",
                    "extra_molecular_weight": 99.17299999999999,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.428Z"
                    },
                    "ncbi_taxonomy_id": 2246
                },
                {
                    "anticodon": "GGC",
                    "organism": "Haloferax volcanii",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCUCGUAGAUCAGGGGUAGAUCACUCCCUUGGCAUGGGAGAGGC??CGGG]PBOAAUCCCGGCGAGUCCACCA",
                    "sequence_bpforms": "GGGCUCGUAGAUCAGGGGUAGAUCACUCCCUUGGCAUGGGAGAGGC{5C}{5C}CGGG{19U}{9U}{0C}{19A}AAUCCCGGCGAGUCCACCA",
                    "sequence_iupac": "GGGCUCGUAGAUCAGGGGUAGAUCACUCCCUUGGCAUGGGAGAGGCCCCGGGUUCAAAUCCCGGCGAGUCCACCA",
                    "length": 75,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 3,
                    "number_of_modified_g": 0,
                    "number_of_modified_u": 2,
                    "formula": "C720H822N291O525P75",
                    "molecular_weight": 24275.04015,
                    "charge": -76,
                    "canonical_formula": "C715H813N292O524P75",
                    "canonical_molecular_weight": 24203.92115,
                    "canonical_charge": -76,
                    "extra_formula": "C5H9N-1O",
                    "extra_molecular_weight": 71.119,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.457Z"
                    },
                    "ncbi_taxonomy_id": 2246
                },
                {
                    "anticodon": "UGC",
                    "organism": "Haloferax volcanii",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCCAUAGCUCAGUGGUAGAGULCCUCCUUUGCAAGGAGGAUGC??AGGG]PBGAAUCCCUGUGGGUCCACCA",
                    "sequence_bpforms": "GGGCCCAUAGCUCAGUGGUAGAGU{2G}CCUCCUUUGCAAGGAGGAUGC{5C}{5C}AGGG{19U}{9U}{0C}GAAUCCCUGUGGGUCCACCA",
                    "sequence_iupac": "GGGCCCAUAGCUCAGUGGUAGAGUGCCUCCUUUGCAAGGAGGAUGCCCAGGGUUCGAAUCCCUGUGGGUCCACCA",
                    "length": 75,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 3,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C718H820N285O528P75",
                    "molecular_weight": 24212.95715,
                    "charge": -76,
                    "canonical_formula": "C713H810N285O528P75",
                    "canonical_molecular_weight": 24142.82215,
                    "canonical_charge": -76,
                    "extra_formula": "C5H10",
                    "extra_molecular_weight": 70.135,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.483Z"
                    },
                    "ncbi_taxonomy_id": 2246
                },
                {
                    "anticodon": "UGC",
                    "organism": "Mycoplasma mycoides",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCCUUAGCUCAGCDGGGAGAGCACCUGCCUUGC=CGCAGGGG7UCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "sequence_bpforms": "GGGCCCUUAGCUCAGC{8U}GGGAGAGCACCUGCCUUGC{6A}CGCAGGGG{7G}UCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "sequence_iupac": "GGGCCCUUAGCUCAGCUGGGAGAGCACCUGCCUUGCACGCAGGGGGUCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 3,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 1,
                    "formula": "C724H832N290O535P76",
                    "molecular_weight": 24510.120909999998,
                    "charge": -76,
                    "canonical_formula": "C722H823N290O535P76",
                    "canonical_molecular_weight": 24477.02691,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9",
                    "extra_molecular_weight": 33.094,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.512Z"
                    },
                    "ncbi_taxonomy_id": 2102
                },
                {
                    "anticodon": "UGC",
                    "organism": "Mycoplasma mycoides",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCCUUAGCUCAGCDGGGAGAGCACCUGCCUUGC=CGCAGGGG7UCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "sequence_bpforms": "GGGCCCUUAGCUCAGC{8U}GGGAGAGCACCUGCCUUGC{6A}CGCAGGGG{7G}UCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "sequence_iupac": "GGGCCCUUAGCUCAGCUGGGAGAGCACCUGCCUUGCACGCAGGGGGUCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 3,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 1,
                    "formula": "C724H832N290O535P76",
                    "molecular_weight": 24510.120909999998,
                    "charge": -76,
                    "canonical_formula": "C722H823N290O535P76",
                    "canonical_molecular_weight": 24477.02691,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9",
                    "extra_molecular_weight": 33.094,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.158Z"
                    },
                    "ncbi_taxonomy_id": 2102
                },
                {
                    "anticodon": "IGC",
                    "organism": "Pichia jadinii",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGCGUGUKGCGUAGDDGGDAGCGCRPUCGCUUIGCOPGCGAAAGGDCUCCGGTPCG\"CUCCGGACUCGUCCACCA",
                    "sequence_bpforms": "GGGCGUGU{1G}GCGUAG{8U}{8U}GG{8U}AGCGC{22G}{9U}UCGCUU{9A}GC{19A}{9U}GCGAAAGG{8U}CUCCGG{5U}{9U}CG{1A}CUCCGGACUCGUCCACCA",
                    "sequence_iupac": "GGGCGUGUGGCGUAGUUGGUAGCGCGUUCGCUUAGCAUGCGAAAGGUCUCCGGUUCGACUCCGGACUCGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 13,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 2,
                    "number_of_modified_u": 8,
                    "formula": "C727H837N282O542P76",
                    "molecular_weight": 24551.13091,
                    "charge": -77,
                    "canonical_formula": "C721H819N284O540P76",
                    "canonical_molecular_weight": 24456.93691,
                    "canonical_charge": -77,
                    "extra_formula": "C6H18N-2O2",
                    "extra_molecular_weight": 94.194,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.281Z"
                    },
                    "ncbi_taxonomy_id": null
                },
                {
                    "anticodon": "IGC",
                    "organism": "Bombyx mori",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGGGCGUALCUCAGADGGUAGAGCRCUCGCJUIGCOP#PGAGAG7UA?CGGGAPCG\"UACCCGGCGCCUCCACCA",
                    "sequence_bpforms": "GGGGGCGUA{2G}CUCAGA{8U}GGUAGAGC{22G}CUCGC{0U}U{9A}GC{19A}{9U}{0G}{9U}GAGAG{7G}UA{5C}CGGGA{9U}CG{1A}UACCCGGCGCCUCCACCA",
                    "sequence_iupac": "GGGGGCGUAGCUCAGAUGGUAGAGCGCUCGCUUAGCAUGUGAGAGGUACCGGGAUCGAUACCCGGCGCCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 13,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 4,
                    "number_of_modified_u": 5,
                    "formula": "C735H845N297O533P76",
                    "molecular_weight": 24721.39691,
                    "charge": -76,
                    "canonical_formula": "C726H824N299O531P76",
                    "canonical_molecular_weight": 24588.14591,
                    "canonical_charge": -77,
                    "extra_formula": "C9H21N-2O2",
                    "extra_molecular_weight": 133.251,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.310Z"
                    },
                    "ncbi_taxonomy_id": 7091
                },
                {
                    "anticodon": "IGC",
                    "organism": "Bombyx mori",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGGGCGUALCUCAGADGGUAGAGCRCUCGCJUIGCOP#PGAGAG7UA?CGGGAPCG\"UACCCGGCGCCUCCACCA",
                    "sequence_bpforms": "GGGGGCGUA{2G}CUCAGA{8U}GGUAGAGC{22G}CUCGC{0U}U{9A}GC{19A}{9U}{0G}{9U}GAGAG{7G}UA{5C}CGGGA{9U}CG{1A}UACCCGGCGCCUCCACCA",
                    "sequence_iupac": "GGGGGCGUAGCUCAGAUGGUAGAGCGCUCGCUUAGCAUGUGAGAGGUACCGGGAUCGAUACCCGGCGCCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 13,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 4,
                    "number_of_modified_u": 5,
                    "formula": "C735H845N297O533P76",
                    "molecular_weight": 24721.39691,
                    "charge": -76,
                    "canonical_formula": "C726H824N299O531P76",
                    "canonical_molecular_weight": 24588.14591,
                    "canonical_charge": -77,
                    "extra_formula": "C9H21N-2O2",
                    "extra_molecular_weight": 133.251,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.345Z"
                    },
                    "ncbi_taxonomy_id": 7091
                },
                {
                    "anticodon": "IGC",
                    "organism": "Homo sapiens",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGGGAUUALCUCAAADGGDAGAGCRCUCGCJUIGCOP#CGAGAG7UAGCGGGAPCG\"UGCCCGCAUCCUCCACCA",
                    "sequence_bpforms": "GGGGGAUUA{2G}CUCAAA{8U}GG{8U}AGAGC{22G}CUCGC{0U}U{9A}GC{19A}{9U}{0G}CGAGAG{7G}UAGCGGGA{9U}CG{1A}UGCCCGCAUCCUCCACCA",
                    "sequence_iupac": "GGGGGAUUAGCUCAAAUGGUAGAGCGCUCGCUUAGCAUGCGAGAGGUAGCGGGAUCGAUGCCCGCAUCCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 12,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 4,
                    "number_of_modified_u": 5,
                    "formula": "C734H844N296O532P76",
                    "molecular_weight": 24678.371909999998,
                    "charge": -76,
                    "canonical_formula": "C726H823N298O530P76",
                    "canonical_molecular_weight": 24557.13191,
                    "canonical_charge": -77,
                    "extra_formula": "C8H21N-2O2",
                    "extra_molecular_weight": 121.24,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.375Z"
                    },
                    "ncbi_taxonomy_id": 9606
                },
                {
                    "anticodon": "IGC",
                    "organism": "Homo sapiens",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGGAAUUALCUCAAADGGDAGAGCRCUCGCJUIGCOP#CGAGAG7UAGCGGGAPCG\"UGCCCGCAUUCUCCACCA",
                    "sequence_bpforms": "GGGGAAUUA{2G}CUCAAA{8U}GG{8U}AGAGC{22G}CUCGC{0U}U{9A}GC{19A}{9U}{0G}CGAGAG{7G}UAGCGGGA{9U}CG{1A}UGCCCGCAUUCUCCACCA",
                    "sequence_iupac": "GGGGAAUUAGCUCAAAUGGUAGAGCGCUCGCUUAGCAUGCGAGAGGUAGCGGGAUCGAUGCCCGCAUUCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 12,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 4,
                    "number_of_modified_u": 5,
                    "formula": "C734H843N295O532P76",
                    "molecular_weight": 24663.35691,
                    "charge": -76,
                    "canonical_formula": "C726H822N297O530P76",
                    "canonical_molecular_weight": 24542.11691,
                    "canonical_charge": -77,
                    "extra_formula": "C8H21N-2O2",
                    "extra_molecular_weight": 121.24,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.403Z"
                    },
                    "ncbi_taxonomy_id": 9606
                },
                {
                    "anticodon": "UGC",
                    "organism": "Neurospora crassa",
                    "organellum": "mitochondrial",
                    "sequence_modomics": "GGGGGUAUAGUAUAADUGGDAGUACAGCAAUCUUGCUCANUGCUUGU?AAGGTPCAAAUCCUUGUAUCUCCACCA",
                    "sequence_bpforms": "GGGGGUAUAGUAUAA{8U}UGG{8U}AGUACAGCAAUCUUGCUCA[id: \"xU\"]UGCUUGU{5C}AAGG{5U}{9U}CAAAUCCUUGUAUCUCCACCA",
                    "sequence_iupac": "GGGGGUAUAGUAUAAUUGGUAGUACAGCAAUCUUGCUCANUGCUUGUCAAGGUUCAAAUCCUUGUAUCUCCACCA",
                    "length": 75,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 0,
                    "number_of_modified_u": 4,
                    "formula": null,
                    "molecular_weight": null,
                    "charge": null,
                    "canonical_formula": null,
                    "canonical_molecular_weight": null,
                    "canonical_charge": null,
                    "extra_formula": null,
                    "extra_molecular_weight": null,
                    "extra_charge": null,
                    "bpforms_errors": "MODOMICS sequence uses monomeric forms xU",
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:56.216Z"
                    },
                    "ncbi_taxonomy_id": 5141
                },
                {
                    "anticodon": "UGC",
                    "organism": "Bos taurus",
                    "organellum": "mitochondrial",
                    "sequence_modomics": "GAGGAUUU\"LCUUAAUUAAAGULGPUGAUUUGCAUPCAAUUGAUGUAAGGUGPAGUCUUGCAAUCCUUACCA",
                    "sequence_bpforms": "GAGGAUUU{1A}{2G}CUUAAUUAAAGU{2G}G{9U}UGAUUUGCAU{9U}CAAUUGAUGUAAGGUG{9U}AGUCUUGCAAUCCUUACCA",
                    "sequence_iupac": "GAGGAUUUAGCUUAAUUAAAGUGGUUGAUUUGCAUUCAAUUGAUGUAAGGUGUAGUCUUGCAAUCCUUACCA",
                    "length": 72,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 2,
                    "number_of_modified_u": 3,
                    "formula": "C687H772N261O512P72",
                    "molecular_weight": 23107.15886,
                    "charge": -73,
                    "canonical_formula": "C684H766N261O512P72",
                    "canonical_molecular_weight": 23065.077859999998,
                    "canonical_charge": -73,
                    "extra_formula": "C3H6",
                    "extra_molecular_weight": 42.081,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:00.410Z"
                    },
                    "ncbi_taxonomy_id": 9913
                },
                {
                    "anticodon": "UGC",
                    "organism": "Lactococcus lactis",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCCU4AGCUCAGCUGGGAGAGCGCCUGCUU5GC6CGCAGGAG7UCAGCGGUUCGAUCCCGCUAGGCUCCA",
                    "sequence_bpforms": "GGGGCCU{74U}AGCUCAGCUGGGAGAGCGCCUGCUU{501U}GC{62A}CGCAGGAG{7G}UCAGCGGUUCGAUCCCGCUAGGCUCCA",
                    "sequence_iupac": "GGGGCCUUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCAGCGGUUCGAUCCCGCUAGGCUCCA",
                    "length": 73,
                    "number_of_modifications": 4,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C701H805N281O518P73S",
                    "molecular_weight": 23747.74463,
                    "charge": -73,
                    "canonical_formula": "C694H790N279O515P73",
                    "canonical_molecular_weight": 23540.47663,
                    "canonical_charge": -74,
                    "extra_formula": "C7H15N2O3S",
                    "extra_molecular_weight": 207.268,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:00.788Z"
                    },
                    "ncbi_taxonomy_id": 1358
                },
                {
                    "anticodon": "GGC",
                    "organism": "Streptomyces griseus",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCUAUAGCUBAGUDGGDAGAGCGCCUGCAUGGCAUGCAGGAG7UCAGGAGUUCA\"UUCUCCUUAGCUCCACAA",
                    "sequence_bpforms": "GGGGCUAUAGCU{0C}AGU{8U}GG{8U}AGAGCGCCUGCAUGGCAUGCAGGAG{7G}UCAGGAGUUCA{1A}UUCUCCUUAGCUCCACAA",
                    "sequence_iupac": "GGGGCUAUAGCUCAGUUGGUAGAGCGCCUGCAUGGCAUGCAGGAGGUCAGGAGUUCAAUUCUCCUUAGCUCCACAA",
                    "length": 76,
                    "number_of_modifications": 5,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C727H832N290O534P76",
                    "molecular_weight": 24530.15491,
                    "charge": -76,
                    "canonical_formula": "C724H819N290O534P76",
                    "canonical_molecular_weight": 24481.01791,
                    "canonical_charge": -77,
                    "extra_formula": "C3H13",
                    "extra_molecular_weight": 49.137,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:02.506Z"
                    },
                    "ncbi_taxonomy_id": 1911
                },
                {
                    "anticodon": "GGC",
                    "organism": "Streptomyces griseus",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCUAUAGCUBAGUDGGDAGAGCGCCUGCAUGGCAUGCAGGAG7UCAGGAGUUCA\"UUCUCCUUAGCUCCA",
                    "sequence_bpforms": "GGGGCUAUAGCU{0C}AGU{8U}GG{8U}AGAGCGCCUGCAUGGCAUGCAGGAG{7G}UCAGGAGUUCA{1A}UUCUCCUUAGCUCCA",
                    "sequence_iupac": "GGGGCUAUAGCUCAGUUGGUAGAGCGCCUGCAUGGCAUGCAGGAGGUCAGGAGUUCAAUUCUCCUUAGCUCCA",
                    "length": 73,
                    "number_of_modifications": 5,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C698H799N277O515P73",
                    "molecular_weight": 23569.57863,
                    "charge": -73,
                    "canonical_formula": "C695H786N277O515P73",
                    "canonical_molecular_weight": 23520.44163,
                    "canonical_charge": -74,
                    "extra_formula": "C3H13",
                    "extra_molecular_weight": 49.137,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:02.535Z"
                    },
                    "ncbi_taxonomy_id": 1911
                },
                {
                    "anticodon": "CGC",
                    "organism": "Streptomyces griseus",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCUGUGGCGCAGUCUGGDAGCGCACCUCGUUCGCAUCGAGGGGGUCUGGGGUUCA\"AUCCCCACAGGUCCA",
                    "sequence_bpforms": "GGGCCUGUGGCGCAGUCUGG{8U}AGCGCACCUCGUUCGCAUCGAGGGGGUCUGGGGUUCA{1A}AUCCCCACAGGUCCA",
                    "sequence_iupac": "GGGCCUGUGGCGCAGUCUGGUAGCGCACCUCGUUCGCAUCGAGGGGGUCUGGGGUUCAAAUCCCCACAGGUCCA",
                    "length": 74,
                    "number_of_modifications": 2,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 0,
                    "number_of_modified_u": 1,
                    "formula": "C704H804N281O523P74",
                    "molecular_weight": 23861.67839,
                    "charge": -75,
                    "canonical_formula": "C703H800N281O523P74",
                    "canonical_molecular_weight": 23845.63539,
                    "canonical_charge": -75,
                    "extra_formula": "CH4",
                    "extra_molecular_weight": 16.043,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:02.564Z"
                    },
                    "ncbi_taxonomy_id": 1911
                },
                {
                    "anticodon": "UGC",
                    "organism": "Streptomyces griseus",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCCUUAGCUBAGUDGGDAGAGCGCUGCCUUUGCAAGGCAGAUGUCAGGAGUUCGAAUCUCCUAGGCUCCA",
                    "sequence_bpforms": "GGGGCCUUAGCU{0C}AGU{8U}GG{8U}AGAGCGCUGCCUUUGCAAGGCAGAUGUCAGGAGUUCGAAUCUCCUAGGCUCCA",
                    "sequence_iupac": "GGGGCCUUAGCUCAGUUGGUAGAGCGCUGCCUUUGCAAGGCAGAUGUCAGGAGUUCGAAUCUCCUAGGCUCCA",
                    "length": 73,
                    "number_of_modifications": 3,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 0,
                    "number_of_modified_u": 2,
                    "formula": "C695H792N275O516P73",
                    "molecular_weight": 23514.47463,
                    "charge": -74,
                    "canonical_formula": "C694H786N275O516P73",
                    "canonical_molecular_weight": 23496.41563,
                    "canonical_charge": -74,
                    "extra_formula": "CH6",
                    "extra_molecular_weight": 18.059,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:02.594Z"
                    },
                    "ncbi_taxonomy_id": 1911
                }
            ]
        }
        result = self.src.build_rna_modification_entity(obj)
        print(result)

    @unittest.skip("passed")
    def test_build_rna_modification_observation(self):
        null = None
        obj = {
            "amino_acid": "Ala",
            "aa_code": "A",
            "aa_name": "Alanine",
            "kegg_orthology_id": "K14218",
            "kegg_gene_name": "tRNA-Ala",
            "definition": "tRNA Ala",
            "kegg_pathway_id": "ko00970 ",
            "kegg_pathway_name": "Aminoacyl-tRNA biosynthesis",
            "modifications": [
                {
                    "anticodon": "VGC",
                    "organism": "Escherichia coli",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCUAUAGCUCAGCDGGGAGAGCGCCUGCUUVGCACGCAGGAG7UCUGCGGTPCGAUCCCGCAUAGCUCCACCA",
                    "sequence_bpforms": "GGGGCUAUAGCUCAGC{8U}GGGAGAGCGCCUGCUU{502U}GCACGCAGGAG{7G}UCUGCGG{5U}{9U}CGAUCCCGCAUAGCUCCACCA",
                    "sequence_iupac": "GGGGCUAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCAUAGCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 5,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 4,
                    "formula": "C726H832N289O538P76",
                    "molecular_weight": 24568.13291,
                    "charge": -77,
                    "canonical_formula": "C722H822N289O535P76",
                    "canonical_molecular_weight": 24462.01191,
                    "canonical_charge": -77,
                    "extra_formula": "C4H10O3",
                    "extra_molecular_weight": 106.12100000000001,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:09.439Z"
                    },
                    "ncbi_taxonomy_id": 562
                },
                {
                    "anticodon": "GGC",
                    "organism": "Escherichia coli",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCUAUAGCUCAGCDGGGAGAGCGCUUGCAUGGCAUGCAAGAG7UCAGCGGTPCGAUCCCGCUUAGCUCCACCA",
                    "sequence_bpforms": "GGGGCUAUAGCUCAGC{8U}GGGAGAGCGCUUGCAUGGCAUGCAAGAG{7G}UCAGCGG{5U}{9U}CGAUCCCGCUUAGCUCCACCA",
                    "sequence_iupac": "GGGGCUAUAGCUCAGCUGGGAGAGCGCUUGCAUGGCAUGCAAGAGGUCAGCGGUUCGAUCCCGCUUAGCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 4,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 3,
                    "formula": "C726H831N293O533P76",
                    "molecular_weight": 24543.157909999998,
                    "charge": -76,
                    "canonical_formula": "C724H822N293O533P76",
                    "canonical_molecular_weight": 24510.06391,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9",
                    "extra_molecular_weight": 33.094,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:09.559Z"
                    },
                    "ncbi_taxonomy_id": 562
                },
                {
                    "anticodon": "AGC",
                    "organism": "Saccharomyces cerevisiae",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGCGUGUKGCGUAGDCGGDAGCGCRCUCCCUUIGCOPGGGAGAGGDCUCCGGTPCGAUUCCGGACUCGUCCACCA",
                    "sequence_bpforms": "GGGCGUGU{1G}GCGUAG{8U}CGG{8U}AGCGC{22G}CUCCCUU{9A}GC{19A}{9U}GGGAGAGG{8U}CUCCGG{5U}{9U}CGAUUCCGGACUCGUCCACCA",
                    "sequence_iupac": "GGGCGUGUGGCGUAGUCGGUAGCGCGCUCCCUUAGCAUGGGAGAGGUCUCCGGUUCGAUUCCGGACUCGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 10,
                    "number_of_modified_a": 2,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 2,
                    "number_of_modified_u": 6,
                    "formula": "C726H834N283O542P76",
                    "molecular_weight": 24550.102909999998,
                    "charge": -77,
                    "canonical_formula": "C721H820N285O540P76",
                    "canonical_molecular_weight": 24471.95191,
                    "canonical_charge": -77,
                    "extra_formula": "C5H14N-2O2",
                    "extra_molecular_weight": 78.15100000000001,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:10.684Z"
                    },
                    "ncbi_taxonomy_id": 4932
                },
                {
                    "anticodon": "UGC",
                    "organism": "Bacillus subtilis",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGAGCCUUAGCUCAGCDGGGAGAGCGCCUGCUU5GC=CGCAGGAG7UCAGCGGTPCGAUCCCGCUAGGCUCCACCA",
                    "sequence_bpforms": "GGAGCCUUAGCUCAGC{8U}GGGAGAGCGCCUGCUU{501U}GC{6A}CGCAGGAG{7G}UCAGCGG{5U}{9U}CGAUCCCGCUAGGCUCCACCA",
                    "sequence_iupac": "GGAGCCUUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCAGCGGUUCGAUCCCGCUAGGCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 4,
                    "formula": "C726H836N290O535P76",
                    "molecular_weight": 24538.174909999998,
                    "charge": -76,
                    "canonical_formula": "C722H823N290O534P76",
                    "canonical_molecular_weight": 24461.02791,
                    "canonical_charge": -77,
                    "extra_formula": "C4H13O",
                    "extra_molecular_weight": 77.14699999999999,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:11.889Z"
                    },
                    "ncbi_taxonomy_id": 1423
                },
                {
                    "anticodon": "UGC",
                    "organism": "Mycoplasma capricolum",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCCU4AGCUCAGCDGGGAGAGCACCUGCCUUGC=CGCAGGGG7UCGACGGUPCGAUCCCGUUAGGGUCCACCA",
                    "sequence_bpforms": "GGGCCCU{74U}AGCUCAGC{8U}GGGAGAGCACCUGCCUUGC{6A}CGCAGGGG{7G}UCGACGGU{9U}CGAUCCCGUUAGGGUCCACCA",
                    "sequence_iupac": "GGGCCCUUAGCUCAGCUGGGAGAGCACCUGCCUUGCACGCAGGGGGUCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 5,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 3,
                    "formula": "C724H832N290O534P76S",
                    "molecular_weight": 24526.18191,
                    "charge": -76,
                    "canonical_formula": "C722H823N290O535P76",
                    "canonical_molecular_weight": 24477.02691,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9O-1S",
                    "extra_molecular_weight": 49.155,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:11.918Z"
                    },
                    "ncbi_taxonomy_id": 2095
                },
                {
                    "anticodon": "GGC",
                    "organism": "Bacillus subtilis",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCCAUAGCUCAGCDGGGAGAGCGCUACGCUGGCAGCGUAGAG7UCAGGGGTPCGAGCCCCCUUGGCUCCACCA",
                    "sequence_bpforms": "GGGGCCAUAGCUCAGC{8U}GGGAGAGCGCUACGCUGGCAGCGUAGAG{7G}UCAGGGG{5U}{9U}CGAGCCCCCUUGGCUCCACCA",
                    "sequence_iupac": "GGGGCCAUAGCUCAGCUGGGAGAGCGCUACGCUGGCAGCGUAGAGGUCAGGGGUUCGAGCCCCCUUGGCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 4,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 3,
                    "formula": "C727H834N298O532P76",
                    "molecular_weight": 24612.228909999998,
                    "charge": -76,
                    "canonical_formula": "C725H825N298O532P76",
                    "canonical_molecular_weight": 24579.13491,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9",
                    "extra_molecular_weight": 33.094,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:11.947Z"
                    },
                    "ncbi_taxonomy_id": 1423
                },
                {
                    "anticodon": "CGC",
                    "organism": "Halobacterium salinarum",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCUCGUAGAUCAGCGGUAGAUCRCUUCCUUCGCAAGGAAGAGGCC?UGGG]PBOAAUCCCAGCGAGUCCACCA",
                    "sequence_bpforms": "GGGCUCGUAGAUCAGCGGUAGAUC{22G}CUUCCUUCGCAAGGAAGAGGCC{5C}UGGG{19U}{9U}{0C}{19A}AAUCCCAGCGAGUCCACCA",
                    "sequence_iupac": "GGGCUCGUAGAUCAGCGGUAGAUCGCUUCCUUCGCAAGGAAGAGGCCCUGGGUUCAAAUCCCAGCGAGUCCACCA",
                    "length": 75,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 2,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C720H823N288O524P75",
                    "molecular_weight": 24218.02815,
                    "charge": -76,
                    "canonical_formula": "C714H812N289O523P75",
                    "canonical_molecular_weight": 24132.88215,
                    "canonical_charge": -76,
                    "extra_formula": "C6H11N-1O",
                    "extra_molecular_weight": 85.146,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.398Z"
                    },
                    "ncbi_taxonomy_id": 2242
                },
                {
                    "anticodon": "CGC",
                    "organism": "Haloferax volcanii",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCUCGUAGAUCAGUGGCAGAUCRCUUCCUUCGCAAGGAAGAGGC??GGGG]PBOAAUCCCCGCGAGUCCACCA",
                    "sequence_bpforms": "GGGCUCGUAGAUCAGUGGCAGAUC{22G}CUUCCUUCGCAAGGAAGAGGC{5C}{5C}GGGG{19U}{9U}{0C}{19A}AAUCCCCGCGAGUCCACCA",
                    "sequence_iupac": "GGGCUCGUAGAUCAGUGGCAGAUCGCUUCCUUCGCAAGGAAGAGGCCCGGGGUUCAAAUCCCCGCGAGUCCACCA",
                    "length": 75,
                    "number_of_modifications": 7,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 3,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C721H826N289O524P75",
                    "molecular_weight": 24247.07015,
                    "charge": -76,
                    "canonical_formula": "C714H813N290O523P75",
                    "canonical_molecular_weight": 24147.89715,
                    "canonical_charge": -76,
                    "extra_formula": "C7H13N-1O",
                    "extra_molecular_weight": 99.17299999999999,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.428Z"
                    },
                    "ncbi_taxonomy_id": 2246
                },
                {
                    "anticodon": "GGC",
                    "organism": "Haloferax volcanii",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCUCGUAGAUCAGGGGUAGAUCACUCCCUUGGCAUGGGAGAGGC??CGGG]PBOAAUCCCGGCGAGUCCACCA",
                    "sequence_bpforms": "GGGCUCGUAGAUCAGGGGUAGAUCACUCCCUUGGCAUGGGAGAGGC{5C}{5C}CGGG{19U}{9U}{0C}{19A}AAUCCCGGCGAGUCCACCA",
                    "sequence_iupac": "GGGCUCGUAGAUCAGGGGUAGAUCACUCCCUUGGCAUGGGAGAGGCCCCGGGUUCAAAUCCCGGCGAGUCCACCA",
                    "length": 75,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 3,
                    "number_of_modified_g": 0,
                    "number_of_modified_u": 2,
                    "formula": "C720H822N291O525P75",
                    "molecular_weight": 24275.04015,
                    "charge": -76,
                    "canonical_formula": "C715H813N292O524P75",
                    "canonical_molecular_weight": 24203.92115,
                    "canonical_charge": -76,
                    "extra_formula": "C5H9N-1O",
                    "extra_molecular_weight": 71.119,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.457Z"
                    },
                    "ncbi_taxonomy_id": 2246
                },
                {
                    "anticodon": "UGC",
                    "organism": "Haloferax volcanii",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCCAUAGCUCAGUGGUAGAGULCCUCCUUUGCAAGGAGGAUGC??AGGG]PBGAAUCCCUGUGGGUCCACCA",
                    "sequence_bpforms": "GGGCCCAUAGCUCAGUGGUAGAGU{2G}CCUCCUUUGCAAGGAGGAUGC{5C}{5C}AGGG{19U}{9U}{0C}GAAUCCCUGUGGGUCCACCA",
                    "sequence_iupac": "GGGCCCAUAGCUCAGUGGUAGAGUGCCUCCUUUGCAAGGAGGAUGCCCAGGGUUCGAAUCCCUGUGGGUCCACCA",
                    "length": 75,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 3,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C718H820N285O528P75",
                    "molecular_weight": 24212.95715,
                    "charge": -76,
                    "canonical_formula": "C713H810N285O528P75",
                    "canonical_molecular_weight": 24142.82215,
                    "canonical_charge": -76,
                    "extra_formula": "C5H10",
                    "extra_molecular_weight": 70.135,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.483Z"
                    },
                    "ncbi_taxonomy_id": 2246
                },
                {
                    "anticodon": "UGC",
                    "organism": "Mycoplasma mycoides",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCCUUAGCUCAGCDGGGAGAGCACCUGCCUUGC=CGCAGGGG7UCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "sequence_bpforms": "GGGCCCUUAGCUCAGC{8U}GGGAGAGCACCUGCCUUGC{6A}CGCAGGGG{7G}UCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "sequence_iupac": "GGGCCCUUAGCUCAGCUGGGAGAGCACCUGCCUUGCACGCAGGGGGUCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 3,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 1,
                    "formula": "C724H832N290O535P76",
                    "molecular_weight": 24510.120909999998,
                    "charge": -76,
                    "canonical_formula": "C722H823N290O535P76",
                    "canonical_molecular_weight": 24477.02691,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9",
                    "extra_molecular_weight": 33.094,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:31:13.512Z"
                    },
                    "ncbi_taxonomy_id": 2102
                },
                {
                    "anticodon": "UGC",
                    "organism": "Mycoplasma mycoides",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCCUUAGCUCAGCDGGGAGAGCACCUGCCUUGC=CGCAGGGG7UCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "sequence_bpforms": "GGGCCCUUAGCUCAGC{8U}GGGAGAGCACCUGCCUUGC{6A}CGCAGGGG{7G}UCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "sequence_iupac": "GGGCCCUUAGCUCAGCUGGGAGAGCACCUGCCUUGCACGCAGGGGGUCGACGGUUCGAUCCCGUUAGGGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 3,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 1,
                    "formula": "C724H832N290O535P76",
                    "molecular_weight": 24510.120909999998,
                    "charge": -76,
                    "canonical_formula": "C722H823N290O535P76",
                    "canonical_molecular_weight": 24477.02691,
                    "canonical_charge": -77,
                    "extra_formula": "C2H9",
                    "extra_molecular_weight": 33.094,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.158Z"
                    },
                    "ncbi_taxonomy_id": 2102
                },
                {
                    "anticodon": "IGC",
                    "organism": "Pichia jadinii",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGCGUGUKGCGUAGDDGGDAGCGCRPUCGCUUIGCOPGCGAAAGGDCUCCGGTPCG\"CUCCGGACUCGUCCACCA",
                    "sequence_bpforms": "GGGCGUGU{1G}GCGUAG{8U}{8U}GG{8U}AGCGC{22G}{9U}UCGCUU{9A}GC{19A}{9U}GCGAAAGG{8U}CUCCGG{5U}{9U}CG{1A}CUCCGGACUCGUCCACCA",
                    "sequence_iupac": "GGGCGUGUGGCGUAGUUGGUAGCGCGUUCGCUUAGCAUGCGAAAGGUCUCCGGUUCGACUCCGGACUCGUCCACCA",
                    "length": 76,
                    "number_of_modifications": 13,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 2,
                    "number_of_modified_u": 8,
                    "formula": "C727H837N282O542P76",
                    "molecular_weight": 24551.13091,
                    "charge": -77,
                    "canonical_formula": "C721H819N284O540P76",
                    "canonical_molecular_weight": 24456.93691,
                    "canonical_charge": -77,
                    "extra_formula": "C6H18N-2O2",
                    "extra_molecular_weight": 94.194,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.281Z"
                    },
                    "ncbi_taxonomy_id": null
                },
                {
                    "anticodon": "IGC",
                    "organism": "Bombyx mori",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGGGCGUALCUCAGADGGUAGAGCRCUCGCJUIGCOP#PGAGAG7UA?CGGGAPCG\"UACCCGGCGCCUCCACCA",
                    "sequence_bpforms": "GGGGGCGUA{2G}CUCAGA{8U}GGUAGAGC{22G}CUCGC{0U}U{9A}GC{19A}{9U}{0G}{9U}GAGAG{7G}UA{5C}CGGGA{9U}CG{1A}UACCCGGCGCCUCCACCA",
                    "sequence_iupac": "GGGGGCGUAGCUCAGAUGGUAGAGCGCUCGCUUAGCAUGUGAGAGGUACCGGGAUCGAUACCCGGCGCCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 13,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 4,
                    "number_of_modified_u": 5,
                    "formula": "C735H845N297O533P76",
                    "molecular_weight": 24721.39691,
                    "charge": -76,
                    "canonical_formula": "C726H824N299O531P76",
                    "canonical_molecular_weight": 24588.14591,
                    "canonical_charge": -77,
                    "extra_formula": "C9H21N-2O2",
                    "extra_molecular_weight": 133.251,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.310Z"
                    },
                    "ncbi_taxonomy_id": 7091
                },
                {
                    "anticodon": "IGC",
                    "organism": "Bombyx mori",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGGGCGUALCUCAGADGGUAGAGCRCUCGCJUIGCOP#PGAGAG7UA?CGGGAPCG\"UACCCGGCGCCUCCACCA",
                    "sequence_bpforms": "GGGGGCGUA{2G}CUCAGA{8U}GGUAGAGC{22G}CUCGC{0U}U{9A}GC{19A}{9U}{0G}{9U}GAGAG{7G}UA{5C}CGGGA{9U}CG{1A}UACCCGGCGCCUCCACCA",
                    "sequence_iupac": "GGGGGCGUAGCUCAGAUGGUAGAGCGCUCGCUUAGCAUGUGAGAGGUACCGGGAUCGAUACCCGGCGCCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 13,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 4,
                    "number_of_modified_u": 5,
                    "formula": "C735H845N297O533P76",
                    "molecular_weight": 24721.39691,
                    "charge": -76,
                    "canonical_formula": "C726H824N299O531P76",
                    "canonical_molecular_weight": 24588.14591,
                    "canonical_charge": -77,
                    "extra_formula": "C9H21N-2O2",
                    "extra_molecular_weight": 133.251,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.345Z"
                    },
                    "ncbi_taxonomy_id": 7091
                },
                {
                    "anticodon": "IGC",
                    "organism": "Homo sapiens",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGGGAUUALCUCAAADGGDAGAGCRCUCGCJUIGCOP#CGAGAG7UAGCGGGAPCG\"UGCCCGCAUCCUCCACCA",
                    "sequence_bpforms": "GGGGGAUUA{2G}CUCAAA{8U}GG{8U}AGAGC{22G}CUCGC{0U}U{9A}GC{19A}{9U}{0G}CGAGAG{7G}UAGCGGGA{9U}CG{1A}UGCCCGCAUCCUCCACCA",
                    "sequence_iupac": "GGGGGAUUAGCUCAAAUGGUAGAGCGCUCGCUUAGCAUGCGAGAGGUAGCGGGAUCGAUGCCCGCAUCCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 12,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 4,
                    "number_of_modified_u": 5,
                    "formula": "C734H844N296O532P76",
                    "molecular_weight": 24678.371909999998,
                    "charge": -76,
                    "canonical_formula": "C726H823N298O530P76",
                    "canonical_molecular_weight": 24557.13191,
                    "canonical_charge": -77,
                    "extra_formula": "C8H21N-2O2",
                    "extra_molecular_weight": 121.24,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.375Z"
                    },
                    "ncbi_taxonomy_id": 9606
                },
                {
                    "anticodon": "IGC",
                    "organism": "Homo sapiens",
                    "organellum": "cytosolic",
                    "sequence_modomics": "GGGGAAUUALCUCAAADGGDAGAGCRCUCGCJUIGCOP#CGAGAG7UAGCGGGAPCG\"UGCCCGCAUUCUCCACCA",
                    "sequence_bpforms": "GGGGAAUUA{2G}CUCAAA{8U}GG{8U}AGAGC{22G}CUCGC{0U}U{9A}GC{19A}{9U}{0G}CGAGAG{7G}UAGCGGGA{9U}CG{1A}UGCCCGCAUUCUCCACCA",
                    "sequence_iupac": "GGGGAAUUAGCUCAAAUGGUAGAGCGCUCGCUUAGCAUGCGAGAGGUAGCGGGAUCGAUGCCCGCAUUCUCCACCA",
                    "length": 76,
                    "number_of_modifications": 12,
                    "number_of_modified_a": 3,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 4,
                    "number_of_modified_u": 5,
                    "formula": "C734H843N295O532P76",
                    "molecular_weight": 24663.35691,
                    "charge": -76,
                    "canonical_formula": "C726H822N297O530P76",
                    "canonical_molecular_weight": 24542.11691,
                    "canonical_charge": -77,
                    "extra_formula": "C8H21N-2O2",
                    "extra_molecular_weight": 121.24,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:47.403Z"
                    },
                    "ncbi_taxonomy_id": 9606
                },
                {
                    "anticodon": "UGC",
                    "organism": "Neurospora crassa",
                    "organellum": "mitochondrial",
                    "sequence_modomics": "GGGGGUAUAGUAUAADUGGDAGUACAGCAAUCUUGCUCANUGCUUGU?AAGGTPCAAAUCCUUGUAUCUCCACCA",
                    "sequence_bpforms": "GGGGGUAUAGUAUAA{8U}UGG{8U}AGUACAGCAAUCUUGCUCA[id: \"xU\"]UGCUUGU{5C}AAGG{5U}{9U}CAAAUCCUUGUAUCUCCACCA",
                    "sequence_iupac": "GGGGGUAUAGUAUAAUUGGUAGUACAGCAAUCUUGCUCANUGCUUGUCAAGGUUCAAAUCCUUGUAUCUCCACCA",
                    "length": 75,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 0,
                    "number_of_modified_u": 4,
                    "formula": null,
                    "molecular_weight": null,
                    "charge": null,
                    "canonical_formula": null,
                    "canonical_molecular_weight": null,
                    "canonical_charge": null,
                    "extra_formula": null,
                    "extra_molecular_weight": null,
                    "extra_charge": null,
                    "bpforms_errors": "MODOMICS sequence uses monomeric forms xU",
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:33:56.216Z"
                    },
                    "ncbi_taxonomy_id": 5141
                },
                {
                    "anticodon": "UGC",
                    "organism": "Bos taurus",
                    "organellum": "mitochondrial",
                    "sequence_modomics": "GAGGAUUU\"LCUUAAUUAAAGULGPUGAUUUGCAUPCAAUUGAUGUAAGGUGPAGUCUUGCAAUCCUUACCA",
                    "sequence_bpforms": "GAGGAUUU{1A}{2G}CUUAAUUAAAGU{2G}G{9U}UGAUUUGCAU{9U}CAAUUGAUGUAAGGUG{9U}AGUCUUGCAAUCCUUACCA",
                    "sequence_iupac": "GAGGAUUUAGCUUAAUUAAAGUGGUUGAUUUGCAUUCAAUUGAUGUAAGGUGUAGUCUUGCAAUCCUUACCA",
                    "length": 72,
                    "number_of_modifications": 6,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 2,
                    "number_of_modified_u": 3,
                    "formula": "C687H772N261O512P72",
                    "molecular_weight": 23107.15886,
                    "charge": -73,
                    "canonical_formula": "C684H766N261O512P72",
                    "canonical_molecular_weight": 23065.077859999998,
                    "canonical_charge": -73,
                    "extra_formula": "C3H6",
                    "extra_molecular_weight": 42.081,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:00.410Z"
                    },
                    "ncbi_taxonomy_id": 9913
                },
                {
                    "anticodon": "UGC",
                    "organism": "Lactococcus lactis",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCCU4AGCUCAGCUGGGAGAGCGCCUGCUU5GC6CGCAGGAG7UCAGCGGUUCGAUCCCGCUAGGCUCCA",
                    "sequence_bpforms": "GGGGCCU{74U}AGCUCAGCUGGGAGAGCGCCUGCUU{501U}GC{62A}CGCAGGAG{7G}UCAGCGGUUCGAUCCCGCUAGGCUCCA",
                    "sequence_iupac": "GGGGCCUUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCAGCGGUUCGAUCCCGCUAGGCUCCA",
                    "length": 73,
                    "number_of_modifications": 4,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C701H805N281O518P73S",
                    "molecular_weight": 23747.74463,
                    "charge": -73,
                    "canonical_formula": "C694H790N279O515P73",
                    "canonical_molecular_weight": 23540.47663,
                    "canonical_charge": -74,
                    "extra_formula": "C7H15N2O3S",
                    "extra_molecular_weight": 207.268,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:00.788Z"
                    },
                    "ncbi_taxonomy_id": 1358
                },
                {
                    "anticodon": "GGC",
                    "organism": "Streptomyces griseus",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCUAUAGCUBAGUDGGDAGAGCGCCUGCAUGGCAUGCAGGAG7UCAGGAGUUCA\"UUCUCCUUAGCUCCACAA",
                    "sequence_bpforms": "GGGGCUAUAGCU{0C}AGU{8U}GG{8U}AGAGCGCCUGCAUGGCAUGCAGGAG{7G}UCAGGAGUUCA{1A}UUCUCCUUAGCUCCACAA",
                    "sequence_iupac": "GGGGCUAUAGCUCAGUUGGUAGAGCGCCUGCAUGGCAUGCAGGAGGUCAGGAGUUCAAUUCUCCUUAGCUCCACAA",
                    "length": 76,
                    "number_of_modifications": 5,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C727H832N290O534P76",
                    "molecular_weight": 24530.15491,
                    "charge": -76,
                    "canonical_formula": "C724H819N290O534P76",
                    "canonical_molecular_weight": 24481.01791,
                    "canonical_charge": -77,
                    "extra_formula": "C3H13",
                    "extra_molecular_weight": 49.137,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:02.506Z"
                    },
                    "ncbi_taxonomy_id": 1911
                },
                {
                    "anticodon": "GGC",
                    "organism": "Streptomyces griseus",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCUAUAGCUBAGUDGGDAGAGCGCCUGCAUGGCAUGCAGGAG7UCAGGAGUUCA\"UUCUCCUUAGCUCCA",
                    "sequence_bpforms": "GGGGCUAUAGCU{0C}AGU{8U}GG{8U}AGAGCGCCUGCAUGGCAUGCAGGAG{7G}UCAGGAGUUCA{1A}UUCUCCUUAGCUCCA",
                    "sequence_iupac": "GGGGCUAUAGCUCAGUUGGUAGAGCGCCUGCAUGGCAUGCAGGAGGUCAGGAGUUCAAUUCUCCUUAGCUCCA",
                    "length": 73,
                    "number_of_modifications": 5,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 1,
                    "number_of_modified_u": 2,
                    "formula": "C698H799N277O515P73",
                    "molecular_weight": 23569.57863,
                    "charge": -73,
                    "canonical_formula": "C695H786N277O515P73",
                    "canonical_molecular_weight": 23520.44163,
                    "canonical_charge": -74,
                    "extra_formula": "C3H13",
                    "extra_molecular_weight": 49.137,
                    "extra_charge": 1,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:02.535Z"
                    },
                    "ncbi_taxonomy_id": 1911
                },
                {
                    "anticodon": "CGC",
                    "organism": "Streptomyces griseus",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGCCUGUGGCGCAGUCUGGDAGCGCACCUCGUUCGCAUCGAGGGGGUCUGGGGUUCA\"AUCCCCACAGGUCCA",
                    "sequence_bpforms": "GGGCCUGUGGCGCAGUCUGG{8U}AGCGCACCUCGUUCGCAUCGAGGGGGUCUGGGGUUCA{1A}AUCCCCACAGGUCCA",
                    "sequence_iupac": "GGGCCUGUGGCGCAGUCUGGUAGCGCACCUCGUUCGCAUCGAGGGGGUCUGGGGUUCAAAUCCCCACAGGUCCA",
                    "length": 74,
                    "number_of_modifications": 2,
                    "number_of_modified_a": 1,
                    "number_of_modified_c": 0,
                    "number_of_modified_g": 0,
                    "number_of_modified_u": 1,
                    "formula": "C704H804N281O523P74",
                    "molecular_weight": 23861.67839,
                    "charge": -75,
                    "canonical_formula": "C703H800N281O523P74",
                    "canonical_molecular_weight": 23845.63539,
                    "canonical_charge": -75,
                    "extra_formula": "CH4",
                    "extra_molecular_weight": 16.043,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:02.564Z"
                    },
                    "ncbi_taxonomy_id": 1911
                },
                {
                    "anticodon": "UGC",
                    "organism": "Streptomyces griseus",
                    "organellum": "prokaryotic cytosol",
                    "sequence_modomics": "GGGGCCUUAGCUBAGUDGGDAGAGCGCUGCCUUUGCAAGGCAGAUGUCAGGAGUUCGAAUCUCCUAGGCUCCA",
                    "sequence_bpforms": "GGGGCCUUAGCU{0C}AGU{8U}GG{8U}AGAGCGCUGCCUUUGCAAGGCAGAUGUCAGGAGUUCGAAUCUCCUAGGCUCCA",
                    "sequence_iupac": "GGGGCCUUAGCUCAGUUGGUAGAGCGCUGCCUUUGCAAGGCAGAUGUCAGGAGUUCGAAUCUCCUAGGCUCCA",
                    "length": 73,
                    "number_of_modifications": 3,
                    "number_of_modified_a": 0,
                    "number_of_modified_c": 1,
                    "number_of_modified_g": 0,
                    "number_of_modified_u": 2,
                    "formula": "C695H792N275O516P73",
                    "molecular_weight": 23514.47463,
                    "charge": -74,
                    "canonical_formula": "C694H786N275O516P73",
                    "canonical_molecular_weight": 23496.41563,
                    "canonical_charge": -74,
                    "extra_formula": "CH6",
                    "extra_molecular_weight": 18.059,
                    "extra_charge": 0,
                    "bpforms_errors": null,
                    "reference": {
                        "doi": "10.1093/nar/gkx1030"
                    },
                    "last_modified": {
                        "$date": "2020-04-28T23:34:02.594Z"
                    },
                    "ncbi_taxonomy_id": 1911
                }
            ]
        }
        results = self.src.build_rna_modification_observation(obj)
        print(results[0])