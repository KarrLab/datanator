import unittest
from datanator_query_python.config import config
from datanator.schema_2 import transform_metabolites_meta


class TestTMM(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        conf = config.SchemaMigration()
        cls.obj = {
                    "InChI_Key": "TYEYBOSBBBHJIV-UHFFFAOYSA-N",
                    "average_molecular_weight": "102.0886",
                    "biocyc_id": "2-OXOBUTANOATE",
                    "cas_registry_number": "600-18-0",
                    "cellular_locations": [
                        {
                            "reference": [
                                "YMDB"
                            ],
                            "cellular_location": "mitochondrion"
                        },
                        {
                            "reference": [
                                "YMDB"
                            ],
                            "cellular_location": "cytoplasm"
                        }
                    ],
                    "chebi_id": "30831",
                    "chemical_formula": "C4H6O3",
                    "chemspider_id": "57",
                    "description": "Alpha-Ketoglutaric acid (2-oxobutanoate, 2-keto-butyrate) is a key intermediate in the Krebs cycle. It is involved in the metabolism of many amino acids (glycine, methionine, valine, leucine, serine, threonine, isoleucine) as well as propanoate metabolism and C-5 branched dibasic acid metabolism. It can be converted to propionyl-CoA and thus enter the citric acid cycle.",
                    "hmdb_id": "HMDB00005",
                    "inchi": "InChI=1S/C4H6O3/c1-2-3(5)4(6)7/h2H2,1H3,(H,6,7)",
                    "kegg_id": "C00109",
                    "kinlaw_id": [],
                    "m2m_id": "M2MDB000001",
                    "name": "2-Ketobutyric acid",
                    "pathways": {
                        "pathway": [
                            {
                                "name": "Propanoate metabolism",
                                "description": None,
                                "pathwhiz_id": None,
                                "kegg_map_id": "00640",
                                "subject": None
                            },
                            {
                                "name": "Cysteine and methionine metabolism",
                                "description": None,
                                "pathwhiz_id": None,
                                "kegg_map_id": "00270",
                                "subject": None
                            },
                            {
                                "name": "Valine, leucine and isoleucine biosynthesis",
                                "description": None,
                                "pathwhiz_id": None,
                                "kegg_map_id": "00290",
                                "subject": None
                            },
                            {
                                "name": "Glycine, serine and threonine metabolism",
                                "description": None,
                                "pathwhiz_id": None,
                                "kegg_map_id": "00260",
                                "subject": None
                            },
                            {
                                "name": "Cysteine Metabolism",
                                "description": "The biosynthesis of cysteine begins with aspartate being phosphorylated into L-aspartyl-4-phosphate through an ATP driven aspartate kinase. L-aspartyl-4-phosphate is then catabolized through an NADPH dependent Aspartic Beta-Semiladehyde dehydrogenase resulting in the release of L-aspartate semialdehyde which is transformed into a homoserine through a Homoserine dehydrogenase. Homeserine in turn is acetylated through a homoserine O-trans-acetylase resulting in the release of O-acetyl-L-homoserine.\nThe latter compound interacts with hydrogen sulfide through a O-acetylhomoserine (thiol)-lyase resulting in the release of L-homocysteine. L-homocysteine reacts with serine through a cystathionine beta synthase resulting in the release of water and L-cystathionine. This compound in turn can be turned into cysteine by reacting with water through a cystathionine gama-lyase. Cysteine can be turned back to L-cystathionine by reacting with a acetyl-L-homoserine spontaneously, thus resulting in L-cystathionine.\nCysteine can also be degraded by reacting with a cystathionine gamma lyase resulting in the release of hydrogen sulfide, a hydrogen ion and 2-aminoprop-2-enoate which can spontaneously be converted into 2-iminopropanoate and further degraded into pyruvic acid.",
                                "pathwhiz_id": "PW002383",
                                "kegg_map_id": None,
                                "subject": "Metabolic"
                            },
                            {
                                "name": "Methionine metabolism and salvage",
                                "description": "The biosynthesis of Methionine begins with partate being phosphorylated into L-aspartyl-4-phosphate through an ATP driven aspartate kinase. L-aspartyl-4-phosphate is then catabolized through an NADPH dependent Aspartic Beta-Semiladehyde dehydrogenase resulting in the release of L-aspartate semialdehyde which is transformed into a homoserine through a Homoserine dehydrogenase. Homeserine in turn is acetylated through a homoserine O-trans-acetylase resulting in the release of O-acetyl-L-homoserine. \nThe latter compound interacts with hydrogen sulfide through a O-acetylhomoserine (thiol)-lyase resulting in the release of L-homocysteineThe latter compound then reacts with 5-methylterahydropteroyltri-L-glutamate through a N5-methyltetrahydropteroyltrigluatamate homocysteine methyl transferase resulting in the release of a tetrahydropterooyltri-l-glutamate and Methionine. \nThe degradation of methionine begins with methionine being used to synthesize S-adenosylmethionine through a S-adenosylmethionine synthetase. The s-adenosylmethionine reacts with a demethylated methyl donor resulting in the release of a methylated methyl donor, a hydrogen ion and a S-adenosylhomocysteine. The latter compound the reacts with a S-adenosyl-L-homocysteine hydrolase resulting in the release of adenosine and homocysteine where the cycle can begin again.\nThe salvage of Methionine begins with S-methyl-5'-thioadenosine (a product of spermine biosynthesis) being phosphorylated through a 5-methylthioadenosine phosphorylase resulting in the release of adenine and S-methyl-5-thio-alpha-D-ribose 1-phosphate. This last compound is isomerized into 5-methylthioribulose 1-phosphate. The latter compound is then dehydrated through a methylthioribulose 1-phosphate dehydratase resulting in 5-(methylthio)-2,3-dioxopentyl 1-phosphate. This resulting compound is then dephosphorylated through a 2,3-dioxomethiopentane-1-phosphate enolase/phosphatase resulting in a  1,2-dihydroxy-5-(methylthio)pent-1-en-3-one. This latter compound can react spontaneously or through a acireductone dioxygenase resulting in the release of a  2-oxo-4-methylthiobutanoate. This latter compound is then turned into methionine through a aromatic amino acid aminotransferase II",
                                "pathwhiz_id": "PW002384",
                                "kegg_map_id": None,
                                "subject": "Metabolic"
                            },
                            {
                                "name": "Selenocompound metabolism",
                                "description": "The the metabolism of selenocompounds starts with Selenocysteine being metabolized by a CGS resulting in the release of Seleno-cystathionine. The resulting compound is metabolized by a CBL resulting in the release of selenohomocysteine. The resulting compound reacts with MET resulting in the release of a seleno-methionine. Selenomethionine can be either metabolized into Seleno-methionyl-tRNA or a Methyl-selenol. \nMethyl-selenol can also be the result of Methyl-selenic acid reacting with a thioredoxin reductase or Se-methyl-selenocysteine reacting through a CTH.",
                                "pathwhiz_id": "PW002472",
                                "kegg_map_id": "00450",
                                "subject": "Metabolic"
                            },
                            {
                                "name": "Sulfur metabolism",
                                "description": None,
                                "pathwhiz_id": "PW002483",
                                "kegg_map_id": "00920",
                                "subject": "Metabolic"
                            },
                            {
                                "name": "isoleucine biosynthesis",
                                "description": "Isoleucine biosynthesis begins with L-threonine from the threonine biosynthesis pathway. L-threonine interacts with a threonine dehydratase biosynthetic releasing water, a hydrogen ion and (2Z)-2-aminobut-2-enoate. This compound is isomerized into a 2-iminobutanoate which interacts with water and a hydrogen ion spontaneously, resulting in the release of ammonium and 2-ketobutyric acid. This compound reacts with pyruvic acid and hydrogen ion through an acetohydroxybutanoate synthase / acetolactate synthase 2 resulting in carbon dioxide and (S)-2-Aceto-2-hydroxybutanoic acid. The latter compound is reduced by an NADPH driven acetohydroxy acid isomeroreductase releasing NADP and acetohydroxy acid isomeroreductase. The latter compound is dehydrated by a dihydroxy acid dehydratase resulting in 3-methyl-2-oxovaleric acid.This compound reacts in a reversible reaction with L-glutamic acid through a Branched-chain-amino-acid aminotransferase resulting in oxoglutaric acid and L-isoleucine.",
                                "pathwhiz_id": "PW002476",
                                "kegg_map_id": None,
                                "subject": "Metabolic"
                            }
                        ]
                    },
                    "property": [
                        {
                            "kind": "logp",
                            "value": "0.77",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "pka_strongest_acidic",
                            "value": "3.19",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "pka_strongest_basic",
                            "value": "-9.7",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "iupac",
                            "value": "2-oxobutanoic acid",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "average_mass",
                            "value": "102.0886",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "mono_mass",
                            "value": "102.031694058",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "smiles",
                            "value": "[H]OC(=O)C(=O)C([H])([H])C([H])([H])[H]",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "formula",
                            "value": "C4H6O3",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "inchi",
                            "value": "InChI=1S/C4H6O3/c1-2-3(5)4(6)7/h2H2,1H3,(H,6,7)",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "inchikey",
                            "value": "TYEYBOSBBBHJIV-UHFFFAOYSA-N",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "polar_surface_area",
                            "value": "54.37",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "refractivity",
                            "value": "22.62",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "polarizability",
                            "value": "9.2",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "rotatable_bond_count",
                            "value": "2",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "acceptor_count",
                            "value": "3",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "donor_count",
                            "value": "1",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "physiological_charge",
                            "value": "-1",
                            "source": "ChemAxon"
                        },
                        {
                            "kind": "formal_charge",
                            "value": "0",
                            "source": "ChemAxon"
                        }
                    ],
                    "pubchem_compound_id": "58",
                    "reaction_participants": [],
                    "similar_compounds": [
                        {
                            "inchikey": "DNOPJXBPONYBLB-UHFFFAOYSA-N",
                            "similarity_score": 0.875
                        },
                        {
                            "inchikey": "NGEWQZIDQIYUNV-UHFFFAOYSA-N",
                            "similarity_score": 0.867
                        },
                        {
                            "inchikey": "AFENDNXGAFYKQO-UHFFFAOYSA-N",
                            "similarity_score": 0.867
                        },
                        {
                            "inchikey": "WLAMNBDJUVNPJU-UHFFFAOYSA-N",
                            "similarity_score": 0.857
                        },
                        {
                            "inchikey": "GWYFCOCPABKNJV-UHFFFAOYSA-N",
                            "similarity_score": 0.857
                        },
                        {
                            "inchikey": "FERIUCNNQQJTOY-UHFFFAOYSA-N",
                            "similarity_score": 0.857
                        },
                        {
                            "inchikey": "PKVVTUWHANFMQC-UHFFFAOYSA-N",
                            "similarity_score": 0.824
                        },
                        {
                            "inchikey": "JVQYSWDUAOAHFM-BYPYZUCNSA-M",
                            "similarity_score": 0.824
                        },
                        {
                            "inchikey": "JVQYSWDUAOAHFM-BYPYZUCNSA-N",
                            "similarity_score": 0.824
                        },
                        {
                            "inchikey": "PKVVTUWHANFMQC-UHFFFAOYSA-M",
                            "similarity_score": 0.824
                        },
                        {
                            "inchikey": "BKAJNAXTPSGJCU-UHFFFAOYSA-N",
                            "similarity_score": 0.824
                        },
                        {
                            "inchikey": "JVQYSWDUAOAHFM-UHFFFAOYSA-N",
                            "similarity_score": 0.824
                        },
                        {
                            "inchikey": "LCTONWCANYUPML-UHFFFAOYSA-N",
                            "similarity_score": 0.786
                        },
                        {
                            "inchikey": "JTEYKUFKXGDTEU-UHFFFAOYSA-N",
                            "similarity_score": 0.765
                        },
                        {
                            "inchikey": "JTEYKUFKXGDTEU-VKHMYHEASA-N",
                            "similarity_score": 0.765
                        },
                        {
                            "inchikey": "NMDWGEGFJUBKLB-YFKPBYRVSA-N",
                            "similarity_score": 0.765
                        },
                        {
                            "inchikey": "JTEYKUFKXGDTEU-VKHMYHEASA-M",
                            "similarity_score": 0.765
                        },
                        {
                            "inchikey": "NMDWGEGFJUBKLB-UHFFFAOYSA-N",
                            "similarity_score": 0.765
                        },
                        {
                            "inchikey": "LOUGYXZSURQALL-PWNYCUMCSA-N",
                            "similarity_score": 0.765
                        },
                        {
                            "inchikey": "UIUJIQZEACWQSV-UHFFFAOYSA-N",
                            "similarity_score": 0.75
                        },
                        {
                            "inchikey": "WHBMMWSBFZVSSR-UHFFFAOYSA-N",
                            "similarity_score": 0.75
                        },
                        {
                            "inchikey": "SKCYVGUCBRYGTE-UHFFFAOYSA-N",
                            "similarity_score": 0.75
                        },
                        {
                            "inchikey": "WHBMMWSBFZVSSR-GSVOUGTGSA-N",
                            "similarity_score": 0.75
                        },
                        {
                            "inchikey": "WDJHALXBUFZDSR-UHFFFAOYSA-N",
                            "similarity_score": 0.75
                        },
                        {
                            "inchikey": "SJZRECIVHVDYJC-UHFFFAOYSA-N",
                            "similarity_score": 0.75
                        },
                        {
                            "inchikey": "YJVOWRAWFXRESP-ZCFIWIBFSA-N",
                            "similarity_score": 0.737
                        },
                        {
                            "inchikey": "KHPXUQMNIQBQEV-UHFFFAOYSA-N",
                            "similarity_score": 0.737
                        },
                        {
                            "inchikey": "VUQLHQFKACOHNZ-LURJTMIESA-N",
                            "similarity_score": 0.722
                        },
                        {
                            "inchikey": "LVRFTAZAXQPQHI-UHFFFAOYSA-N",
                            "similarity_score": 0.722
                        },
                        {
                            "inchikey": "OTOIIPJYVQJATP-UHFFFAOYSA-M",
                            "similarity_score": 0.722
                        },
                        {
                            "inchikey": "OTOIIPJYVQJATP-BYPYZUCNSA-N",
                            "similarity_score": 0.722
                        },
                        {
                            "inchikey": "JRHWHSJDIILJAT-UHFFFAOYSA-N",
                            "similarity_score": 0.722
                        },
                        {
                            "inchikey": "RILPIWOPNGRASR-RFZPGFLSSA-N",
                            "similarity_score": 0.722
                        },
                        {
                            "inchikey": "ROWKJAVDOGWPAT-GSVOUGTGSA-N",
                            "similarity_score": 0.714
                        },
                        {
                            "inchikey": "ROWKJAVDOGWPAT-UHFFFAOYSA-N",
                            "similarity_score": 0.714
                        },
                        {
                            "inchikey": "NQPDZGIKBAWPEJ-UHFFFAOYSA-N",
                            "similarity_score": 0.706
                        },
                        {
                            "inchikey": "KDYFGRWQOYBRFD-UHFFFAOYSA-N",
                            "similarity_score": 0.706
                        },
                        {
                            "inchikey": "LBUDVZDSWKZABS-UHFFFAOYSA-M",
                            "similarity_score": 0.706
                        },
                        {
                            "inchikey": "HFKQINMYQUXOCH-UHFFFAOYSA-M",
                            "similarity_score": 0.7
                        },
                        {
                            "inchikey": "HHDDCCUIIUWNGJ-UHFFFAOYSA-N",
                            "similarity_score": 0.688
                        },
                        {
                            "inchikey": "UYTRITJAZOPLCZ-BYPYZUCNSA-N",
                            "similarity_score": 0.684
                        },
                        {
                            "inchikey": "UYTRITJAZOPLCZ-SCSAIBSYSA-N",
                            "similarity_score": 0.684
                        },
                        {
                            "inchikey": "UYTRITJAZOPLCZ-UHFFFAOYSA-N",
                            "similarity_score": 0.684
                        },
                        {
                            "inchikey": "JVTAAEKCZFNVCJ-UWTATZPHSA-N",
                            "similarity_score": 0.667
                        },
                        {
                            "inchikey": "NOXRYJAWRSNUJD-UHFFFAOYSA-N",
                            "similarity_score": 0.667
                        },
                        {
                            "inchikey": "RMHHUKGVZFVHED-UHFFFAOYSA-N",
                            "similarity_score": 0.667
                        },
                        {
                            "inchikey": "JVTAAEKCZFNVCJ-REOHCLBHSA-N",
                            "similarity_score": 0.667
                        },
                        {
                            "inchikey": "BJEPYKJPYRNKOW-UWTATZPHSA-N",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "XFTRTWQBIOMVPK-UHFFFAOYSA-N",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "BJEPYKJPYRNKOW-REOHCLBHSA-L",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "PDGXJDXVGMHUIR-UJURSFKZSA-M",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "FEWJPZIEWOKRBE-JCYAYHJZSA-N",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "NPYQJIHHTGFBLN-STHAYSLISA-N",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "PDGXJDXVGMHUIR-UJURSFKZSA-N",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "BJEPYKJPYRNKOW-UHFFFAOYSA-N",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "FEWJPZIEWOKRBE-UHFFFAOYSA-N",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "BJEPYKJPYRNKOW-REOHCLBHSA-N",
                            "similarity_score": 0.65
                        },
                        {
                            "inchikey": "XBDQKXXYIPTUBI-UHFFFAOYSA-N",
                            "similarity_score": 0.643
                        },
                        {
                            "inchikey": "KQNPFQTWMSNSAP-UHFFFAOYSA-N",
                            "similarity_score": 0.643
                        },
                        {
                            "inchikey": "KPGXRSRHYNQIFN-UHFFFAOYSA-N",
                            "similarity_score": 0.636
                        },
                        {
                            "inchikey": "FGSBNBBHOZHUBO-UHFFFAOYSA-N",
                            "similarity_score": 0.636
                        },
                        {
                            "inchikey": "KPGXRSRHYNQIFN-UHFFFAOYSA-L",
                            "similarity_score": 0.636
                        },
                        {
                            "inchikey": "HIIZAGQWABAMRR-BYPYZUCNSA-N",
                            "similarity_score": 0.636
                        },
                        {
                            "inchikey": "KZSNJWFQEVHDMF-BYPYZUCNSA-N",
                            "similarity_score": 0.632
                        },
                        {
                            "inchikey": "JOOXCMJARBKPKM-UHFFFAOYSA-N",
                            "similarity_score": 0.632
                        },
                        {
                            "inchikey": "UUIQMZJEGPQKFD-UHFFFAOYSA-N",
                            "similarity_score": 0.632
                        },
                        {
                            "inchikey": "WRBRCYPPGUCRHW-UHFFFAOYSA-M",
                            "similarity_score": 0.632
                        },
                        {
                            "inchikey": "KZSNJWFQEVHDMF-SCSAIBSYSA-N",
                            "similarity_score": 0.632
                        },
                        {
                            "inchikey": "QWCKQJZIFLGMSD-VKHMYHEASA-N",
                            "similarity_score": 0.632
                        },
                        {
                            "inchikey": "UQIGQRSJIKIPKZ-GSVOUGTGSA-N",
                            "similarity_score": 0.609
                        },
                        {
                            "inchikey": "ALFQPWXBAWHVDP-UHFFFAOYSA-N",
                            "similarity_score": 0.609
                        },
                        {
                            "inchikey": "AQYCMVICBNBXNA-UHFFFAOYSA-N",
                            "similarity_score": 0.6
                        },
                        {
                            "inchikey": "BTCSSZJGUNDROE-UHFFFAOYSA-N",
                            "similarity_score": 0.6
                        },
                        {
                            "inchikey": "FUZZWVXGSFPDMH-UHFFFAOYSA-M",
                            "similarity_score": 0.6
                        },
                        {
                            "inchikey": "JFCQEDHGNNZCLN-UHFFFAOYSA-N",
                            "similarity_score": 0.6
                        },
                        {
                            "inchikey": "FUZZWVXGSFPDMH-UHFFFAOYSA-N",
                            "similarity_score": 0.6
                        },
                        {
                            "inchikey": "WDRISBUVHBMJEF-MROZADKFSA-N",
                            "similarity_score": 0.6
                        },
                        {
                            "inchikey": "RBNPOMFGQQGHHO-UHFFFAOYSA-N",
                            "similarity_score": 0.588
                        },
                        {
                            "inchikey": "QWBAFPFNGRFSFB-UHFFFAOYSA-N",
                            "similarity_score": 0.588
                        },
                        {
                            "inchikey": "RBNPOMFGQQGHHO-UWTATZPHSA-N",
                            "similarity_score": 0.588
                        },
                        {
                            "inchikey": "QIQXTHQIDYTFRH-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "SAUCHDKDCUROAO-VKHMYHEASA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "TUNFSRHWOTWDNC-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "AYFVYJQAPQTCCC-PWNYCUMCSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "GHVNFZFCNZKVNT-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "BSABBBMNWQWLLU-VKHMYHEASA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "AYFVYJQAPQTCCC-GBXIJSLDSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "TUNFSRHWOTWDNC-UHFFFAOYSA-M",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "POULHZVOKOAJMA-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "BDJRBEYXGGNYIS-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "QIQXTHQIDYTFRH-UHFFFAOYSA-M",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "GHVNFZFCNZKVNT-UHFFFAOYSA-M",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "BSABBBMNWQWLLU-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "FBUKVWPVBMHYJY-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "QZZGJDVWLFXDLK-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "WWZKQHOCKIZLMA-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "XMHIUKTWLZUKEX-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "MNWFXJYAOYHMED-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "POULHZVOKOAJMA-UHFFFAOYSA-M",
                            "similarity_score": 0.571
                        },
                        {
                            "inchikey": "ZDPHROOEEOARMN-UHFFFAOYSA-N",
                            "similarity_score": 0.571
                        }
                    ],
                    "smiles": "[H]OC(=O)C(=O)C([H])([H])C([H])([H])[H]",
                    "synonyms": [
                        "2-Ketobutanoate",
                        "2-Ketobutanoic acid",
                        "2-ketobutyrate",
                        "2-ketobutyric acid",
                        "2-Oxo-Butanoate",
                        "2-Oxo-Butanoic acid",
                        "2-oxo-Butyrate",
                        "2-oxo-Butyric acid",
                        "2-Oxo-n-butyrate",
                        "2-Oxo-n-butyric acid",
                        "2-Oxobutanoate",
                        "2-Oxobutanoic acid",
                        "2-Oxobutyrate",
                        "2-Oxobutyric acid",
                        "3-Methylpyruvate",
                        "3-Methylpyruvic acid",
                        "a-keto-n-Butyrate",
                        "a-keto-n-Butyric acid",
                        "a-Ketobutyrate",
                        "a-Ketobutyric acid",
                        "a-Oxo-n-butyrate",
                        "a-Oxo-n-butyric acid",
                        "a-Oxobutyrate",
                        "a-Oxobutyric acid",
                        "alpha-Keto-n-butyrate",
                        "alpha-Keto-n-butyric acid",
                        "alpha-Ketobutric acid",
                        "alpha-Ketobutyrate",
                        "alpha-Ketobutyric acid",
                        "alpha-Oxo-n-butyrate",
                        "alpha-Oxo-n-butyric acid",
                        "alpha-Oxobutyrate",
                        "alpha-Oxobutyric acid",
                        "Butanoic acid, 2-oxo-",
                        "Butyric acid, 2-oxo-",
                        "Formic acid, propionyl-",
                        "Ketobutyrate",
                        "methyl-Pyruvate",
                        "methyl-Pyruvic acid",
                        "Oxobutyrate",
                        "propionyl-formate",
                        "propionyl-formic acid",
                        "Pyruvic acid, methyl-"
                    ],
                    "ymdb_id": "YMDB00071"
                }
        cls.src = transform_metabolites_meta.TransformMetabolitesMeta(MongoDB=conf.SERVER,
                                                                      db="test",
                                                                      username=conf.USERNAME,
                                                                      password=conf.PASSWORD,
                                                                      max_entries=20)

    @classmethod
    def tearDownClass(cls):
        pass

    @unittest.skip("for now")
    def test_process_docs(self):
        self.src.process_docs()
    
    @unittest.skip("for now")
    def test_build_entity(self):
        entity = self.src.build_entity(self.obj)
        print(entity)

    @unittest.skip("for now")
    def test_build_obs(self):
        a, b = self.src.build_obs(self.obj)
        print(b)