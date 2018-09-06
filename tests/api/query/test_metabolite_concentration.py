""" Tests of metabolite concentration queries

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-09-27
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.api.query import metabolite_concentrations
from kinetic_datanator.core import models, common_schema
import tempfile
import shutil
import unittest
import mock

class TestMetaboliteConcentrationQuery(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cache_dirname = tempfile.mkdtemp()
        flk = common_schema.CommonSchema(cache_dirname=self.cache_dirname)
        self.proline = flk.session.query(models.Metabolite).filter_by(metabolite_name = 'L-Proline').first()
        self.uridine_tp = flk.session.query(models.Metabolite).filter_by(metabolite_name = 'Uridine triphosphate').first()

        self.q = metabolite_concentrations.MetaboliteConcentrationQuery(cache_dirname=self.cache_dirname, include_variants=True)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.cache_dirname)

    def test_filter_observed_results(self):

        obs = self.q.run(self.proline).observed_results

        self.assertEqual(set(c.value for c in obs), set([385.0, 451.0, 361.0, 143.0, 550.0, 531.67]))
        self.assertEqual(set(c.error for c in obs), set([0.0,0.0,0.0,0.0,45.0,11.37]))
        self.assertEqual(set(c.units for c in obs), set(['uM']))
        self.assertEqual(set(c.metadata.genetics.taxon for c in obs) , set(['Escherichia coli']))
        self.assertEqual(set(c.metadata.genetics.variation for c in obs), set(['BL21 DE3', 'K12 NCM3722', 'BW25113']))
        self.assertEqual(set(c.metadata.environment.temperature for c in obs), set([37.0]))
        self.assertEqual(set(c.metadata.environment.ph for c in obs), set([None]))
        self.assertEqual(set(c.metadata.environment.media for c in obs), set(['Gutnick minimal complete medium (4.7 g/L KH2PO4; 13.5 g/L K2HPO4; 1 g/L K2SO4; 0.1 g/L MgSO4-7H2O; 10 mM NH4Cl) with 4 g/L glucose', 'Luria-Bertani (LB) media', 'Gutnick minimal complete medium (4.7 g/L KH2PO4; 13.5 g/L K2HPO4; 1 g/L K2SO4; 0.1 g/L MgSO4-7H2O; 10 mM NH4Cl) with 4 g/L glycerol', '48 mM Na2HPO4, 22 mM KH2PO4, 10 mM NaCl, 45 mM (NH4)2SO4, supplemented with 1 mM MgSO4, 1 mg/l thiamine·HCl, 5.6 mg/l CaCl2, 8 mg/l FeCl3, 1 mg/l MnCl2·4H2O, 1.7 mg/l ZnCl2, 0.43 mg/l CuCl2·2H2O, 0.6 mg/l CoCl2·2H2O and 0.6 mg/l Na2MoO4·2H2O.  4 g/L Gluco', 'Gutnick minimal complete medium (4.7 g/L KH2PO4; 13.5 g/L K2HPO4; 1 g/L K2SO4; 0.1 g/L MgSO4-7H2O; 10 mM NH4Cl) with 4 g/L acetate']))
        self.assertEqual(set(c.observable.specie.name for c in obs), set(['L-Proline']))
        self.assertEqual(set(c.observable.specie.structure for c in obs), set(['InChI=1S/C5H9NO2/c7-5(8)4-2-1-3-6-4/h4,6H,1-3H2,(H,7,8)/t4-/m0/s1']))

        references = [c.observable.specie.cross_references for c in obs]
        self.assertEqual(references[0], references[1])


        obs = self.q.run(self.uridine_tp).observed_results

        self.assertEqual(set(c.value for c in obs), set([2370.0, 3990.0, 8290.0, 663.0]))
        self.assertEqual(set(c.error for c in obs), set([0.0]))
        self.assertEqual(set(c.units for c in obs), set(['uM']))
        self.assertEqual(set(c.metadata.genetics.taxon for c in obs), set(['Escherichia coli']))
        self.assertEqual(set(c.metadata.genetics.variation for c in obs), set(['K12 NCM3722', 'BW25113']))
        self.assertEqual(set(c.metadata.environment.temperature for c in obs), set([37.0]))
        self.assertEqual(set(c.metadata.environment.ph for c in obs), set([None]))
        self.assertEqual(set(c.metadata.environment.media for c in obs), set(['Gutnick minimal complete medium (4.7 g/L KH2PO4; 13.5 g/L K2HPO4; 1 g/L K2SO4; 0.1 g/L MgSO4-7H2O; 10 mM NH4Cl) with 4 g/L glycerol', 'Gutnick minimal complete medium (4.7 g/L KH2PO4; 13.5 g/L K2HPO4; 1 g/L K2SO4; 0.1 g/L MgSO4-7H2O; 10 mM NH4Cl) with 4 g/L glucose', '48 mM Na2HPO4, 22 mM KH2PO4, 10 mM NaCl, 45 mM (NH4)2SO4, supplemented with 1 mM MgSO4, 1 mg/l thiamine·HCl, 5.6 mg/l CaCl2, 8 mg/l FeCl3, 1 mg/l MnCl2·4H2O, 1.7 mg/l ZnCl2, 0.43 mg/l CuCl2·2H2O, 0.6 mg/l CoCl2·2H2O and 0.6 mg/l Na2MoO4·2H2O.  4 g/L Gluco', 'Gutnick minimal complete medium (4.7 g/L KH2PO4; 13.5 g/L K2HPO4; 1 g/L K2SO4; 0.1 g/L MgSO4-7H2O; 10 mM NH4Cl) with 4 g/L acetate']))
        self.assertEqual(set(c.observable.specie.name for c in obs), set(['Uridine triphosphate']))
        self.assertEqual(set(c.observable.specie.structure for c in obs), set(['InChI=1S/C9H15N2O15P3/c12-5-1-2-11(9(15)10-5)8-7(14)6(13)4(24-8)3-23-28(19,20)26-29(21,22)25-27(16,17)18/h1-2,4,6-8,13-14H,3H2,(H,19,20)(H,21,22)(H,10,12,15)(H2,16,17,18)/t4-,6-,7-,8-/m1/s1']))

        references = [c.observable.specie.cross_references for c in obs]
        self.assertEqual(references[0], references[1])

    def test_get_concentration_by_structure(self):

        concentrations = self.q.get_concentration_by_structure(self.proline.structure._value_inchi, only_formula_and_connectivity=False).all()

        self.assertEqual(set(c.value for c in concentrations), set([385.0, 451.0, 361.0, 143.0, 550.0, 531.67]))
        self.assertEqual(set(c._metadata.taxon[0].name for c in concentrations), set(['Escherichia coli']))
        self.assertEqual(set(c._metadata.cell_compartment[0].name for c in concentrations), set(['Cytosol']))

        concentrations = self.q.get_concentration_by_structure(self.uridine_tp.structure._value_inchi, only_formula_and_connectivity=True)
        self.assertEqual(set(c.value for c in concentrations), set([8290.0, 3990.0, 2370.0, 663.0]))
