""" Tests of Ezyme

:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-04
:Copyright: 2017, Karr Lab
:License: MIT
"""

from attrdict import AttrDict
from kinetic_datanator.data_source import ezyme
from kinetic_datanator.core import data_model
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import warning_util
import unittest

warning_util.disable_warnings()


class TestEzyme(unittest.TestCase):
    molecules = {
        'adp': 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1',
        'amet': 'C[S+](CCC([NH3+])C([O-])=O)CC1OC(C(O)C1O)N1C=NC2=C(N)N=CN=C12',
        'atp': 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1',
        'dr1p': 'OCC1OC(CC1O)OP([O-])([O-])=O',
        'dr5p': 'OC1CC(O)C(COP([O-])([O-])=O)O1',
        'e1dG': 'CCN1C(N)=NC2=C(N=CN2C2CC(O)C(CO)O2)C1=O',
        'e1dGMP': 'CCN1C(N)=NC2=C(N=CN2C2CC(O)C(COP([O-])([O-])=O)O2)C1=O',
        'f1p': 'OCC1OC(O)(COP([O-])([O-])=O)C(O)C1O',
        'gly': '[NH3+]CC([O-])=O',
        'glyceraldehyde': 'OCC(O)C=O',
        'h': '[H+]',
        'h2o': 'O',
        'ileile': 'CCC(C)C([NH3+])C(=O)NC(C(C)CC)C([O-])=O',
        'leu': 'CC(C)CC([NH3+])C([O-])=O',
        'leuleu': 'CC(C)CC([NH3+])C(=O)NC(CC(C)C)C([O-])=O',
        'met': 'CSCCC([NH3+])C([O-])=O',
        'metthf': 'NC1=NC(=O)C2=C(NCC3CN(CN23)C2=CC=C(C=C2)C(=O)NC(CCC([O-])=O)C([O-])=O)N1',
        'pep': '[O-]C(=O)C(=C)OP([O-])([O-])=O',
        'pi': 'OP([O-])([O-])=O',
        'ppi': '[O-]P([O-])(=O)OP([O-])([O-])=O',
        'pyr': 'CC(=O)C([O-])=O',
        'r5p': 'OC(COP([O-])([O-])=O)C(O)C(O)C=O',
        's7p': 'OCC(=O)C(O)C(O)C(O)C(O)COP([O-])([O-])=O',
        'ser': '[NH3+]C(CO)C([O-])=O',
        't3p1': 'OC(COP([O-])([O-])=O)C=O',
        't3p2': 'OCC(=O)COP([O-])([O-])=O',
        'thf': 'NC1=NC(=O)C2=C(NCC(CNC3=CC=C(C=C3)C(=O)NC(CCC([O-])=O)C([O-])=O)N2)N1',
        'udp': 'OC1C(O)C(OC1COP([O-])(=O)OP([O-])([O-])=O)N1C=CC(=O)NC1=O',
        'utp': 'OC1C(O)C(OC1COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)N1C=CC(=O)NC1=O',
        'x5p': 'OCC(=O)C(O)C(O)COP([O-])([O-])=O',
        'no_structure': '',
    }

    def test__run(self):
        m = AttrDict()
        for name, structure in self.molecules.items():
            if structure:
                m[name] = molecule_util.Molecule(structure=structure).to_mol()
            else:
                m[name] = None

        # 1 --> 1: dr1p <==> d5rp
        results = ezyme.Ezyme()._run([m.dr1p], [m.dr5p])
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].__dict__, {'ec_number': '5.4.2', 'score': 16.0})

        # 1 --> 2: f1p ==> glyceraldehyde + t3p2
        self.assertEqual(ezyme.Ezyme()._run([m.f1p], [m.glyceraldehyde, m.t3p2])[0].ec_number, "2.2.1")

        # 2 --> 1 and example with stoichiometric coefficient != 1: H2O + LeuLeu ==> (2) LEU
        self.assertEqual(ezyme.Ezyme()._run([m.leuleu, m.h2o], [m.leu])[0].ec_number, "3.5.1")

        # 2 --> 2: t3p1 + s7p <==> x5p + r5p
        self.assertEqual(ezyme.Ezyme()._run([m.t3p1, m.s7p], [m.x5p, m.r5p])[0].ec_number, "4.1.2")

        # 3 --> 2: UDP + H + PEP ==> UTP + PYR
        self.assertEqual(ezyme.Ezyme()._run([m.utp, m.pep, m.h], [m.utp, m.pyr])[0].ec_number, '2.7.3')

        # 3 --> 2: metthf + gly + h2o <==> thf + ser
        self.assertEqual(ezyme.Ezyme()._run([m.metthf, m.gly, m.h2o], [m.thf, m.ser])[0].ec_number, '2.1.2')

        # 3 --> 4: IleIle[e] + ATP[c] + H2O[c] ==> IleIle[c] + PI[c] + H[c] + ADP[c]
        self.assertEqual(ezyme.Ezyme()._run([m.ileile, m.atp, m.h2o], [m.ileile, m.pi, m.adp])[0].ec_number, "3.6.1")

        # 4 --> 3: ATP + MET + H2O ==> AMET + PPI + PI + H
        self.assertEqual(ezyme.Ezyme()._run([m.met, m.atp, m.h2o], [m.amet, m.ppi, m.pi, m.h]), [])

        """ test that Ezyme predicts the same EC number for a reversible reaction when the substrates and products and swapped """
        self.assertEqual(ezyme.Ezyme()._run([m.dr5p], [m.dr1p])[0].ec_number, '5.4.2')

        """ test that Ezyme returns None for partiicpants with no defined structure """
        self.assertEqual(ezyme.Ezyme()._run([m.dr5p], [m.no_structure]), None)

        # todo: additional tests:
        #
        #  * cases where the reaction is 3 and 3, but simply looping through each three with sub_inchi_smiles = sub_inchi_smiles[-1:] + sub_inchi_smiles[:-1]
        #    because the position of the 3 is wrong, and I would have to change the position of the middle, so the side, or soemthing

    def test_run(self):
        rxn = data_model.Reaction(participants=[
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['atp'], id='atp'),
                compartment=data_model.Compartment(id='c'),
                coefficient=-1,
                order=0),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['h2o'], id='h2o'),
                compartment=data_model.Compartment(id='c'),
                coefficient=-1,
                order=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['adp'], id='adp'),
                compartment=data_model.Compartment(id='c'),
                coefficient=1,
                order=2),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['ppi'], id='ppi'),
                compartment=data_model.Compartment(id='c'),
                coefficient=1,
                order=3),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['h'], id='h'),
                compartment=data_model.Compartment(id='c'),
                coefficient=1,
                order=4),
        ])
        self.assertEqual(rxn.get_reactants()[0].specie.id, 'atp')
        self.assertEqual(rxn.get_products()[0].specie.id, 'adp')
        self.assertEqual(rxn.get_reactant_product_pairs()[0][0].specie.id, 'atp')
        self.assertEqual(rxn.get_reactant_product_pairs()[0][1].specie.id, 'adp')
        self.assertEqual(rxn.get_reactant_product_pairs()[1][0].specie.id, 'h2o')
        self.assertEqual(rxn.get_reactant_product_pairs()[1][1].specie.id, 'ppi')
        result = ezyme.Ezyme().run(rxn)
        self.assertEqual(result[0].ec_number, '3.6.1')  # true EC is 3.6.1.3

        # example where Ezyme predicts no EC number when the order is swapped
        rxn = data_model.Reaction(participants=[
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['t3p1'], id='atp'),
                coefficient=-1,
                order=0),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['s7p'], id='h2o'),
                coefficient=-1,
                order=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['x5p'], id='adp'),
                coefficient=1,
                order=2),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['r5p'], id='ppi'),
                coefficient=1,
                order=3),
        ])
        result = ezyme.Ezyme().run(rxn)
        self.assertEqual(result[0].ec_number, "4.1.2")  # true EC is 2.2.1.1
        self.assertEqual(result[1].ec_number, "2.2.1")

        rxn = data_model.Reaction(participants=[
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['s7p'], id='atp'),
                coefficient=-1,
                order=0),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['t3p1'], id='h2o'),
                coefficient=-1,
                order=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['x5p'], id='adp'),
                coefficient=1,
                order=2),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['r5p'], id='ppi'),
                coefficient=1,
                order=3),
        ])
        result = ezyme.Ezyme().run(rxn)
        self.assertEqual(result, [])

        # example where Ezyme predicts different EC numbers when the order is swapped
        rxn = data_model.Reaction(participants=[
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['e1dGMP'], id='atp'),
                coefficient=-1,
                order=0),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['h2o'], id='h2o'),
                coefficient=-1,
                order=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['e1dG'], id='adp'),
                coefficient=1,
                order=2),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['pi'], id='ppi'),
                coefficient=1,
                order=3),
        ])
        result = ezyme.Ezyme().run(rxn)
        self.assertEqual(result[0].ec_number, "3.1.3")  # true EC is 3.1.3.89

        rxn = data_model.Reaction(participants=[
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['h2o'], id='atp'),
                coefficient=-1,
                order=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['e1dGMP'], id='h2o'),
                coefficient=-1,
                order=0),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['e1dG'], id='adp'),
                coefficient=1,
                order=2),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['pi'], id='ppi'),
                coefficient=1,
                order=3),
        ])
        result = ezyme.Ezyme().run(rxn)
        self.assertEqual(result[0].ec_number, "3.1.3")  # true EC is 3.1.3.89

        # example where a structure is not defined
        rxn = data_model.Reaction(participants=[
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['h2o'], id='atp'),
                coefficient=-1,
                order=1),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['e1dGMP'], id='h2o'),
                coefficient=-1,
                order=0),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure=self.molecules['e1dG'], id='adp'),
                coefficient=1,
                order=2),
            data_model.ReactionParticipant(
                specie=data_model.Specie(structure='', id='ppi'),
                coefficient=1,
                order=3),
        ])
        result = ezyme.Ezyme().run(rxn)
        self.assertEqual(result, None)
