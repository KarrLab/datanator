""" Tests of the reaction utilities

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-04-12
:Copyright: 2017, Karr Lab
:License: MIT
"""

from attrdict import AttrDict
from kinetic_datanator.util import compartment_util
from kinetic_datanator.util import molecule_util
from kinetic_datanator.util import reaction_util
import unittest


class TestReaction(unittest.TestCase):
    adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    atp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
    h = '[H+]'
    h2o = 'O'
    pi = 'OP([O-])([O-])=O'

    def make_reaction(self):
        return reaction_util.Reaction([
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.h2o, name='h2o'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.adp, name='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.pi, name='pi'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.h, name='h'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
        ])

    def test_init(self):
        rxn = self.make_reaction()
        self.assertEqual(len(rxn.participants), 5)

    def test_name(self):
        rxn = reaction_util.Reaction([], name='rxn')
        self.assertEqual(rxn.name, 'rxn')

    def test_normalize(self):
        rxn = reaction_util.Reaction([
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=0),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.atp, name='atp'),
                compartment=compartment_util.Compartment('e'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.adp, name='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.adp, name='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=3),
        ])

        rxn.normalize()

        self.assertEqual(len(rxn.participants), 4)

        part = list(filter(lambda p: p.molecule.name == 'atp' and p.compartment.name == 'c' and p.coefficient < 0, rxn.participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, -2)

        part = list(filter(lambda p: p.molecule.name == 'atp' and p.compartment.name == 'c' and p.coefficient > 0, rxn.participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, 1)

        part = list(filter(lambda p: p.molecule.name == 'atp' and p.compartment.name == 'e', rxn.participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, -1)

        part = list(filter(lambda p: p.molecule.name == 'adp', rxn.participants))
        self.assertEqual(len(part), 1)
        self.assertEqual(part[0].coefficient, 4)

    def test_get_reactants(self):
        rxn = self.make_reaction()
        self.assertEqual(rxn.get_reactants(), rxn.participants[0:2])

    def test_get_products(self):
        rxn = self.make_reaction()
        self.assertEqual(rxn.get_products(), rxn.participants[2:])

    def test_get_reactant_product_pairs(self):
        rxn = self.make_reaction()
        pairs = rxn.get_reactant_product_pairs()

        self.assertEqual(pairs[0][0].molecule.name, 'atp')
        self.assertEqual(pairs[0][1].molecule.name, 'adp')

        self.assertEqual(pairs[1][0].molecule.name, 'h2o')
        self.assertIn(pairs[1][1].molecule.name, ['h', 'pi'])

        self.assertEqual(pairs[2][0], None)
        self.assertIn(pairs[2][1].molecule.name, ['h', 'pi'])


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
    }

    def test__run(self):
        m = AttrDict({name: molecule_util.Molecule(structure).to_mol() for name, structure in self.molecules.items()})

        # 1 --> 1: dr1p <==> d5rp
        results = reaction_util.Ezyme._run([m.dr1p], [m.dr5p])
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].__dict__, {'ec_number': '5.4.2', 'score': 16.0})

        # 1 --> 2: f1p ==> glyceraldehyde + t3p2
        self.assertEqual(reaction_util.Ezyme._run([m.f1p], [m.glyceraldehyde, m.t3p2])[0].ec_number, "2.2.1")

        # 2 --> 1 and example with stoichiometric coefficient != 1: H2O + LeuLeu ==> (2) LEU
        self.assertEqual(reaction_util.Ezyme._run([m.leuleu, m.h2o], [m.leu])[0].ec_number, "3.5.1")

        # 2 --> 2: t3p1 + s7p <==> x5p + r5p
        self.assertEqual(reaction_util.Ezyme._run([m.t3p1, m.s7p], [m.x5p, m.r5p])[0].ec_number, "4.1.2")

        # 3 --> 2: UDP + H + PEP ==> UTP + PYR
        self.assertEqual(reaction_util.Ezyme._run([m.utp, m.pep, m.h], [m.utp, m.pyr])[0].ec_number, '2.7.3')

        # 3 --> 2: metthf + gly + h2o <==> thf + ser
        self.assertEqual(reaction_util.Ezyme._run([m.metthf, m.gly, m.h2o], [m.thf, m.ser])[0].ec_number, '2.1.2')

        # 3 --> 4: IleIle[e] + ATP[c] + H2O[c] ==> IleIle[c] + PI[c] + H[c] + ADP[c]
        self.assertEqual(reaction_util.Ezyme._run([m.ileile, m.atp, m.h2o], [m.ileile, m.pi, m.adp])[0].ec_number, "3.6.1")

        # 4 --> 3: ATP + MET + H2O ==> AMET + PPI + PI + H
        self.assertEqual(reaction_util.Ezyme._run([m.met, m.atp, m.h2o], [m.amet, m.ppi, m.pi, m.h]), [])

        """ test that Ezyme predicts the same EC number for a reversible reaction when the substrates and products and swapped """
        self.assertEqual(reaction_util.Ezyme._run([m.dr5p], [m.dr1p])[0].ec_number, '5.4.2')

        # todo: additional tests:
        #
        #  * cases where the reaction is 3 and 3, but simply looping through each three with sub_inchi_smiles = sub_inchi_smiles[-1:] + sub_inchi_smiles[:-1]
        #    because the position of the 3 is wrong, and I would have to change the position of the middle, so the side, or soemthing

    def test_run(self):
        rxn = reaction_util.Reaction([
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['atp'], name='atp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['h2o'], name='h2o'),
                compartment=compartment_util.Compartment('c'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['adp'], name='adp'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['ppi'], name='ppi'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['h'], name='h'),
                compartment=compartment_util.Compartment('c'),
                coefficient=1),
        ])
        self.assertEqual(rxn.get_reactants()[0].molecule.name, 'atp')
        self.assertEqual(rxn.get_products()[0].molecule.name, 'adp')
        self.assertEqual(rxn.get_reactant_product_pairs()[0][0].molecule.name, 'atp')
        self.assertEqual(rxn.get_reactant_product_pairs()[0][1].molecule.name, 'adp')
        self.assertEqual(rxn.get_reactant_product_pairs()[1][0].molecule.name, 'h2o')
        self.assertEqual(rxn.get_reactant_product_pairs()[1][1].molecule.name, 'ppi')
        result = reaction_util.Ezyme.run(rxn)
        self.assertEqual(result[0].ec_number, '3.6.1')

        rxn.participants.reverse()
        self.assertEqual(rxn.get_reactants()[0].molecule.name, 'h2o')
        self.assertEqual(rxn.get_products()[0].molecule.name, 'h')
        self.assertEqual(rxn.get_reactant_product_pairs()[0][0].molecule.name, 'atp')
        self.assertEqual(rxn.get_reactant_product_pairs()[0][1].molecule.name, 'adp')
        self.assertEqual(rxn.get_reactant_product_pairs()[1][0].molecule.name, 'h2o')
        self.assertEqual(rxn.get_reactant_product_pairs()[1][1].molecule.name, 'ppi')
        result = reaction_util.Ezyme.run(rxn)
        self.assertEqual(result[0].ec_number, '3.6.1')

        # example where Ezyme predicts no EC number when the order is swapped
        rxn = reaction_util.Reaction([
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['t3p1'], name='atp'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['s7p'], name='h2o'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['x5p'], name='adp'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['r5p'], name='ppi'),
                coefficient=1),
        ])
        result = reaction_util.Ezyme.run(rxn)
        self.assertEqual(result[0].ec_number, "4.1.2")

        rxn = reaction_util.Reaction([
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['s7p'], name='atp'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['t3p1'], name='h2o'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['x5p'], name='adp'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['r5p'], name='ppi'),
                coefficient=1),
        ])
        result = reaction_util.Ezyme.run(rxn)
        self.assertEqual(result[0].ec_number, "4.1.2")

        # example where Ezyme predicts different EC numbers when the order is swapped
        rxn = reaction_util.Reaction([
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['e1dGMP'], name='atp'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['h2o'], name='h2o'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['e1dG'], name='adp'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['pi'], name='ppi'),
                coefficient=1),
        ])
        result = reaction_util.Ezyme.run(rxn)
        self.assertEqual(result[0].ec_number, "3.1.3")

        rxn = reaction_util.Reaction([
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['h2o'], name='atp'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['e1dGMP'], name='h2o'),
                coefficient=-1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['e1dG'], name='adp'),
                coefficient=1),
            reaction_util.ReactionParticipant(
                molecule=molecule_util.Molecule(self.molecules['pi'], name='ppi'),
                coefficient=1),
        ])
        result = reaction_util.Ezyme.run(rxn)
        self.assertEqual(result[0].ec_number, "3.1.3")
