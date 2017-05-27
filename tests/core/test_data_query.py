""" Test of observation filtering and ordering

:Author: Jonathan Karr <jonrkarr@gmail.com>
:Author: Yosef Roth <yosefdroth@gmail.com>
:Date: 2017-04-10
:Copyright: 2017, Karr Lab
:License: MIT
"""

from kinetic_datanator.core import data_model
from kinetic_datanator.core import data_query
from kinetic_datanator.util import warning_util
import copy
import math
import numpy
import scipy.stats
import unittest


warning_util.disable_warnings()


class TestDataQueryGenerator(unittest.TestCase):

    class ConcreteDataQueryGenerator(data_query.DataQueryGenerator):

        def get_observed_values(self):
            pass

    def test_filter_observed_values(self):
        gen = self.ConcreteDataQueryGenerator()
        gen.filters = [
            data_query.TemperatureRangeFilter(min=36., max=38.),
        ]
        observed_values = [
            data_model.ObservedValue(observation=data_model.Observation(environment=data_model.Environment(temperature=36.5))),
            data_model.ObservedValue(observation=data_model.Observation(environment=data_model.Environment(temperature=35.0))),
            data_model.ObservedValue(observation=data_model.Observation(environment=data_model.Environment(temperature=37.0))),
        ]
        result = gen.filter_observed_values(None, observed_values)
        self.assertEqual(set(result.observed_values), set([observed_values[0], observed_values[2]]))

    @unittest.skip('implement me')
    def test_get_consensus(self):
        gen = self.ConcreteDataQueryGenerator()
        # gen.get_consensus()


class TestFilters(unittest.TestCase):

    def test_TaxonomicDistanceFilter(self):
        # example 1
        f = data_query.TaxonomicDistanceFilter('Mycoplasma pneumoniae M129')
        self.assertEqual(f.max, 8.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma pneumoniae M129')))
        self.assertEqual(f.score(None, ov), 1.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma pneumoniae')))
        self.assertEqual(f.scale, 8./5)
        self.assertEqual(f.score(None, ov), math.exp(-1/f.scale))

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Escherichia coli')))
        self.assertEqual(f.score(None, ov), math.exp(-5))

        # example 2
        f = data_query.TaxonomicDistanceFilter('Mycoplasma pneumoniae M129', max=5)
        self.assertEqual(f.max, 5.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma pneumoniae M129')))
        self.assertEqual(f.score(None, ov), 1.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma pneumoniae')))
        self.assertEqual(f.score(None, ov), math.exp(-1/f.scale))

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Escherichia coli')))
        self.assertEqual(f.score(None, ov), -1)

        # example 3
        f = data_query.TaxonomicDistanceFilter('Mycoplasma genitalium')
        self.assertEqual(f.max, 7.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma genitalium G37')))
        self.assertEqual(f.score(None, ov), 1.)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma genitalium')))
        self.assertEqual(f.score(None, ov), 1)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(taxon='Mycoplasma')))
        self.assertEqual(f.score(None, ov), math.exp(-1/f.scale))

    def test_OptionsFilter(self):
        f = data_query.OptionsFilter(('observation', 'genetics', 'variation', ), [''])
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='')))
        self.assertEqual(f.score(None, ov), 1)
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='wildtype')))
        self.assertEqual(f.score(None, ov), -1)
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='Delta gene-01')))
        self.assertEqual(f.score(None, ov), -1)

        f = data_query.OptionsFilter(('observation', 'genetics', 'variation', ), ['wildtype'])
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='')))
        self.assertEqual(f.score(None, ov), -1)
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='wildtype')))
        self.assertEqual(f.score(None, ov), 1)
        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='Delta gene-01')))
        self.assertEqual(f.score(None, ov), -1)

    def test_WildtypeFilter(self):
        f = data_query.WildtypeFilter()

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='')))
        self.assertEqual(f.score(None, ov), 1)

        ov = data_model.ObservedValue(
            observation=data_model.Observation(genetics=data_model.Genetics(variation='Delta gene-01')))
        self.assertEqual(f.score(None, ov), -1)

    def test_SpecieStructuralSimilarityFilter(self):
        adp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
        atp = 'NC1=C2N=CN(C3OC(COP([O-])(=O)OP([O-])(=O)OP([O-])([O-])=O)C(O)C3O)C2=NC=N1'
        h2o = 'O'

        # min_similarity = 0
        f = data_query.SpecieStructuralSimilarityFilter(min_similarity=0.)

        ov = data_model.ObservedValue(observable=data_model.Observable(specie=data_model.Specie(structure=adp)))
        numpy.testing.assert_almost_equal(f.score(data_model.Specie(structure=atp), ov), 0.955, decimal=3)

        ov = data_model.ObservedValue(observable=data_model.Observable(specie=data_model.Specie(structure=h2o)))
        numpy.testing.assert_almost_equal(f.score(data_model.Specie(structure=atp), ov), 0, decimal=3)

        # min_similarity = 0.75
        f = data_query.SpecieStructuralSimilarityFilter(min_similarity=0.75)

        ov = data_model.ObservedValue(observable=data_model.Observable(specie=data_model.Specie(structure=adp)))
        numpy.testing.assert_almost_equal(f.score(data_model.Specie(structure=atp), ov), 0.955, decimal=3)

        ov = data_model.ObservedValue(observable=data_model.Observable(specie=data_model.Specie(structure=h2o)))
        numpy.testing.assert_almost_equal(f.score(data_model.Specie(structure=atp), ov), -1., decimal=3)

    def test_SpecieSequenceSimilarityFilter(self):
        seq1 = 'CTAACTCTACCTCGTATGTATGGAAGTTCGTCTATCTCTGGTCGGTTGCT'
        seq2 = 'CTAACTCTACCTCGTATTATGGAAGTTCGTCTATCTTCTGGTCGGTTGCT'

        ov = data_model.ObservedValue(observable=data_model.Observable(specie=data_model.PolymerSpecie(sequence=seq1)))

        # min_similarity = 0
        f = data_query.SpecieSequenceSimilarityFilter(min_similarity=0.)
        self.assertEqual(f.score(data_model.PolymerSpecie(sequence=seq2), ov), 0.96)

        # min_similarity = 0.75
        f = data_query.SpecieSequenceSimilarityFilter(min_similarity=0.98)
        self.assertEqual(f.score(data_model.PolymerSpecie(sequence=seq2), ov), -1)

    def test_ReactionSimilarityFilter(self):
        atp_structure = (
            'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20'
            '/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-4/t4-,6-,7-,10-/m1/s1'
        )
        h2o_structure = 'InChI=1S/H2O/h1H2'
        adp_structure = (
            'InChI=1S/C10H15N5O10P2/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(24-10)1-23-27(21,22)25-26(18,19)20'
            '/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H2,11,12,13)(H2,18,19,20)/p-3/t4-,6-,7-,10-/m1/s1'
        )
        pi_structure = 'InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-2'
        h_structure = 'InChI=1S/p+1/i/hH'
        def get_reaction(pi_structure=pi_structure, ec='1.1.1.1'):
            return data_model.Reaction(
                participants=[
                    data_model.ReactionParticipant(coefficient=-1, specie=data_model.Specie(structure=atp_structure)),
                    data_model.ReactionParticipant(coefficient=-1, specie=data_model.Specie(structure=h2o_structure)),
                    data_model.ReactionParticipant(coefficient=1, specie=data_model.Specie(structure=adp_structure)),
                    data_model.ReactionParticipant(coefficient=1, specie=data_model.Specie(structure=pi_structure)),
                    data_model.ReactionParticipant(coefficient=1, specie=data_model.Specie(structure=h_structure)),
                ],
                cross_references=[
                    data_model.Resource(namespace='ec-code', id=ec)
                ])

        rxn = get_reaction()
        f = data_query.ReactionSimilarityFilter(min_ec_level=3, scale=1)

        # same participants
        ov = data_model.ObservedValue(observable=data_model.Observable(interaction=get_reaction()))
        numpy.testing.assert_almost_equal(f.score(rxn, ov), 1, decimal=3)

        # similiar participants
        ov = data_model.ObservedValue(observable=data_model.Observable(interaction=get_reaction(pi_structure='InChI=1S/H3O4P/c1-5(2,3)4')))
        numpy.testing.assert_almost_equal(f.score(rxn, ov), 1, decimal=3)

        # different participants, same 4-digit EC
        ov = data_model.ObservedValue(observable=data_model.Observable(interaction=get_reaction(pi_structure='InChI=1S/H4O4P/c1-5(2,3)4')))
        numpy.testing.assert_almost_equal(f.score(rxn, ov), math.exp(-1), decimal=3)

        # different participants, same 3-digit EC
        ov = data_model.ObservedValue(observable=data_model.Observable(
            interaction=get_reaction(pi_structure='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1.2')))
        numpy.testing.assert_almost_equal(f.score(rxn, ov), math.exp(-2), decimal=3)

        ov = data_model.ObservedValue(observable=data_model.Observable(
            interaction=get_reaction(pi_structure='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1')))
        numpy.testing.assert_almost_equal(f.score(rxn, ov), math.exp(-2), decimal=3)

        ov = data_model.ObservedValue(observable=data_model.Observable(
            interaction=get_reaction(pi_structure='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1.')))
        numpy.testing.assert_almost_equal(f.score(rxn, ov), math.exp(-2), decimal=3)

        ov = data_model.ObservedValue(observable=data_model.Observable(
            interaction=get_reaction(pi_structure='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1.-')))
        numpy.testing.assert_almost_equal(f.score(rxn, ov), math.exp(-2), decimal=3)

        # different participants, same 2-digit EC
        ov = data_model.ObservedValue(observable=data_model.Observable(
            interaction=get_reaction(pi_structure='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.2.1')))
        numpy.testing.assert_almost_equal(f.score(rxn, ov), -1, decimal=3)

        # target reaction only has 3 digits
        rxn_1_1_1 = get_reaction(ec='1.1.1')
        rxn_1_1_1_1 = get_reaction(ec='1.1.1.1')
        f1 = data_query.ReactionSimilarityFilter(min_ec_level=3, scale=1)
        f2 = data_query.ReactionSimilarityFilter(min_ec_level=3, scale=1)

        ov = data_model.ObservedValue(observable=data_model.Observable(
            interaction=get_reaction(pi_structure='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1')))
        numpy.testing.assert_almost_equal(f1.score(rxn_1_1_1, ov), math.exp(-2), decimal=3)
        numpy.testing.assert_almost_equal(f2.score(rxn_1_1_1_1, ov), math.exp(-2), decimal=3)

        ov = data_model.ObservedValue(observable=data_model.Observable(
            interaction=get_reaction(pi_structure='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.1.1')))
        numpy.testing.assert_almost_equal(f1.score(rxn_1_1_1, ov), math.exp(-2), decimal=3)
        numpy.testing.assert_almost_equal(f2.score(rxn_1_1_1_1, ov), math.exp(-1), decimal=3)

        ov = data_model.ObservedValue(observable=data_model.Observable(
            interaction=get_reaction(pi_structure='InChI=1S/H4O4P/c1-5(2,3)4', ec='1.1.2.1')))
        numpy.testing.assert_almost_equal(f1.score(rxn_1_1_1, ov), -1, decimal=3)
        numpy.testing.assert_almost_equal(f2.score(rxn_1_1_1_1, ov), -1, decimal=3)

        # reverse direction, different numbers of reactants/products
        f = data_query.ReactionSimilarityFilter(min_ec_level=3, scale=1)

        for_rxn = get_reaction()
        rev_rxn = get_reaction()
        for part in rev_rxn.participants:
            part.coefficient = -1 * part.coefficient

        ov = data_model.ObservedValue(observable=data_model.Observable(interaction=for_rxn))
        self.assertEqual(f.score(for_rxn, ov), 1.)

        ov = data_model.ObservedValue(observable=data_model.Observable(interaction=rev_rxn))
        self.assertEqual(f.score(for_rxn, ov), -1.)

        # reverse direction, same numbers of reactants/products
        for_rxn = data_model.Reaction(
            participants=[
                data_model.ReactionParticipant(coefficient=-1, specie=data_model.Specie(structure=atp_structure)),
                data_model.ReactionParticipant(coefficient=1, specie=data_model.Specie(structure=adp_structure)),
            ],
            cross_references=[
                data_model.Resource(namespace='ec-code', id='1.1.1.1')
            ])
        rev_rxn = data_model.Reaction(
            participants=[
                data_model.ReactionParticipant(coefficient=1, specie=data_model.Specie(structure=atp_structure)),
                data_model.ReactionParticipant(coefficient=-1, specie=data_model.Specie(structure=adp_structure)),
            ],
            cross_references=[
                data_model.Resource(namespace='ec-code', id='1.1.1.1')
            ])

        ov = data_model.ObservedValue(observable=data_model.Observable(interaction=for_rxn))
        self.assertEqual(f.score(for_rxn, ov), 1.)

        ov = data_model.ObservedValue(observable=data_model.Observable(interaction=rev_rxn))
        self.assertEqual(f.score(for_rxn, ov), -1.)

    def test_ReactionParticipantFilter(self):
        atp = data_model.Specie(structure=(
            'InChI=1S/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20'
            '/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/p-4/t4-,6-,7-,10-/m1/s1'
        ))
        h2o = data_model.Specie(structure='InChI=1S/H2O/h1H2')
        adp = data_model.Specie(structure=(
            'InChI=1S/C10H15N5O10P2/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(24-10)1-23-27(21,22)25-26(18,19)20'
            '/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H2,11,12,13)(H2,18,19,20)/p-3/t4-,6-,7-,10-/m1/s1'
        ))
        pi = data_model.Specie(structure='InChI=1S/H3O4P/c1-5(2,3)4/h(H3,1,2,3,4)/p-2')
        h = data_model.Specie(structure='InChI=1S/p+1/i/hH')

        glc = data_model.Specie(structure='InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1')
        glc_2 = data_model.Specie(structure='InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2')
        gtp = data_model.Specie(structure=(
            'InChI=1S/C10H16N5O14P3/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(27-9)1-26-31(22,23)29-32(24,25)28-30(19,20)21'
            '/h2-3,5-6,9,16-17H,1H2,(H,22,23)(H,24,25)(H2,19,20,21)(H3,11,13,14,18)/p-4/t3-,5-,6-,9-/m1/s1'
        ))
        lactate = data_model.Specie(structure='InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/p-1')

        def get_reaction(ntp=atp, glc=glc):
            return data_model.Reaction(
                participants=[
                    data_model.ReactionParticipant(coefficient=-1, specie=ntp),
                    data_model.ReactionParticipant(coefficient=-1, specie=h2o),
                    data_model.ReactionParticipant(coefficient=1, specie=adp),
                    data_model.ReactionParticipant(coefficient=1, specie=pi),
                    data_model.ReactionParticipant(coefficient=1, specie=h),
                    data_model.ReactionParticipant(coefficient=0, specie=glc),
                ])

        f = data_query.ReactionParticipantFilter()

        rxn_atp = get_reaction(ntp=atp)
        rxn_gtp = get_reaction(ntp=gtp)
        rxn_lac = get_reaction(ntp=lactate)
        rxn_glc_2 = get_reaction(glc=glc_2)

        # identical reactant
        ov = data_model.ObservedValue(observable=data_model.Observable(property='Km', specie=atp, interaction=rxn_atp))
        self.assertEqual(f.score(rxn_atp, ov), 1)

        # identical product
        ov = data_model.ObservedValue(observable=data_model.Observable(property='Ki', specie=adp, interaction=rxn_atp))
        self.assertEqual(f.score(rxn_atp, ov), 1)

        # similar modifier
        ov = data_model.ObservedValue(observable=data_model.Observable(property='Ki', specie=glc_2, interaction=rxn_glc_2))
        self.assertEqual(f.score(rxn_atp, ov), 1)

        # similar species
        ov = data_model.ObservedValue(observable=data_model.Observable(property='Km', specie=gtp, interaction=rxn_gtp))
        numpy.testing.assert_almost_equal(f.score(rxn_atp, ov), 0.767, decimal=3)

        # different species
        ov = data_model.ObservedValue(observable=data_model.Observable(property='Km', specie=lactate, interaction=rxn_lac))
        self.assertEqual(f.score(rxn_atp, ov), -1)

        # property without species
        ov = data_model.ObservedValue(observable=data_model.Observable(property='kcat', interaction=rxn_atp))
        self.assertEqual(f.score(rxn_atp, ov), 1)

    def test_RangeFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = data_query.RangeFilter(('observation', 'environment', 'temperature', ), min=15., max=30.)
        self.assertEqual(f.score(None, obs_val(float('nan'))), -1)
        self.assertEqual(f.score(None, obs_val(10)), -1)
        self.assertEqual(f.score(None, obs_val(15.0)), 1)
        self.assertEqual(f.score(None, obs_val(15.01)), 1)
        self.assertEqual(f.score(None, obs_val(29.99)), 1)
        self.assertEqual(f.score(None, obs_val(30.0)), 1)
        self.assertEqual(f.score(None, obs_val(31)), -1)

        f = data_query.RangeFilter(('observation', 'environment', 'temperature', ), min=float('nan'), max=30.)
        self.assertEqual(f.score(None, obs_val(float('nan'))), -1)
        self.assertEqual(f.score(None, obs_val(10)), 1)
        self.assertEqual(f.score(None, obs_val(15.0)), 1)
        self.assertEqual(f.score(None, obs_val(15.01)), 1)
        self.assertEqual(f.score(None, obs_val(29.99)), 1)
        self.assertEqual(f.score(None, obs_val(30.0)), 1)
        self.assertEqual(f.score(None, obs_val(31)), -1)

        f = data_query.RangeFilter(('observation', 'environment', 'temperature', ), min=15., max=float('nan'))
        self.assertEqual(f.score(None, obs_val(float('nan'))), -1)
        self.assertEqual(f.score(None, obs_val(10)), -1)
        self.assertEqual(f.score(None, obs_val(15.0)), 1)
        self.assertEqual(f.score(None, obs_val(15.01)), 1)
        self.assertEqual(f.score(None, obs_val(29.99)), 1)
        self.assertEqual(f.score(None, obs_val(30.0)), 1)
        self.assertEqual(f.score(None, obs_val(31)), 1)

        f = data_query.RangeFilter(('observation', 'environment', 'temperature', ))
        self.assertEqual(f.score(None, obs_val(float('nan'))), 1)
        self.assertEqual(f.score(None, obs_val(10)), 1)
        self.assertEqual(f.score(None, obs_val(15.0)), 1)
        self.assertEqual(f.score(None, obs_val(15.01)), 1)
        self.assertEqual(f.score(None, obs_val(29.99)), 1)
        self.assertEqual(f.score(None, obs_val(30.0)), 1)
        self.assertEqual(f.score(None, obs_val(31)), 1)

    def test_TemperatureRangeFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = data_query.TemperatureRangeFilter(min=15., max=30.)
        self.assertEqual(f.score(None, obs_val(10)), -1)
        self.assertEqual(f.score(None, obs_val(15.0)), 1)
        self.assertEqual(f.score(None, obs_val(15.01)), 1)
        self.assertEqual(f.score(None, obs_val(29.99)), 1)
        self.assertEqual(f.score(None, obs_val(30.0)), 1)
        self.assertEqual(f.score(None, obs_val(31)), -1)

    def test_PhRangeFilter(self):
        def obs_val(ph):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(ph=ph)))

        f = data_query.PhRangeFilter(min=5., max=9.)
        self.assertEqual(f.score(None, obs_val(3)), -1)
        self.assertEqual(f.score(None, obs_val(5.0)), 1)
        self.assertEqual(f.score(None, obs_val(5.01)), 1)
        self.assertEqual(f.score(None, obs_val(8.99)), 1)
        self.assertEqual(f.score(None, obs_val(9.0)), 1)
        self.assertEqual(f.score(None, obs_val(10)), -1)

    def test_NormalFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = data_query.NormalFilter(('observation', 'environment', 'temperature', ), mean=37, std=1)
        self.assertEqual(f.score(None, obs_val(37)), 1)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(36)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(38)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(35)), 2 * scipy.stats.norm.cdf(-2), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(39)), 2 * scipy.stats.norm.cdf(-2), decimal=5)

    def test_TemperatureNormalFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = data_query.TemperatureNormalFilter(mean=37, std=1)
        self.assertEqual(f.score(None, obs_val(37)), 1)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(36)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(38)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(35)), 2 * scipy.stats.norm.cdf(-2), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(39)), 2 * scipy.stats.norm.cdf(-2), decimal=5)

    def test_PhNormalFilter(self):
        def obs_val(ph):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(ph=ph)))

        f = data_query.PhNormalFilter(mean=7, std=1)
        self.assertEqual(f.score(None, obs_val(7)), 1)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(6)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(8)), 2 * scipy.stats.norm.cdf(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(5)), 2 * scipy.stats.norm.cdf(-2), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(9)), 2 * scipy.stats.norm.cdf(-2), decimal=5)

    def test_ExponentialFilter(self):
        def obs_val(temperature):
            return data_model.ObservedValue(
                observation=data_model.Observation(environment=data_model.Environment(temperature=temperature)))

        f = data_query.ExponentialFilter(('observation', 'environment', 'temperature', ), center=1., scale=1.)
        self.assertEqual(f.score(None, obs_val(1)), 1)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(100)), 0, decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(2)), math.exp(-1), decimal=5)
        numpy.testing.assert_almost_equal(f.score(None, obs_val(3)), math.exp(-2), decimal=5)


class TestFilterResult(unittest.TestCase):

    def test(self):
        ov = [
            data_model.ObservedValue(value=1),
            data_model.ObservedValue(value=2),
        ]
        f = [
            data_query.RangeFilter(('environment', 'temperature', )),
            data_query.RangeFilter(('environment', 'ph', )),
        ]
        s = [0.3, 0.5]

        filter_result = data_query.FilterResult(copy.copy(ov), copy.copy(s), [1, 2], ov, s)
        self.assertEqual(filter_result.observed_values, ov)
        self.assertEqual(filter_result.scores, s)
        self.assertEqual(filter_result.observed_value_indices, [1, 2])
        self.assertEqual(filter_result.all_observed_values, ov)
        self.assertEqual(filter_result.all_scores, s)


class TestFilterRunner(unittest.TestCase):

    def test_score(self):
        ov = [
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
        ]
        f = [
            data_query.TemperatureRangeFilter(min=36.5, max=37.5),
            data_query.PhRangeFilter(min=6.5, max=7.5),
            data_query.TemperatureNormalFilter(mean=37., std=1),
        ]

        runner = data_query.FilterRunner(f[0])
        self.assertEqual(list(runner.score(None, ov).ravel()), [1., 1., -1., -1., 1., 1., -1., -1.])

        runner = data_query.FilterRunner(f[1])
        self.assertEqual(list(runner.score(None, ov).ravel()), [1., -1., 1., -1., 1., -1., 1., -1.])

        runner = data_query.FilterRunner(f[2])
        s = 2 * scipy.stats.norm.cdf(-1)
        numpy.testing.assert_almost_equal(list(runner.score(None, ov).ravel()), [1., 1., s, s, 1., 1., s, s], decimal=5)

    def test_filter(self):
        ov = [
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
        ]

        # one filter
        s = numpy.array([1, -1, 0.5, 1, -1, 0, -1, 0.3], ndmin=2).transpose()
        runner = data_query.FilterRunner([])
        filt_o, filt_s, i_filt_o = runner.filter(ov, s)

        self.assertEqual(filt_o, [ov[0], ov[2], ov[3], ov[5], ov[7]])
        numpy.testing.assert_equal(filt_s, numpy.array([1, 0.5, 1, 0, 0.3], ndmin=2).transpose())
        self.assertEqual(i_filt_o, [0, 2, 3, 5, 7])

        # multiple filters
        s = numpy.array([[1, -1, 0.5, 1, -1, 0, -1, 0.3], [0, 0, 0, 0, 0, 0, 0, -1]], ndmin=2).transpose()
        runner = data_query.FilterRunner([])
        filt_o, filt_s, i_filt_o = runner.filter(ov, s)

        self.assertEqual(filt_o, [ov[0], ov[2], ov[3], ov[5]])
        numpy.testing.assert_equal(filt_s, numpy.array([[1, 0.5, 1, 0], [0, 0, 0, 0]], ndmin=2).transpose())
        self.assertEqual(i_filt_o, [0, 2, 3, 5])

    def test_order(self):
        ov = [
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
        ]

        # one filter
        s = numpy.array([0.9, -1, 0.5, 1, -1.1, 0, -1.3, 0.3], ndmin=2).transpose()
        runner = data_query.FilterRunner([])
        order_o, order_s, i_order_o = runner.order(ov, s)

        self.assertEqual(order_o, [ov[3], ov[0], ov[2], ov[7], ov[5], ov[1], ov[4], ov[6]])
        numpy.testing.assert_equal(order_s, numpy.array([1, 0.9, 0.5, 0.3, 0, -1, -1.1, -1.3], ndmin=2).transpose())
        self.assertEqual(i_order_o, [3, 0, 2, 7, 5, 1, 4, 6])

        # multiple filters
        s = numpy.array([[0.9, -1, 0.5, 1, -1.1, 0, -1.3, 0.3], [0, 0, 0, 0, 0, 0, 0, 0]], ndmin=2).transpose()
        runner = data_query.FilterRunner([])
        order_o, order_s, i_order_o = runner.order(ov, s)

        self.assertEqual(order_o, [ov[3], ov[0], ov[2], ov[7], ov[5], ov[1], ov[4], ov[6]])
        numpy.testing.assert_equal(order_s, numpy.array(
            [[1, 0.9, 0.5, 0.3, 0, -1, -1.1, -1.3], [0, 0, 0, 0, 0, 0, 0, 0]], ndmin=2).transpose())
        self.assertEqual(i_order_o, [3, 0, 2, 7, 5, 1, 4, 6])

    def test_run(self):
        ov = [
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37, ph=6))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=7))),
            data_model.ObservedValue(value=1, observation=data_model.Observation(
                environment=data_model.Environment(temperature=36, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37.1, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37.1, ph=6))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37.01, ph=7))),
            data_model.ObservedValue(value=2, observation=data_model.Observation(
                environment=data_model.Environment(temperature=37.01, ph=6))),
        ]
        f = [
            data_query.TemperatureNormalFilter(mean=37, std=1),
            data_query.PhRangeFilter(min=6.5, max=7.5),
        ]
        runner = data_query.FilterRunner(f)

        # return_info=True
        result = runner.run(None, ov, return_info=True)
        self.assertEqual(result.all_observed_values, ov)
        s1 = 2 * scipy.stats.norm.cdf(-1)
        s01 = 2 * scipy.stats.norm.cdf(-0.1)
        s001 = 2 * scipy.stats.norm.cdf(-0.01)
        numpy.testing.assert_almost_equal(result.all_scores, numpy.array([[1., 1., s1, s1, s01, s01, s001, s001],
                                                                          [1., -1., 1., -1., 1., -1., 1., -1.]], ndmin=2).transpose())

        self.assertEqual(result.observed_values, [ov[0], ov[6], ov[4], ov[2]])
        numpy.testing.assert_almost_equal(result.scores, numpy.array([[1., s001, s01, s1], [1., 1., 1., 1.]], ndmin=2).transpose())
        self.assertEqual(result.observed_value_indices, [0, 6, 4, 2])

        # return_info=False
        ordered_observed_values = runner.run(None, ov, return_info=False)
        self.assertEqual(ordered_observed_values, [ov[0], ov[6], ov[4], ov[2]])


class TestConsensusGenerator(unittest.TestCase):

    def test_calc_average(self):
        gen = data_query.ConsensusGenerator()

        # mean
        values = [1, 3]
        weights = [1, 10]

        value, error, method = gen.calc_average(values, weights=None, method='mean')
        self.assertEqual(value, 2)
        self.assertEqual(method, data_model.ConsensusMethod.mean)

        value, error, method = gen.calc_average(values, weights=weights, method='mean')
        self.assertEqual(value, (1. * 1. + 3. * 10.) / (1. + 10.))
        self.assertEqual(method, data_model.ConsensusMethod.weighted_mean)

        # median
        values = [1, 3, 2]
        weights = [1, 10, 9]

        value, error, method = gen.calc_average(values, weights=None, method='median')
        self.assertEqual(value, 2)
        self.assertEqual(method, data_model.ConsensusMethod.median)

        value, error, method = gen.calc_average(values, weights=weights, method='median')
        self.assertEqual(value, 2.5)
        self.assertEqual(method, data_model.ConsensusMethod.weighted_median)

        # mode
        values = [1, 1, 3, 2]
        weights = [1, 3, 10, 5]

        value, error, method = gen.calc_average(values, weights=None, method='mode')
        self.assertEqual(value, 1)
        self.assertEqual(method, data_model.ConsensusMethod.mode)

        value, error, method = gen.calc_average(values, weights=weights, method='mode')
        self.assertEqual(value, 3)
        self.assertEqual(method, data_model.ConsensusMethod.weighted_mode)

        # error
        values = [1, 3]

        value, error, method = gen.calc_average(values, weights=None, method='mean')
        self.assertEqual(error, numpy.std(values))

        value, error, method = gen.calc_average(values, weights=[1, 1], method='mean')
        self.assertEqual(error, numpy.std(values))

        value, error, method = gen.calc_average(values, weights=[2, 2], method='mean')
        self.assertEqual(error, numpy.std(values))

        # handling nan values
        value, error, method = gen.calc_average([1, numpy.nan], weights=None, method='mean')
        self.assertEqual(value, 1)
        self.assertTrue(numpy.isnan(error))
        self.assertEqual(method, data_model.ConsensusMethod.mean)

        value, error, method = gen.calc_average([1, 2, numpy.nan], weights=None, method='mean')
        self.assertEqual(value, 1.5)
        self.assertEqual(error, numpy.std([1, 2]))
        self.assertEqual(method, data_model.ConsensusMethod.mean)

        value, error, method = gen.calc_average([numpy.nan, numpy.nan], weights=None, method='mean')
        self.assertTrue(numpy.isnan(value))
        self.assertTrue(numpy.isnan(error))
        self.assertEqual(method, None)

        value, error, method = gen.calc_average([1, 1, 2], weights=[1, 1, numpy.nan], method='mean')
        self.assertEqual(value, numpy.mean([1, 1, 2]))
        self.assertTrue(error, numpy.std([1, 1, 2]))
        self.assertEqual(method, data_model.ConsensusMethod.mean)

        value, error, method = gen.calc_average([1, 1, 2], weights=[1, numpy.nan, numpy.nan], method='mean')
        self.assertEqual(value, numpy.mean([1, 1, 2]))
        self.assertTrue(error, numpy.std([1, 1, 2]))
        self.assertEqual(method, data_model.ConsensusMethod.mean)
