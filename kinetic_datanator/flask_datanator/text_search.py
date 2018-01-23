#-*- coding: utf-8 -*-

""" Test of common schema

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-12-18
:Copyright: 2017, Karr Lab
:License: MIT
"""


from kinetic_datanator.data_query import dna_protein_interactions, metabolite_concentrations, protein_concentrations, protein_protein_interactions, reaction_kinetics
from kinetic_datanator.flask_datanator import flask_common_schema, models
import flask_whooshalchemy
import tempfile


class TextSearchSession(object):
    """
    Represents a text search through the database grabbing relevant information

    """

    def __init__(self):

        flaskdb = flask_common_schema.FlaskCommonSchema()

        self.metabolite_concentration_query = metabolite_concentrations.FlaskMetaboliteConcentrationsQueryGenerator()
        self.protein_concentration_query = protein_concentrations.FlaskProteinConcentrationsQueryGenerator()
        self.reaction_kinetics_query = reaction_kinetics.ReactionKineticsQueryGenerator()
        self.protein_interaction_query = protein_protein_interactions.ProteintoProteinInteractionQueryGenerator()
        self.protein_dna_query = dna_protein_interactions.ProteintoDNAInteractionQueryGenerator()
        self.dna_protein_query = dna_protein_interactions.DNAtoProteinInteractionQueryGenerator()

        for item in flaskdb.text_indicies:
            flask_whooshalchemy.whoosh_index(models.app, item)


    def return_information(self, string):
        """
        Args:
            string (:obj:`str`): string of what is to be queried

        Returns:
            query_models (:obj:`dict` of many :obj:``): Related Metadata ID

        """

        query_models = {}
        db_models = self.collect_data_objects(string)
        query_models['MetaboliteConcentration'] = self.text_return_metabolite_concentration()
        query_models['ProteinConcentration'] = self.text_return_protein_concentration()

        return (db_models, query_models)

    def collect_data_objects(self, string):
        self.ans = {}

        self.ans['Compound'] = models.Compound.query.whoosh_search(string).all()
        self.ans['ProteinComplex'] = models.ProteinComplex.query.whoosh_search(string).all()
        self.ans['ProteinInteractions'] = models.ProteinInteractions.query.whoosh_search(string).all()
        self.ans['Taxon'] = models.Taxon.query.whoosh_search(string).all()
        self.ans['Synonym'] = models.Synonym.query.whoosh_search(string).all()
        self.ans['CellLine'] = models.CellLine.query.whoosh_search(string).all()
        self.ans['CellCompartment'] = models.CellCompartment.query.whoosh_search(string).all()
        self.ans['ProteinSubunit'] = models.ProteinSubunit.query.whoosh_search(string).all()

        return self.ans

    def text_return_metabolite_concentration(self):
        ans = []
        for items in self.ans['Compound']:
            ans.append(
                self.metabolite_concentration_query.get_observed_values(items))

        return ans

    def text_return_protein_concentration(self):
        ans = []
        for items in self.ans['ProteinSubunit']:
            ans.append(
                self.protein_concentration_query.get_observed_values(items))

        return ans

    def text_return_reaction_kinetics(self):
        pass
