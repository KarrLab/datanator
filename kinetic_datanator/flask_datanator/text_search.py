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


    def return_all_queries(self, string):
        """
        Returns all the collected information from the given string
        Args:
            string (:obj:`str`): string of what is to be queried

        Returns:
            query_models (:obj:`dict` of many :obj:``): Related Metadata ID

        """
        db_models = collect_objects(string)

        query_models = {}
        query_models['MetaboliteConcentration'] = self.text_return_metabolite_concentration()
        query_models['ProteinConcentration'] = self.text_return_protein_concentration()
        query_models['ReactionKinetics'] = self.text_return_reaction_kinetics()
        query_models['ProteinInteractions'] = self.text_return_protein_protein_interactions()
        query_models['ProteinDNA'] = self.text_return_protein_dna_interactions()
        query_models['DNAProtein'] = self.text_dna_protein_interactions()

        return (db_models, query_models)

    def collect_objects(self, string):
        """
        Collects and creates dictionary of all Database objects coming up in a full text search for related string

        """


        self.db_models = {}
        self.db_models['Compound'] = models.Compound.query.whoosh_search(string).all()
        self.db_models['ProteinComplex'] = models.ProteinComplex.query.whoosh_search(string).all()
        self.db_models['ProteinInteractions'] = models.ProteinInteractions.query.whoosh_search(string).all()
        self.db_models['Taxon'] = models.Taxon.query.whoosh_search(string).all()
        self.db_models['Synonym'] = models.Synonym.query.whoosh_search(string).all()
        self.db_models['CellLine'] = models.CellLine.query.whoosh_search(string).all()
        self.db_models['CellCompartment'] = models.CellCompartment.query.whoosh_search(string).all()
        self.db_models['ProteinSubunit'] = models.ProteinSubunit.query.whoosh_search(string).all()

        return self.db_models


    def text_return_metabolite_concentration(self):
        ans = []
        for items in self.db_models['Compound']:
            ans.append(self.metabolite_concentration_query.get_observed_values(items))

        return ans

    def text_return_protein_concentration(self):
        ans = []
        for items in self.db_models['ProteinSubunit']:
            ans.append(self.protein_concentration_query.get_observed_values(items))

        return ans

    def text_return_reaction_kinetics(self):
        pass

    def text_return_protein_protein_interactions(self):
        interaction = []
        complex = []
        for items in self.db_models['ProteinSubunit']:
            interact, plex = self.protein_interaction_query.get_observable_interactions_and_complex(items)
            interaction.append(interact)
            complex.append(plex)
        return interaction, complex

    def text_return_protein_dna_interactions(self):
        pass

    def text_return_dna_protein_interactions(self):
        pass
