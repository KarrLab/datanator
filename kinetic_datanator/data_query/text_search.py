#-*- coding: utf-8 -*-

""" Test of common schema

:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-12-18
:Copyright: 2017, Karr Lab
:License: MIT
"""


from kinetic_datanator.data_query import dna_protein_interactions, metabolite_concentrations, protein_protein_interactions, reaction_kinetics
from kinetic_datanator.flask_datanator import flask_common_schema, models
import flask_whooshalchemy

class TextSearchSession(object):
    """
    Represents a text search through the database grabbing relevant information

    """

    def __init__(self):
        self.db = flask_common_schema.FlaskCommonSchema()

        flask_whooshalchemy.whoosh_index(models.app, models.Compound)
        flask_whooshalchemy.whoosh_index(models.app, models.ProteinComplex)
        flask_whooshalchemy.whoosh_index(models.app, models.ProteinInteractions)
        flask_whooshalchemy.whoosh_index(models.app, models.Taxon)
        flask_whooshalchemy.whoosh_index(models.app, models.Synonym)
        flask_whooshalchemy.whoosh_index(models.app, models.CellLine)
        flask_whooshalchemy.whoosh_index(models.app, models.CellCompartment)
        flask_whooshalchemy.whoosh_index(models.app, models.ProteinSubunit)



    def text_search(self, string):
        ans = {}

        ans['Compound'] = models.Compound.query.whoosh_search(string).all()
        ans['ProteinComplex'] = models.ProteinComplex.query.whoosh_search(string).all()
        ans['ProteinInteractions'] = models.ProteinInteractions.query.whoosh_search(string).all()
        ans['Taxon'] = models.Taxon.query.whoosh_search(string).all()
        ans['Synonym'] = models.Synonym.query.whoosh_search(string).all()
        ans['CellLine'] = models.CellLine.query.whoosh_search(string).all()
        ans['CellCompartment'] = models.CellCompartment.query.whoosh_search(string).all()
        ans['ProteinSubunit'] = models.ProteinSubunit.query.whoosh_search(string).all()


        return ans
