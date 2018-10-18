"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2018-08-13
:Copyright: 2018, Karr Lab
:License: MIT
"""

from flask_restplus import Api, Resource, reqparse
import json
from flask import  Blueprint, Response, render_template, make_response
from datanator.core import common_schema, models
from datanator.api.lib.search.manager import search_manager
from datanator.api.lib.metabolite.manager import metabolite_manager
from datanator.api.lib.subunit.manager import subunit_manager
from datanator.api.lib.complex.manager import complex_manager
from datanator.api.lib.reaction.manager import reaction_manager
from datanator.api.serializer import *
from datanator.util.constants import DATA_CACHE_DIR
import json
import os



api_blueprint = Blueprint('api', __name__, url_prefix='/api')
api = Api(api_blueprint, version='0.0', title='Datanator API',
    description='Providing Data for Modelers', doc='/docs/')

parser = reqparse.RequestParser()
parser.add_argument('download', type=bool, default=False)

# @api.representation('text/html')
# def output_html(data, code, headers=None):
#     resp = make_response(render_template('api/api.html', content = json.dumps(data, sort_keys=True, indent=4)), code)
#     resp.headers.extend(headers or {})
#     return resp

@api.representation('application/json')
def output_json(data, code, headers=None):
    resp = make_response(json.dumps(data, sort_keys=True, indent=4), code)
    resp.headers.extend(headers or {})
    return resp


class Search(Resource):

    @api.doc(params={'value': 'Value to search for over the database',
                    'download': 'Boolean option to download content'})
    def get(self, value):
        search_dict = search_manager.search(value)

        serialized_metabolites = MetaboliteSerializer().dump(search_dict['Metabolite'], many=True)
        serialized_complexes = ProteinComplexSerializer().dump(search_dict['ProteinComplex'], many=True)
        serialized_subunits = ProteinSubunitSerializer().dump(search_dict['ProteinSubunit'], many=True)
        serialized_reactions = ReactionSerializer().dump(search_dict['Reaction'], many=True)

        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}.json".format(value)

        return {'metabolites': serialized_metabolites.data,
                'complexes': serialized_complexes.data,
                'subunits': serialized_subunits.data,
                'reactions': serialized_reactions.data}, 200, headers

class MetaboliteSearch(Resource):
    @api.doc(params={'value': 'Value to search over in the metabolite space',
                    'download': 'Boolean option to download content'})
    def get(self,value):
        metabolite_search = metabolite_manager._search_complex(value)
        serialized_metabolites = MetaboliteSerializer().dump(metabolite_search, many=True)

        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}.json".format(value)

        return {'metabolites': serialized_metabolites.data}, 200, headers

class ProteinSubunitSearch(Resource):
    @api.doc(params={'value': 'Value to search over in the protein subunit space',
                    'download': 'Boolean option to download content'})
    def get(self,value):
        subunit_search = subunit_manager._search_complex(value)
        serialized_subunits = ProteinSubunitSerializer().dump(subunit_search, many=True)

        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}.json".format(value)

        return {'subunits': serialized_subunits.data}, 200, headers

class ProteinComplexSearch(Resource):
    @api.doc(params={'value': 'Value to search over in the protein complex space',
                    'download': 'Boolean option to download content'})
    def get(self,value):
        complex_search  = complex_manager._search(value)
        serialized_complexes = ProteinComplexSerializer().dump(complex_search, many=True)
        return {'complexes': serialized_complexes.data}, 200, headers

class Metabolite(Resource):

    @api.doc(params={'id': 'Metabolite ID to find information for',
                    'download': 'Boolean option to download content'})
    def get(self, id):
        metabolite = metabolite_manager.get_metabolite_by_id(id)
        observed_concentrations = metabolite_manager.get_observed_concentrations(metabolite)
        reactions = reaction_manager.get_reaction_by_metabolite(metabolite)

        serialized_metabolite = MetaboliteSerializer().dump(metabolite)
        serialized_concentrations = ObservedValueSerializer().dump(observed_concentrations, many=True)
        serialized_reactions =  ReactionSerializer().dump(reactions, many=True)

        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}.json".format(metabolite.metabolite_name)

        return {'object':serialized_metabolite.data,
                'concentrations': serialized_concentrations.data,
                'reactions' : serialized_reactions.data}, 200, headers

class ProteinSubunit(Resource):

    @api.doc(params={'id': 'Protein Subunit ID to find information for',
                    'download': 'Boolean option to download content'})
    def get(self, id):
        subunit = subunit_manager.get_subunit_by_id(id)
        observed_abundances = subunit_manager.get_observed_abundances(subunit)
        observed_interactions = subunit_manager.get_observable_interactions(subunit)
        observed_complexes = subunit_manager.get_observable_complex(subunit)

        serialized_subunit = ProteinSubunitSerializer().dump(subunit)
        serialized_abundances = ObservedValueSerializer().dump(observed_abundances, many=True)
        serialized_interactions = ObservedInteractionSerializer().dump(observed_interactions, many=True)
        serialized_complexes = ObservedComplexSpecieSerializer().dump(observed_complexes, many=True)

        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}.json".format(subunit.uniprot_id)


        return {'object': serialized_subunit.data,
                'abundances': serialized_abundances.data,
                'interactions': serialized_interactions.data,
                'complexes': serialized_complexes.data}, 200, headers

class ProteinComplex(Resource):

    @api.doc(params={'id': 'Protein Complex ID to find information for',
                    'download': 'Boolean option to download content'})
    def get(self, id):
        complex = complex_manager.get_complex_by_id(id)
        observed_subunits = complex_manager.get_observable_subunits(complex)

        serialized_complex = ProteinComplexSerializer().dump(complex)
        serialized_subunits = ObservedProteinSpecieSerializer().dump(observed_subunits, many=True)


        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}.json".format(complex.complex_name)


        return {'object':serialized_complex.data,
                'subunits': serialized_subunits.data}, 200, headers

class Reaction(Resource):

    @api.doc(params={'id': 'Reaction ID to find information for',
                    'download': 'Boolean option to download content'})
    def get(self, id):
        reaction = reaction_manager.get_reaction_by_kinetic_law_id(id)
        observed_parameters = reaction_manager.get_observed_parameter_value(reaction)

        serialized_reaction = ReactionSerializer().dump(reaction)
        serialized_parameters = ObservedValueSerializer().dump(observed_parameters, many=True)

        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}_reaction.json".format(reaction.id)


        return {'object': serialized_reaction.data,
                'parameters': serialized_parameters.data}, 200, headers

class MetaboliteConcentration(Resource):

    @api.doc(params={'id': 'Metabolite ID to find concentration information for',
                    'download': 'Boolean option to download content'})
    def get(self, id):
        metabolite = metabolite_manager.get_metabolite_by_id(id)
        observed_concentrations = metabolite_manager.get_observed_concentrations(metabolite)
        serialized_concentrations = ObservedValueSerializer().dump(observed_concentrations, many=True)

        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}_concentrations.json".format(metabolite.metabolite_name)


        return {'concentrations': serialized_concentrations}, 200, headers

class ProteinAbundance(Resource):

    @api.doc(params={'id': 'Protein Subunit ID to find abundance information for',
                    'download': 'Boolean option to download content'})
    def get(self, id):
        subunit = subunit_manager.get_subunit_by_id(id)
        observed_abundances = subunit_manager.get_observed_abundances(subunit)
        serialized_abundances =  ObservedValueSerializer().dump(observed_abundances, many=True).data

        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}_abundances.json".format(subunit.uniprot_id)


        return {'abundances':serialized_abundances}, 200, headers

class ProteinInteraction(Resource):

    def get(self,id):
        pass

class ReactionParameter(Resource):

    @api.doc(params={'id': 'Reaction ID to find reaction paramater information for',
                    'download': 'Boolean option to download content'})
    def get(self,id):
        reaction = reaction_manager.get_reaction_by_kinetic_law_id(id)
        observed_parameters = reaction_manager.get_observed_parameter_value(reaction)
        serialized_parameters = ObservedValueSerializer().dump(observed_parameters, many=True).data

        headers = {}
        if parser.parse_args()['download'] == True:
            headers['Content-Disposition'] = "attachment; filename={0}_rxn_parameters.json".format(reaction.id)

        return {'parameters': serialized_parameters}, 200, headers
