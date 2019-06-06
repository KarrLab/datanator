# """
# :Author: Saahith Pochiraju <saahith116@gmail.com>
# :Date: 2018-08-13
# :Copyright: 2018, Karr Lab
# :License: MIT
# """

# from flask_restplus import Api, Resource, reqparse
# import json
# from flask import  Blueprint, Response, render_template, make_response
# import json
# import os



# api_blueprint = Blueprint('api', __name__, url_prefix='/api')
# api = Api(api_blueprint, version='0.0', title='Datanator API',
#     description='Providing Data for Modelers', doc='/docs/')

# parser = reqparse.RequestParser()
# parser.add_argument('download', type=bool, default=False)

# # @api.representation('text/html')
# # def output_html(data, code, headers=None):
# #     resp = make_response(render_template('api/api.html', content = json.dumps(data, sort_keys=True, indent=4)), code)
# #     resp.headers.extend(headers or {})
# #     return resp

# @api.representation('application/json')
# def output_json(data, code, headers=None):
#     pass

# class Search(Resource):

#     pass

# class MetaboliteSearch(Resource):
#     pass

#     return {'metabolites': serialized_metabolites.data}, 200, headers

# class ProteinSubunitSearch(Resource):
#     pass

#     return {'subunits': serialized_subunits.data}, 200, headers

# class ProteinComplexSearch(Resource):
#     pass
#     return {'complexes': serialized_complexes.data}, 200, headers

# class Metabolite(Resource):

#     pass
#     return {'object':serialized_metabolite.data,
#                 'concentrations': serialized_concentrations.data,
#                 'reactions' : serialized_reactions.data}, 200, headers

# class ProteinSubunit(Resource):

#     pass


#     return {'object': serialized_subunit.data,
#                 'abundances': serialized_abundances.data,
#                 'interactions': serialized_interactions.data,
#                 'complexes': serialized_complexes.data}, 200, headers

# class ProteinComplex(Resource):

#     pass
#     return {'object':serialized_complex.data,
#                 'subunits': serialized_subunits.data}, 200, headers

# class Reaction(Resource):

#     pass
#     return {'object': serialized_reaction.data,
#                 'parameters': serialized_parameters.data}, 200, headers

# class MetaboliteConcentration(Resource):

#     pass

#     return {'concentrations': serialized_concentrations}, 200, headers

# class ProteinAbundance(Resource):

#     pass
#     return {'abundances':serialized_abundances}, 200, headers

# class ProteinInteraction(Resource):

#     def get(self,id):
#         pass

# class ReactionParameter(Resource):

#     pass
#     return {'parameters': serialized_parameters}, 200, headers
