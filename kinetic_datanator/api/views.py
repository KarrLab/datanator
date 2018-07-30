"""
:Author: Saahith Pochiraju <saahith116@gmail.com>
:Date: 2017-02-28
:Copyright: 2018, Karr Lab
:License: MIT
"""

from flask_restful import Api, Resource, reqparse
from flask import  Blueprint, jsonify, Response
from kinetic_datanator.core import common_schema, models
import json
import os

cachedir = os.path.join(os.path.abspath(os.path.dirname(__file__)),'..', 'cache')
# flk = common_schema.CommonSchema(cache_dirname = cachedir)

class Documentation(Resource):
    """
    Represents the documentation page of the API

    """
    def get(self):
        return {
            'observation_info':  '/api/v1.0/Observation',
            'physical_entity_info': '/api/v1.0/PhysicalEntity',
            'physical_property_info': '/api/v1.0/PhysicalProperty',
            'compound_table_info': '/api/v1.0/Compound',
            'compound_structure_info': '/api/v1.0/Structure',
            'protein_subunit_info': '/api/v1.0/ProteinSubunit',
            'protein_complex_info': '/api/v1.0/ProteinComplex',
            'concentration_info':  '/api/v1.0/Concentration',
            'protein_interaction_info':  '/api/v1.0/ProteinInteraction'
        }

# class DataDump(Resource):
#     """
#     Represents the API data dump from the common schema database
#
#     """
#     #TODO: Account for other filetypes/Cache data tables
#
#     def get(self, table):
#         type_mappings = {'INTEGER': int, 'VARCHAR(255)': str, 'BOOLEAN': bool, 'FLOAT': float, 'VARCHAR': str}
#
#         table_obj = getattr(model, table)
#
#         parser = reqparse.RequestParser()
#         parser.add_argument('download', type=bool)
#         parser.add_argument('format', type=str)
#         for c in table_obj.__table__.columns:
#             parser.add_argument(str(c.name), type=type_mappings[str(c.type)])
#
#         ## Full Table Return
#         if all(value == None for value in parser.parse_args().values()):
#             with open(cachedir+'/'+str(table)+'.json') as json_data:
#                 return json.load(json_data)
#
#         ##Filtered Return and Download
#         results = flk.session.query(table_obj)
#         for key,value in parser.parse_args().items():
#             if value == None or key=='download' or key=='format':
#                 continue
#             condition = getattr(table_obj, key) == value
#             results = results.filter(condition)
#
#         ans = {}
#         for item in results.all():
#             ans[item.id] = item.serialize()
#
#         json_return = jsonify({table: ans})
#
#         if parser.parse_args()['download'] == True and parser.parse_args()['format'] == 'json' :
#             json_return.headers['Content-Disposition'] = "attachment; filename="+str(table)+'.json'
#
#         return json_return
