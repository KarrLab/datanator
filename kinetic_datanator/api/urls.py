from kinetic_datanator.api.views import *
from flask_restful import Api
from flask import  Blueprint, jsonify, Response

api_blueprint = Blueprint('api', __name__,)
api = Api(api_blueprint)

api.add_resource(Documentation, '/v0/docs')
# api.add_resource(DataDump, '/v0/data/<table>')
