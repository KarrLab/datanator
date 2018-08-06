from kinetic_datanator.api.views import *

V0_ENDPOINT = '/v0'


# Queries
api.add_resource(TextSearch, V0_ENDPOINT+'/search/<value>')
