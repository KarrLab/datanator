from datanator.api.views import *
from datanator.util.constants import CURRENT_VERSION_ENDPOINT

def get_version_endpoint(endpoint):
    return CURRENT_VERSION_ENDPOINT+endpoint

# Text Search
api.add_resource(Search, get_version_endpoint('/search/<value>'))
api.add_resource(MetaboliteSearch, get_version_endpoint('/search/metabolite/<value>'))
api.add_resource(ProteinSubunitSearch, get_version_endpoint('/search/subunit/<value>'))
api.add_resource(ProteinComplexSearch, get_version_endpoint('/search/complex/<value>'))

# Object Specific Queries
api.add_resource(Metabolite, get_version_endpoint('/metabolite/<id>'))
api.add_resource(ProteinSubunit, get_version_endpoint('/subunit/<id>'))
api.add_resource(ProteinComplex, get_version_endpoint('/complex/<id>'))
api.add_resource(Reaction, get_version_endpoint('/reaction/<id>'))

# Data Specific Queries
api.add_resource(MetaboliteConcentration, get_version_endpoint('/concentrations/<id>'))
api.add_resource(ProteinAbundance, get_version_endpoint('/abundances/<id>'))
api.add_resource(ProteinInteraction, get_version_endpoint('/interactions/<id>'))
api.add_resource(ReactionParameter, get_version_endpoint('/parameters/<id>'))
