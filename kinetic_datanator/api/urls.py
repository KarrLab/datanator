from kinetic_datanator.api.views import *
from kinetic_datanator.util.constants import CURRENT_VERSION_ENDPOINT

def get_version_endpoint(endpoint):
    return CURRENT_VERSION_ENDPOINT+endpoint

# Text Search
api.add_resource(Search, get_version_endpoint('/search/<value>'))
api.add_resource(Metabolite, get_version_endpoint('/search/metabolite/<value>'))
api.add_resource(ProteinSubunit, get_version_endpoint('/search/subunit/<value>'))
api.add_resource(ProteinComplex, get_version_endpoint('/search/complex/<value>'))

api.add_resource(Concentration, get_version_endpoint('/concentrations/<id>'))

#
# """
# Compound, /compound
# Protein Subunit, /subunit
# Protein Complex, /complex
# Reaction, /reaction
# ProteinInteractions, /interactions
# Abundance, /abundance
#
# """
