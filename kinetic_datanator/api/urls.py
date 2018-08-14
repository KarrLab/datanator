from kinetic_datanator.api.views import *

V0_ENDPOINT = '/v0'


# Text Search
api.add_resource(Search, V0_ENDPOINT+'/search/<value>')


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
