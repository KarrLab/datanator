from kinetic_datanator.api.views import *

V0_ENDPOINT = '/v0'
# Text Search
api.add_resource(Search, V0_ENDPOINT+'/search/<value>')
api.add_resource(Metabolite, V0_ENDPOINT+'/search/metabolite/<value>')
api.add_resource(ProteinSubunit, V0_ENDPOINT+'/search/subunit/<value>')
api.add_resource(Concentration, V0_ENDPOINT+'/concentrations/<id>')

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
