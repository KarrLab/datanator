class Super:


	def __init__(self, b="", g=""):
		self.blue = b
		self.green = g




class Sub(Super):
	def __init__(self):
		Super.__init__(self, b="blue", g="green")

a = Sub()
print(a.blue)
print(a.green)




import openbabel
import pybel
import re
from kinetic_datanator.new_database_schema.util import observable_util

class Molecule(observable_util.Observable):
	""" Represents a molecule

	Attributes:
		id (:obj:`str`): identifier
		name (:obj:`str`): name
		structure (:obj:`str`): structure in InChI, MOL, or canonical SMILES format
	"""

	# todo: integrate with :obj:`InchiMolecule`

	def __init__(self, name='', id={}):#, structure=''):
		"""
		Args:
			id (:obj:`str`, optional): identifier
			name (:obj:`str`, optional): name
			structure (:obj:`str`, optional): structure in InChI, MOL, or canonical SMILES format
			cross_references (:obj:`list` of :obj:`CrossReference`, optional): list of cross references

		Raises:
			:obj:`ValueError`: if the structure is not valid
		"""
		print(name)
		#print(obs_type)
		observable_util.Observable.__init__(self, name, id, "")#self, name, id)#, "")
		#self.id = id
		#self.name = name
		#self.structure = structure
a = observable_util.Observable()
print(a.obs_type)
a = Molecule()
#a = Molecule('ATP', {"KEGG": "543543"})
