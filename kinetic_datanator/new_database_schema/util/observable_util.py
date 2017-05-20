"""
:Author: Yosef Roth <yosefdroth@gmail.com>
:Author: Jonathan Karr <jonrkarr@gmail.com>
:Date: 2017-05-17
:Copyright: 2017, Karr Lab
:License: MIT
"""


#What I should do is combine molecule and reaction and factor out whatever I can into observable.m



class Observable:
	""" Represents all cell components. This includes molecules(e.g. metabolites, DNA, RNA, protein) and events (e.g. reaction). 
			Observable is the superclass from which the more specific subclasses inherit from (e.g. metabolite, reaction)

	Attributes:
		name (:obj:`string`): a canonical name for the observable (e.g. ATP)
		id (:obj:'dict' of :obj:`str`: :obj:`str`): dictionary of ids (keys = name of the database (e.g. KEGG), values = identifier)
		type (:obj: 'string'): the observables type (e.g. protein, reaction, etc)
	"""

	RECOGNIZED_TYPES = ['metabolite', 'DNA', 'RNA', 'protein', 'reaction']

	def __init__(self, name='', id={}, obs_type=''):
		"""
		Args:
			name (:obj:`string`): a canonical name for the observable (e.g. ATP)
			id (:obj:'dict' of :obj:`str`: :obj:`str`): dictionary of ids (keys = name of the database (e.g. KEGG), values = identifier)
			type (:obj: 'string'): the observables type (e.g. protein, reaction, etc)servations from mutants
		
		Raises:
			:obj:`ValueError`: if a type is entered that is not also in the RECOGNIZED_TYPES list

		""" 
		self.name = name 
		self.id = id
		if obs_type in self.RECOGNIZED_TYPES or obs_type=='':
			self.obs_type = obs_type
		else:
			raise ValueError

	def get_id(self, id_name):
		"""Return an identifier

		Args:
			id_name (:obj:'str'): identification name. Usually associated with a database (e.g. KEGG, PUBCHEM, CHEBI)

		Returns:
			:obj:'str':the idenfifier (usually a number, e.g. 30128)

		"""
		try:
			return self.id[id_name]
		except:
			return ""
