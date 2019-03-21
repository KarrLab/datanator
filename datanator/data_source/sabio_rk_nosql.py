'''Parse SabioRk json files into MongoDB documents
'''

import os
import json
from pymongo import MongoClient

class SabioRkNoSQL():
	def __init__(self, directory):
		self.directory = directory