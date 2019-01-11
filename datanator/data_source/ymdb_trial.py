
import datetime
import io
import json

class Ymdb(object):
	"""load YMDB JSON"""
	def __init__(self, arg):
		super(Ymdb, self).__init__()
		self.arg = arg
		
	def load_content(self):
	    """ Download the content of YMDB and store it to a local sqlite database. """
	    db_session = self.session
	    req_session = self.requests_session

	    # download content from server
	    if self.verbose:
	        print('Downloading compound IDs ...')

	    response = req_session.get(self.DOWNLOAD_FULL_DB_URL)
	    response.raise_for_status()

	    if self.verbose:
	        print('  done')

	    # unzip and parse content
	    if self.verbose:
	        print('Parsing compound IDs ...')

	    with zipfile.ZipFile(io.BytesIO(response.content), 'r') as zip_file:
	        with zip_file.open('ymdb.json', 'r') as json_file:
	            entries = json.load(json_file)

	    if self.verbose:
	        print('  found {} compounds'.format(len(entries)))

