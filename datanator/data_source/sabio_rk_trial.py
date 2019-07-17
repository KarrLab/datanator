import datanator.config.core
from datanator.util import mongo_util
import six

class SabioRk:

    def __init__(self, cache_dirname=None, MongoDB=None, replicaSet=None, db=None,
                 verbose=False, max_entries=float('inf'), username=None,
                 password=None, authSource='admin'):
        self.cache_dirname = cache_dirname
        self.MongoDB = MongoDB
        self.replicaSet = replicaSet
        self.db = db
        self.verbose = verbose
        self.max_entries = max_entries
        self.client, self.db_obj, self.collection = mongo_util.MongoUtil(
            MongoDB=MongoDB, db=db, username=username, password=password, authSource=authSource)
        ENDPOINT_DOMAINS = {
            'sabio_rk': 'http://sabiork.h-its.org',
            'uniprot': 'http://www.uniprot.org',
        }
        ENDPOINT_KINETIC_LAWS_SEARCH = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/sabioRestWebServices/searchKineticLaws/entryIDs'
        ENDPOINT_WEBSERVICE = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/sabioRestWebServices/kineticLaws'
        ENDPOINT_EXCEL_EXPORT = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/entry/exportToExcelCustomizable'
        ENDPOINT_COMPOUNDS_PAGE = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/compdetails.jsp'
        ENDPOINT_KINETIC_LAWS_PAGE = ENDPOINT_DOMAINS['sabio_rk'] + \
            '/kindatadirectiframe.jsp'
        PUBCHEM_MAX_TRIES = 10
        PUBCHEM_TRY_DELAY = 0.25


    def load_kinetic_law_ids(self):
        """ Download the IDs of all of the kinetic laws stored in SABIO-RK

        Returns:
            :obj:`list` of :obj:`int`: list of kinetic law IDs

        Raises:
            :obj:`Error`: if an HTTP request fails or the expected number of kinetic laws is not returned
        """
        # create session
        session = self.requests_session
        response = session.get(self.ENDPOINT_KINETIC_LAWS_SEARCH, params={
            'q': 'DateSubmitted:01/01/2000',
        })
        response.raise_for_status()

        # get IDs of kinetic laws
        root = etree.ElementTree.fromstring(response.text)
        ids = [int(float(node.text)) for node in root.findall('SabioEntryID')]

        # sort ids
        ids.sort()

        # return IDs
        return ids