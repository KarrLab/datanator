from datanator.core import data_model, data_query, common_schema, models
from datanator.api.lib.data_manager import BaseManager
from datanator.util.constants import DATA_CACHE_DIR

class FilterManager(BaseManager):
    """ Manages filtering of information for API. Filters objects """

    def __init__(self, cache_dirname = DATA_CACHE_DIR, params=None, data=None):
        self.data_source = common_schema.CommonSchema(cache_dirname=cache_dirname)
        self.params = params
        self.data = data

        self.run()


    def run(self):
        pass
