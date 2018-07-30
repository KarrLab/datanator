import os


## Paths
BASEDIR = os.path.abspath(os.path.dirname(__file__))
DATA_CACHE_DIR = os.path.join(BASEDIR, '..', 'data', 'cache')



## Batching Test Constants
PAX_TEST_BATCH = 2
INTACT_INTERACTION_TEST_BATCH = 10
ARRAY_EXPRESS_TEST_BATCH = 10
SABIO_TEST_BATCH = 100




## Batching Build Constants
PAX_BUILD_BATCH = 5
INTACT_INTERACTION_BUILD_BATCH = 10000
ARRAY_EXPRESS_BUILD_BATCH = 20
SABIO_BUILD_BATCH = 10000
