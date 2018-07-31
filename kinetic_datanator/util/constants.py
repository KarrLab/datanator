import os

## Paths
BASEDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_CACHE_DIR = os.path.join(BASEDIR, 'kinetic_datanator', 'data', 'cache')

# Common Schema Constants
DATA_DUMP_PATH = os.path.join(DATA_CACHE_DIR , 'CommonSchema.dump')


## Batching Test Constants
PAX_TEST_BATCH = 2
INTACT_INTERACTION_TEST_BATCH = 10
ARRAY_EXPRESS_TEST_BATCH = 5
SABIO_TEST_BATCH = 100

## Batching Build Constants
PAX_BUILD_BATCH = 5
INTACT_INTERACTION_BUILD_BATCH = 10000
ARRAY_EXPRESS_BUILD_BATCH = 20
SABIO_BUILD_BATCH = 10000
