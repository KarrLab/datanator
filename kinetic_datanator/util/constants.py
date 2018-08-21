import os

## Paths
BASEDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_CACHE_DIR = os.path.join(BASEDIR, 'kinetic_datanator', 'data', 'cache')

# Common Schema Constants
DATA_DUMP_PATH = os.path.join(DATA_CACHE_DIR , 'CommonSchema.dump')
PAX_NAME = 'Pax'
PAX_INITIAL_AMOUNT = 1
SABIO_NAME = 'Sabio'
SABIO_INITIAL_AMOUNT = 1
ARRAY_EXPRESS_NAME = 'Array Express'
ARRAY_EXPRESS_INITIAL_AMOUNT = 1
INTACT_NAME = 'IntAct'
INTACT_INITIAL_AMOUNT = 0

## Batching Test Constants
PAX_TEST_BATCH = 2
INTACT_INTERACTION_TEST_BATCH = 10
ARRAY_EXPRESS_TEST_BATCH = 5
SABIO_TEST_BATCH = 100

## Batching Build Constants
PAX_BUILD_BATCH = 20
INTACT_INTERACTION_BUILD_BATCH = 200000
ARRAY_EXPRESS_BUILD_BATCH = 40
SABIO_BUILD_BATCH = 40000
