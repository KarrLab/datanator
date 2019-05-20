#!/bin/bash

MONGO_DATABASE="datanator"
APP_NAME="datanator"

MONGO_HOST="rs0/mongo:27017,mongo_secondary:27017,mongo_tertiary:27017"
TIMESTAMP=`date +%F-%H%M`
MONGORESTORE_PATH="/usr/bin/mongorestore"
BACKUPS_DIR="/root/host/karr_lab/datanator/datanator/data_source/cache/$MONGO_DATABASE"

$MONGORESTORE_PATH -d $MONGO_DATABASE BACKUPS_DIR