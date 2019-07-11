#!/bin/bash

MONGO_DATABASE="datanator"
APP_NAME="datanator"

MONGO_HOST="mongo"
TIMESTAMP=`date +%F-%H%M`
MONGORESTORE_PATH="/usr/bin/mongorestore"
BACKUPS_DIR="/root/host/karr_lab/datanator.20190701.archive"

$MONGORESTORE_PATH -d $MONGO_DATABASE --host $MONGO_HOST --archive=BACKUPS_DIR