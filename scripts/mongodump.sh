#!/bin/bash
 
MONGO_DATABASE="datanator"
APP_NAME="datanator"

MONGO_HOST="rs0/mongo:27017,mongo_secondary:27017,mongo_tertiary:27017"
TIMESTAMP=`date +%F-%H%M`
MONGODUMP_PATH="/usr/bin/mongodump"
BACKUPS_DIR="/root/host/karr_lab/datanator/datanator/data_source/cache"
mkdir -p $BACKUPS_DIR
BACKUP_NAME="/root/host/karr_lab/datanator/datanator/data_source/cache/$APP_NAME-$TIMESTAMP"
 
$MONGODUMP_PATH -d $MONGO_DATABASE --host $MONGO_HOST -o $BACKUP_NAME

