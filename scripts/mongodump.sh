#!/bin/bash
 
MONGO_DATABASE="datanator"
APP_NAME="datanator"

MONGO_HOST="rs0/mongo:27017,mongo_secondary:27017,mongo_tertiary:27017"
TIMESTAMP=`date +%F-%H%M`
MONGODUMP_PATH="/usr/bin/mongodump"
BACKUPS_DIR="/tmp"
BACKUP_NAME="$APP_NAME-$TIMESTAMP"
 
$MONGODUMP_PATH -d $MONGO_DATABASE --host $MONGO_HOST 
 
mkdir -p $BACKUPS_DIR
mv dump $BACKUP_NAME
tar -zcvf $BACKUPS_DIR/$BACKUP_NAME.tgz $BACKUP_NAME
rm -rf $BACKUP_NAME