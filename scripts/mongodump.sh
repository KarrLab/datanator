#!/bin/bash
 
MONGO_DATABASE="datanator"
APP_NAME="datanator"
source <(grep = ~/.wc/datanator.cfg | tr -d ' ')

MONGO_HOST="localhost"
TIMESTAMP=`date +%F-%H%M`
MONGODUMP_PATH="/usr/bin/mongodump"
BACKUPS_DIR="/data/mongodump_$TIMESTAMP"
mkdir -p $BACKUPS_DIR
BACKUP_NAME="$BACKUPS_DIR/$APP_NAME"
 
$MONGODUMP_PATH -d $MONGO_DATABASE -u $user -p $password --authenticationDatabase admin -o $BACKUP_NAME
aws s3 cp $BACKUP_NAME s3://mongo-dbdump/ --recursive
rm -rf $BACKUPS_DIR