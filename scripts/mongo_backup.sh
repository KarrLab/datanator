#/bin/bash
cd /opt/mongodb/mongodb/bin/
echo `date` >>/mongo_data/backup/backup.log
APP_NAME="datanator"
MONGO_HOST="mongo"
MONGO_PORT="27017"
TIMESTAMP=`date +%F-%H%M`
MONGODUMP_PATH="/usr/bin/mongodump"
BACKUPS_DIR="/mongo_data/backup/$APP_NAME-$TIMESTAMP"
BACKUP_NAME="/mongo_data/backup/$APP_NAME-$TIMESTAMP"
mkdir -p $BACKUPS_DIR
cd /opt/mongodb/mongodb/bin/
#Delete all backups older than 30 days from /mongo_data/backup
echo "Deleting following backup files older than 30 days:" >>
    /mongo_data/backup/backup.log
find /mongo_data/backup/ -type d -name 'app1-*' -mtime +30 >>
    /mongo_data/backup/backup.log
find /mongo_data/backup/ -type d -name 'app1-*' -mtime +30 -exec rm -rf {}
    +
#Run the daily backup 'local' database only.
for databaseName in local
do
echo "Starting daily backup of $databaseName ...." >>
    /mongo_data/backup/backup.log
./mongodump --ssl --sslCAFile ../cert/mongo.server.trust-certs.pem
    --sslPEMKeyPassword password123 --host pre-mongo01.ibmcloud.com:27017 --db
    $databaseName >>/mongo_data/backup/backup.log
#Run the daily backup of remaining databases.
echo "Starting daily backup of all databases...." >>
    /mongo_data/backup/backup.log
./mongodump --ssl --sslCAFile ../cert/mongo.server.trust-certs.pem
    --sslPEMKeyPassword password123 --host pre-mongo01.ibmcloud.com:27017
    >>/mongo_data/backup/backup.log
if [ $? != 0 ]; then
echo "Failed to make backup of $databaseName on `date +%F_%T`"|mailx -s
    "MongoDB backup failed" amolbarsagade@in.ibm.com
fi
done
mv /opt/mongodb/mongodb/bin/dump $BACKUP_NAME
echo `date` >> /mongo_data/backup/backup.log
echo "End of backup run" >> /mongo_data/backup/backup.log
echo "----------------------------------" >>
    /mongo_data/backup/backup.log