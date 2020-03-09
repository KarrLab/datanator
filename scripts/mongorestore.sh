#!/bin/bash
FILE=/root/karr_lab/datanator.archive
if [ ! -f "$FILE" ]; then
	curl -o $FILE https://mongo-dbdump.s3.amazonaws.com/datanator.20190701.archive
fi
mongorestore -d datanator --host mongo:27017 --archive=$FILE