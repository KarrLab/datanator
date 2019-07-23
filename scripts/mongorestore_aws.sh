#!/bin/bash
FILE=~/karr_lab/datanator-dump
source <(grep = ~/.wc/datanator.cfg | tr -d ' ')
if [ ! -d "$FILE" ]; then
	mkdir $FILE
	aws s3 cp https://mongo-dbdump.s3.amazonaws.com/datanator $FILE --recursive
fi
mongorestore -d datanator -u $user -p $password --authenticationDatabase admin "$FILE/datanator"