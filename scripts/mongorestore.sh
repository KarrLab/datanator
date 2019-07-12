#!/bin/bash
curl -o /root/host/karr_lab/datanator.archive https://mongo-dbdump.s3.amazonaws.com/datanator.20190701.archive
mongorestore -d datanator --host mongo:27017 --archive=/root/host/karr_lab/datanator.archive