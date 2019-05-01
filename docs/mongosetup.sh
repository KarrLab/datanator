#!/bin/bash

mongo --host mongo:27017 <<EOF
   var cfg = {
        "_id": "rs0",
        "version": 1,
        protocolVersion: 1,
        "members": [
            {
                "_id": 0,
                "host": "mongo:27017",
                "priority": 2
            },
            {
                "_id": 1,
                "host": "mongo_secondary:27017",
                "priority": 0.5
            },
            {
                "_id": 2,
                "host": "mongo_tertiary:27017",
                "priority": 0.5
            }
        ]
    };
    rs.initiate(cfg, { force: true });
    rs.reconfig(cfg, { force: true });
    rs.slaveOk();
    db.getMongo().setReadPref('nearest');
    db.getMongo().setSlaveOk(); 
EOF
