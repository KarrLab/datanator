
var http = require('http');
var mongoose = require('mongoose');

mongoose.connect(
  'mongodb://mongo-datanator-1:27017,mongo-datanator-2:27017,mongo-datanator-3:27017/test?replicaSet=datanator'
);

//create a server object:
http
  .createServer(function(req, res) {
    res.write('echo');
    res.end();
  })
.listen(3001); //the server object listens on port 3001