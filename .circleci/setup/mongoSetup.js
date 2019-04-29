
rs.initiate();
cfg = {
  _id: 'datanator',
  members: [
    { _id: 0, host: 'mongo-datanator-1:27017' },
    { _id: 1, host: 'mongo-datanator-2:27018' },
    { _id: 2, host: 'mongo-datanator-3:27019' }
  ]
};
cfg.protocolVersion = 1;
rs.reconfig(cfg, { force: true });