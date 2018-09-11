from kinetic_datanator.core import common_schema
cs = common_schema.CommonSchema(load_content=True, clear_content=True, verbose=True, max_entries=20, load_entire_small_DBs=True)
cs.dump_database()
cs.upload_backup()
