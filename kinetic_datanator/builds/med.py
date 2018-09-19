import sys
sys.path.append("/Users/pochis01/Desktop/GitHub/kinetic_datanator")
import datetime

old_stdout = sys.stdout
log_file = open("kinetic_datanator/kinetic_datanator/builds/logs/{}.txt".format(str(datetime.datetime.now())),"w")
sys.stdout = log_file

from kinetic_datanator.core import common_schema
cs = common_schema.CommonSchema(load_content=True, clear_content=True, verbose=True, max_entries=20, load_entire_small_DBs=True)
# cs.dump_database()
# cs.upload_backup()

sys.stdout = old_stdout
log_file.close()
