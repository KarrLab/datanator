import sys
sys.path.append("/Users/pochis01/Desktop/GitHub/datanator")
import datetime
from datanator.core import common_schema
old_stdout = sys.stdout
log_file = open("datanator/builds/logs/{}.txt".format(str(datetime.datetime.now())),"w")
sys.stdout = log_file

cs = common_schema.CommonSchema(load_content=True, verbose=True, load_entire_small_dbs=True)
cs.upload_backup()

sys.stdout = old_stdout
log_file.close()
