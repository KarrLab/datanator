import sys
sys.path.append("/Users/pochis01/Desktop/GitHub/datanator")
import datetime

# old_stdout = sys.stdout
# log_file = open("datanator/builds/logs/{}.txt".format(str(datetime.datetime.now())),"w")
# sys.stdout = log_file


from datanator.core import common_schema
cs = common_schema.CommonSchema(load_content=True, clear_content=True, verbose=True, test=True, max_entries=10)

# sys.stdout = old_stdout
# log_file.close()
