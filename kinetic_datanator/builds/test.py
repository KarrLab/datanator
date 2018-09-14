import sys
sys.path.append("/Users/pochis01/Desktop/GitHub/kinetic_datanator")
import datetime

# old_stdout = sys.stdout
# log_file = open("kinetic_datanator/builds/logs/{}.txt".format(str(datetime.datetime.now())),"w")
# sys.stdout = log_file


from kinetic_datanator.core import common_schema
cs = common_schema.CommonSchema(load_content=True, clear_content=True, verbose=True, test=True, max_entries=10)

# sys.stdout = old_stdout
# log_file.close()
