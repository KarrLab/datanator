import sys
sys.path.append("/Users/pochis01/Desktop/GitHub/datanator")
import datetime

old_stdout = sys.stdout
<<<<<<< HEAD:kinetic_datanator/builds/med.py
log_file = open("kinetic_datanator/builds/logs/{}.txt".format(str(datetime.datetime.now())),"w")
sys.stdout = log_file

from kinetic_datanator.core import common_schema
cs = common_schema.CommonSchema(load_content=True, max_entries=20, load_entire_small_dbs=False, verbose=True)
=======
log_file = open("datanator/datanator/builds/logs/{}.txt".format(str(datetime.datetime.now())),"w")
sys.stdout = log_file

from datanator.core import common_schema
cs = common_schema.CommonSchema(load_content=True, clear_content=True,
                                max_entries=20, load_entire_small_dbs=True, verbose=True)
>>>>>>> 650d74fcc6d11695e9a12ba7db13ea71ad681d51:datanator/builds/med.py
# cs.upload_backup()

sys.stdout = old_stdout
log_file.close()
