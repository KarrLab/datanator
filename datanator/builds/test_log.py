import sys
import datetime
old_stdout = sys.stdout
log_file = open("datanator/builds/logs/{}.txt".format(str(datetime.datetime.now())),"w")
sys.stdout = log_file


sys.stdout = old_stdout
log_file.close()
