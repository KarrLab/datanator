from datanator.core import common_schema
import datetime
import pkg_resources
import sys


def build():
    old_stdout = sys.stdout
    log_filename = pkg_resources.resource_filename(
        'datanator', "builds/logs/{}.txt".format(str(datetime.datetime.now())))
    with open(log_filename, "w") as log_file:
        sys.stdout = log_file
        cs = common_schema.CommonSchema(load_content=True, verbose=True,
                                        load_entire_small_dbs=True)
        cs.upload_backup()
    sys.stdout = old_stdout
