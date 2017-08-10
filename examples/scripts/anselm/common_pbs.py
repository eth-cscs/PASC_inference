#
# 
#

# include some funny functions
import os, shutil
from subprocess import call

# define function for writing a fun into file
def write_pbs(problem_name, host_string, batch_filename, exec_path, exec_name, modules_path):
    "this function prints a fun into batch script file, the fun is based on parameters"
    myfile = open(batch_filename, 'w+');
    # write some funny stuff into file
    myfile.write("#!/bin/bash --login\n")
    myfile.write("#PBS -A IT4I-7-5\n")
    myfile.write("#PBS -q qexp\n")
    myfile.write("#PBS -N %s\n" %(problem_name) ) 
    myfile.write("#PBS -l %s\n" %(host_string) )
    myfile.write("#PBS -j oe\n\n")
    myfile.write("cd %s\n" %(exec_path) )
    myfile.write("source %s\n\n" %(modules_path) )
    myfile.write("%s\n" %(exec_name))
    return

def commit_pbs(batchfile_list):
    "this function commits PBS batch files"
    # say what we are doing now:
    print "Commiting PBS scripts: "
	# send every batch file from the list
    for batchfile_name in batchfile_list:
        print  " - %s" % (batchfile_name);
        call(["qsub", batchfile_name])
    return

def show_jobs(account):
    "this function shows PBS queue"
    print "--------------------------------------- MY JOBS: -----------------------------------"
    call(["qstat", "-u %s" % (account)])
    return






