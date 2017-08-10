#
# 
#

# include some funny functions
import os, shutil
from subprocess import call

# define function for writing a fun into file
def write_batch(problem_name, nnodes, ntaskspernode, nthreads, time, library_path, build_path, exec_list, module_name):
    "this function prints a fun into batch script file, the fun is based on parameters"
    problem_name_full = "%s" % (problem_name)
    batchfile_name = "%s/%s.batch" % (build_path,problem_name_full);
    myfile = open(batchfile_name, 'w+');
    # write some funny stuff into file
    myfile.write("#!/bin/bash -l\n")
    myfile.write("\n## sbatch settings\n")
    myfile.write("#SBATCH --nodes=%d\n" % (nnodes))
    myfile.write("#SBATCH --ntasks-per-node=%d\n" % (ntaskspernode))
    myfile.write("#SBATCH --ntasks-per-core=1\n")
    myfile.write("#SBATCH --threads-per-core=1\n")
    myfile.write("#SBATCH --time=%s\n" % (time))
    myfile.write("#SBATCH --partition=normal\n")
    myfile.write("#SBATCH --constraint=gpu\n")
    myfile.write("#SBATCH --output=batch_out/%s.%%j.o\n" % (problem_name_full))
    myfile.write("#SBATCH --error=batch_out/%s.%%j.e\n" % (problem_name_full))
    myfile.write("\n## load modules\n")
    myfile.write("source %s/util/%s\n" % (library_path,module_name))
    myfile.write("\n## set number of threads\n")
    myfile.write("export OMP_NUM_THREADS=%d\n" % (nthreads))
    myfile.write("\n## run the job\n")
    for exec_name in exec_list:
        myfile.write("%s\n" %(exec_name))
    return

def commit_batch(batchfile_list, additional_parameters):
    "this function commits batch files"
    # say what we are doing now:
    print "Commiting batch scripts: "
	# send every batch file from the list
    for batchfile_name in batchfile_list:
        print  " - %s %s" % (batchfile_name, additional_parameters);
#        call(["sbatch",additional_parameters, batchfile_name])
        call(["sbatch", additional_parameters, batchfile_name])
    return

def show_jobs(account):
    "this function shows queue"
    print "--------------------------------------- MY JOBS: -----------------------------------"
    call(["squeue", "--account=%s" % (account)])
    return

def delete_directory_content(folder):
    "this function deletes recursively the content of folder"
    print "- delete content of folder : %s" % (folder)
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)

def delete_file(file_path):
    "this function deletes the file with given path"
    print "- delete file              : %s" % (file_path)
    try:
        if os.path.isfile(file_path):
           os.unlink(file_path)
    except Exception as e:
        print(e)


