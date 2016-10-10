#
# 
#

# include some funny functions
import os, shutil
from subprocess import call

# define function for writing a fun into file
def write_batchfile(problem_name, problem_name_full, problem_time, problem_parameters, library_path, architecture, N, Nthreads, Ngpu):
    "this function prints a fun into batch script file, the fun is based on parameters"
    batchfile_name = "./batch/%s.batch" % (problem_name_full);
    myfile = open(batchfile_name, 'w+');
    # write some funny stuff into file
    myfile.write("#!/bin/bash -l\n")
    myfile.write("\n## sbatch settings\n")
    myfile.write("#SBATCH --nodes=%d\n" % (N))
    myfile.write("#SBATCH --ntasks-per-core=1\n")
    myfile.write("#SBATCH --ntasks=%d\n" %(N))
    myfile.write("#SBATCH --gres=gpu:%d\n" % (Ngpu))
    myfile.write("#SBATCH --time=%s\n" % (problem_time))
    myfile.write("#SBATCH --partition=normal\n")
    myfile.write("#SBATCH --output=batch_out/%%j.%s.o\n" % (problem_name_full))
    myfile.write("#SBATCH --error=batch_out/%%j.%s.e\n" % (problem_name_full))
    myfile.write("\n## load modules\n")
    myfile.write("source %s/util/module_load_daint\n" % (library_path))
    myfile.write("source %s/util/set_petsc_daint\n" % (library_path))
    myfile.write("\n## set number of threads\n")
    myfile.write("export OMP_NUM_THREADS=%d\n" % (Nthreads))
    myfile.write("\n## run the job\n")
    myfile.write("srun -N %d -n %d -T 1 ./%s %s\n" % (N,N,problem_name,problem_parameters))
    return batchfile_name

def commit_batchfiles(batchfile_list, account, partition):
    "this function commits batch files"
    # say what we are doing now:
    print "Commiting batch scripts: "
	# send every batch file from the list
    for batchfile_name in batchfile_list:
        print  " - %s" % (batchfile_name);
        call(["sbatch", batchfile_name, "--account=%s" % (account), "--partition=%s" % (partition)])
    return

def show_jobs(account):
    "this function shows queue"
    print "--------------------------------------- MY JOBS: -----------------------------------"
    call(["squeue", "--account=%s" % (account)])
    return

def delete_directory_content(folder):
    "this function deletes recursively the content of folder"
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
    try:
        if os.path.isfile(file_path):
           os.unlink(file_path)
    except Exception as e:
        print(e)


