#
# 
#

# include some funny functions
import os

from subprocess import call

# create folder
if not os.path.exists("batch"):
    os.makedirs("batch")

if not os.path.exists("shortinfo"):
    os.makedirs("shortinfo")

# define function for writing a lot of batchscripts
def write_batchfiles(Ts,Ns,problem_name, problem_time, problem_parameters, library_path, architecture, n_gpu, n_threads):
    "this function creates a lot of bash scripts"
    # say hello
    print "Preparing batch scripts: %s (%d threads, %d gpus)" % (architecture,n_threads,n_gpu)
    for t in Ts:
        for n in Ns:
            print " - N = %d, T = %d" % (n,t);
            write_batchfile(problem_name, problem_time, problem_parameters, library_path, architecture, n_gpu*n, n_threads, n, t)
    return

# define function for writing a fun into file
def write_batchfile(problem_name, problem_time, problem_parameters, library_path, architecture, n_gpu, n_threads, n, t):
    "this function prints a fun into batch script file, the fun is based on parameters"
    myfile = open("./batch/%s_%s_N%dT%d.batch" % (problem_name,architecture,n,t), 'w+');
    # write some funny stuff into file
    myfile.write("#!/bin/bash -l\n")
    myfile.write("\n## sbatch settings\n")
    myfile.write("#SBATCH --nodes=%d\n" % (n))
    myfile.write("#SBATCH --ntasks-per-core=1\n")
    myfile.write("#SBATCH --ntasks=%d\n" %(n))
    myfile.write("#SBATCH --gres=gpu:%d\n" % (n_gpu))
    myfile.write("#SBATCH --time=%s\n" % (problem_time))
    myfile.write("#SBATCH --partition=normal\n")
    myfile.write("#SBATCH --output=%s.%%j.o\n" % (problem_name))
    myfile.write("#SBATCH --error=%s.%%j.e\n" % (problem_name))
    myfile.write("\n## load modules\n")
    myfile.write("source %s/util/module_load_daint\n" % (library_path))
    myfile.write("source %s/util/set_petsc_daint\n" % (library_path))
    myfile.write("\n## set number of threads\n")
    myfile.write("export OMP_NUM_THREADS=%d\n" % (n_threads))
    myfile.write("\n## run the job\n")
    myfile.write("srun -N %d -n %d -T 1 ./test_kmeans_load" % (n,n))
    myfile.write(" --test_data_filename='data/data_kmeans_T%d.bin' --test_gamma0_filename='data/gamma0_kmeans_T%dK3.bin'" % (t,t))
    myfile.write(" --test_shortinfo_header='architecture,nodes,threads,T,K,' --test_shortinfo_values='%s,%d,%d,%d,3,' --test_shortinfo_filename='shortinfo/kmeans_%s_N%d_threads%d_T%d_K3.txt'" % (architecture,n,n_threads,t,architecture,n,n_threads,t))
    myfile.write(" %s\n" % (problem_parameters))
    return

def commit_batchfiles(Ts,Ns,problem_name, architecture, account, partition):
    "this function commits batch files"
    # say hello
    print "Commiting batch scripts: %s" % (architecture)
    for t in Ts:
        for n in Ns:
            filename = "batch/%s_%s_N%dT%d.batch" % (problem_name,architecture,n,t)
            print  " - %s" % (filename);
            call(["sbatch", filename, "--account=%s" % (account), "--partition=%s" % (partition)])
    return

def show_jobs(account):
    "this function shows queue"
    print "--------------------------------------- MY JOBS: -----------------------------------"
    call(["squeue", "--account=%s" % (account)])
    return
    	
##sbatch scripts/test.batch --account=c11 --partition=normal
##call(["ls", "-l"])

