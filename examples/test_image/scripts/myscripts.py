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
def write_batchfiles(image_dir, image_name, dimensions, noises, epssqrs, Ks, Ns, problem_name, problem_time, problem_parameters, library_path, architecture, Nthreads, Ngpu):
    "this function creates a lot of bash scripts"
    # say hello
    print "Preparing batch scripts: %s (Nthreads=%d, Ngpu=%d)" % (architecture,Nthreads,Ngpu)
    batchfile_list = [];
    for dimension in dimensions:
        for noise in noises:
            for epssqr in epssqrs:
                for K in Ks:
                    for N in Ns:
                        image_path = "%s/%s_%s_%s_%s.bin" % (image_dir,image_name,dimension[0],dimension[1],noise);
                        problem_name_full = "%s_%s_w%s_h%s_noise%s_epssqr%f_K%s_arch%s_N%s_Nthreads%s_Ngpu%s" % (problem_name,image_name,dimension[0],dimension[1],noise,epssqr,K,architecture,N,Nthreads,Ngpu)
                        print " - %s: %s" % (problem_name, problem_name_full);
                        problem_parameters_full = "%s --test_image_filename=\"%s\" --test_image_out=\"%s\" --test_width=%s --test_height=%s --test_epssqr=%f --test_K=%s " % (problem_parameters, image_path, problem_name_full, dimension[0], dimension[1], epssqr, K);
                        batchfile_name = write_batchfile(problem_name, problem_name_full, problem_time, problem_parameters_full, library_path, architecture, N, Nthreads, Ngpu);
                        batchfile_list.append(batchfile_name);
    return batchfile_list

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
    myfile.write("#SBATCH --output=%s.%%j.o\n" % (problem_name_full))
    myfile.write("#SBATCH --error=%s.%%j.e\n" % (problem_name_full))
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
    # say hello
    print "Commiting batch scripts: "
    for batchfile_name in batchfile_list:
        print  " - %s" % (batchfile_name);
        call(["sbatch", batchfile_name, "--account=%s" % (account), "--partition=%s" % (partition)])
    return

def show_jobs(account):
    "this function shows queue"
    print "--------------------------------------- MY JOBS: -----------------------------------"
    call(["squeue", "--account=%s" % (account)])
    return
    	
##sbatch scripts/test.batch --account=c11 --partition=normal
##call(["ls", "-l"])

