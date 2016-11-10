#
#
#

# import my stuff
from common import write_pbs
from common import commit_pbs
from common import show_jobs
import os, shutil

# path to exec folder
username = "pos220"
main_folder = "/home_lustre/pos220/cuda_scalability";
batch_path = "%s/batch/" %(main_folder);
exec_name = "./sample_strong"
source_path = "%s/util/module_load_anselm" %(main_folder);
T = 100000000
nmb_of_tests = 1000
N = [1,2,3,4];

gpu_exec_path = "%s/build_gpu/" %(main_folder);
gpu_host_string = ["select=1:ncpus=1:mpiprocs=1:host=cn200,walltime=00:20:00",\
                   "select=1:ncpus=1:mpiprocs=1:host=cn200+1:ncpus=1:mpiprocs=1:host=cn201,walltime=00:20:00",\
                   "select=1:ncpus=1:mpiprocs=1:host=cn200+1:ncpus=1:mpiprocs=1:host=cn201+1:ncpus=1:mpiprocs=1:host=cn202,walltime=00:20:00",\
                   "select=1:ncpus=1:mpiprocs=1:host=cn200+1:ncpus=1:mpiprocs=1:host=cn201+1:ncpus=1:mpiprocs=1:host=cn202+1:ncpus=1:mpiprocs=1:host=cn203,walltime=00:20:00"];
gpu_problem_name = "strong_GPU";

cpu_exec_path = "%s/build_cpu/" %(main_folder);
cpu_host_string = gpu_host_string;
cpu_problem_name = "strong_CPU";

# GPU: generate bash scripts
batchfile_list = [];
for index in range(len(N)):
    print "GPU: Preparing batch scripts: %s/%s" % (index+1,len(N))
    problem_name = "%s%d" %(gpu_problem_name,N[index])
    host_string = gpu_host_string[index]
    exec_path = gpu_exec_path
    exec_name_full = "mpiexec -n %d %s %d %d" %(N[index], exec_name, T, nmb_of_tests)
    batch_filename = os.path.join(batch_path, "%s.pbs" % (problem_name))
    write_pbs(problem_name, host_string, batch_filename, exec_path, exec_name_full, source_path)
    batchfile_list.append(batch_filename);

for index in range(len(N)):
    print "CPU: Preparing batch scripts: %s/%s" % (index+1,len(N))
    problem_name = "%s%d" %(cpu_problem_name,N[index])
    host_string = cpu_host_string[index]
    exec_path = cpu_exec_path
    exec_name_full = "mpiexec -n %d %s %d %d" %(N[index], exec_name, T, nmb_of_tests)
    batch_filename = os.path.join(batch_path, "%s.pbs" % (problem_name))
    write_pbs(problem_name, host_string, batch_filename, exec_path, exec_name_full, source_path)
    batchfile_list.append(batch_filename);


# run bash scripts
#commit_pbs(batchfile_list)

# show the queue
#show_jobs(username)

