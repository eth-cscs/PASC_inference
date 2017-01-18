#
#
#

# import my stuff
from common_batch import write_batch
from common_batch import commit_batch
from common_batch import show_jobs
import os, shutil, sys, getopt

# parse input arguments
if len(sys.argv) < 2:
    print 'loremipsum_strong.py <epssqr>'
    sys.exit()

epssqr = sys.argv[1];
spgqpsolver_eps = "1e-6";
problem_time = "00:30:00";

# path to exec folder
username = "pospisil"
library_path = "~/soft/PASC_inference/";
build_path = "%s/" % (os.getenv( "SCRATCH"));
exec_name = "./test_image"
mpiexec = "srun"
N = [1,2,3,4,5,6,7,8];
Ntaskspernode = 24;

# define console parameters
params_list = [];
params_list.append("--test_image_filename=data/loremipsum.bin --test_width=1000 --test_height=1000")
params_list.append("--test_K=2 --test_Theta=0.5 --test_Theta=0.6")
params_list.append("--test_cutdata=true --test_scaledata=false")
params_list.append("--test_epssqr=%s" %(epssqr))
params_list.append("--spgqpsolver_eps=%s --spgqpsolver_monitor=true --test_annealing=1 --tssolver_maxit=1 --tssolver_debugmode=0 --spgqpsolver_maxit=10000 --spgqpsolver_debugmode=0 --spgqpsolver_stop_difff=false --spgqpsolver_stop_normgp=true" %(spgqpsolver_eps))

params_list.append("--test_shortinfo=true --test_shortinfo_header=imagename, --test_shortinfo_values=loremipsum, --test_shortinfo_filename=shortinfo/epssqr_G4.txt")

params = ' '.join(params_list);


# GPU
gpu_problem_name = "strong_G";
gpu_exec_path = "%s/build_gpu/" %(build_path);
gpu_batch_path = "%s/batch/" %(gpu_exec_path);

# GPU: generate bash scripts
batchfile_list = [];
for index in range(len(N)):
    print "GPU: Preparing batch scripts: %s/%s" % (index+1,len(N))
    problem_name = "%s%d" %(gpu_problem_name,N[index])
    exec_path = gpu_exec_path
    params2 = "--test_image_out=%s --test_shortinfo_header=ngpus, --test_shortinfo_values=%d, --test_shortinfo_filename=shortinfo/%s.txt" % (problem_name, N[index], problem_name)
    exec_name_full = "%s -n %d %s %s %s > batch_out/%s.log" %(mpiexec, N[index], exec_name, params, params2, problem_name)
    batch_filename = os.path.join(gpu_batch_path, "%s.batch" % (problem_name))
    write_batch(problem_name, N[index], Ntaskspernode, 1, problem_time, library_path, gpu_batch_path, exec_name_full)
    batchfile_list.append(batch_filename);


# CPU
cpu_problem_name = "strong_C";
cpu_exec_path = "%s/build_cpu/" %(build_path);
cpu_batch_path = "%s/batch/" %(cpu_exec_path);

# CPU: generate bash scripts
batchfile_list = [];
for index in range(len(N)):
    print "CPU: Preparing batch scripts: %s/%s" % (index+1,len(N))
    problem_name = "%s%d" %(cpu_problem_name,N[index])
    exec_path = cpu_exec_path
    params2 = "--test_image_out=%s --test_shortinfo_header=ncpus, --test_shortinfo_values=%d, --test_shortinfo_filename=shortinfo/%s.txt" % (problem_name, N[index], problem_name)
    exec_name_full = "%s -n %d %s %s %s > batch_out/%s.log" %(mpiexec, N[index], exec_name, params, params2, problem_name)
    batch_filename = os.path.join(cpu_batch_path, "%s.batch" % (problem_name))
    write_batch(problem_name, N[index], Ntaskspernode, 1, problem_time, library_path, cpu_batch_path, exec_name_full)
    batchfile_list.append(batch_filename);


# CPUT
cput_problem_name = "strong_CT";
cput_exec_path = "%s/build_cput/" %(build_path);
cput_batch_path = "%s/batch/" %(cput_exec_path);

# CPUT: generate bash scripts
batchfile_list = [];
for index in range(len(N)):
    print "CPUT: Preparing batch scripts: %s/%s" % (index+1,len(N))
    problem_name = "%s%d" %(cput_problem_name,N[index])
    exec_path = cput_exec_path
    params2 = "--test_image_out=%s --test_shortinfo_header=ngpus, --test_shortinfo_values=%d, --test_shortinfo_filename=shortinfo/%s.txt" % (problem_name, N[index], problem_name)
    exec_name_full = "%s -n %d %s %s %s > batch_out/%s.log" %(mpiexec, N[index]*Ntaskspernode, exec_name, params, params2, problem_name)
    batch_filename = os.path.join(cput_batch_path, "%s.batch" % (problem_name))
    write_batch(problem_name, N[index], Ntaskspernode, 1, problem_time, library_path, cput_batch_path, exec_name_full)
    batchfile_list.append(batch_filename);


# run bash scripts
#commit_batch(batchfile_list,"-C gpu --account=c11 --constraint=gpu")

# show the queue
#show_jobs(username)

