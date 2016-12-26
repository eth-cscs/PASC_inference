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
    print 'daint_batch_strong.py <inputfile>'
    sys.exit()

inputfile = sys.argv[1];
spgqpsolver_eps = "1e-4";
problem_time = "00:20:00";

# path to exec folder
username = "pospisil"
library_path = "~/soft/PASC_inference/";
build_path = "%s/" % (os.getenv( "SCRATCH"));
exec_name = "./test_signal1D"
mpiexec = "srun"
N = [1,2,3,4,5,6,7,8];

# define console parameters
params_list = [];
params_list.append("--test_filename=data/%s_data.bin" %(inputfile))
params_list.append("--test_filename_solution=data/%s_solution.bin" %(inputfile))
params_list.append("--test_filename_gamma0=data/%s_gamma0.bin" %(inputfile))
params_list.append("--test_cutdata=false --test_scaledata=false")
params_list.append("--test_epssqr=1e-2 --test_annealing=1")
params_list.append("--tssolver_maxit=1 --tssolver_debugmode=0")
params_list.append("--spgqpsolver_maxit=10000 --spgqpsolver_debugmode=0 --spgqpsolver_stop_difff=false --spgqpsolver_stop_normgp=true")
params_list.append("--test_shortinfo=true")
params_list.append("--test_K=2 --test_Theta=1.0 --test_Theta=2.0")
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
    params2 = "--test_filename_out=%s --test_shortinfo_header=ngpus, --test_shortinfo_values=%d, --test_shortinfo_filename=shortinfo/%s.txt --spgqpsolver_eps=%s" % (problem_name, N[index], problem_name, spgqpsolver_eps)
    exec_name_full = "%s -n %d %s %s %s > batch_out/%s.log" %(mpiexec, N[index], exec_name, params, params2, problem_name)
    batch_filename = os.path.join(gpu_batch_path, "%s.batch" % (problem_name))
    write_batch(problem_name, N[index], 1, 1, N[index], 1, problem_time, library_path, gpu_batch_path, exec_name_full)
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
    params2 = "--test_filename_out=%s --test_shortinfo_header=ncpus, --test_shortinfo_values=%d, --test_shortinfo_filename=shortinfo/%s.txt --spgqpsolver_eps=%s" % (problem_name, N[index], problem_name, spgqpsolver_eps)
    exec_name_full = "%s -n %d %s %s %s > batch_out/%s.log" %(mpiexec, N[index], exec_name, params, params2, problem_name)
    batch_filename = os.path.join(cpu_batch_path, "%s.batch" % (problem_name))
    write_batch(problem_name, N[index], 1, 1, N[index], 1, problem_time, library_path, cpu_batch_path, exec_name_full)
    batchfile_list.append(batch_filename);


# run bash scripts
#commit_batch(batchfile_list,"-C gpu --account=c11 --constraint=gpu")

# show the queue
#show_jobs(username)

