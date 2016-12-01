#
#
#

# import my stuff
from common_pbs import write_pbs
from common_pbs import commit_pbs
from common_pbs import show_jobs
import os, shutil, sys, getopt

# parse input arguments
if len(sys.argv) < 3:
    print 'anselm_pbs_weak.py <inputfile> <spgqpsolver_eps>'
    sys.exit()

inputfile = sys.argv[1];
spgqpsolver_eps = sys.argv[2];

# path to exec folder
username = "pos220"
main_folder = "/apps/permon/Lugano/PASC_inference/";
exec_name = "./test_signal1D"
mpiexec = "/apps/permon/petsc/petsc-3.7.4/bin/petscmpiexec"
N = [1,2,3,4];

# define console parameters
params_list = [];
params_list.append("--test_cutdata=false --test_scaledata=false")
params_list.append("--test_epssqr=1e-2 --test_annealing=1")
params_list.append("--tssolver_maxit=1 --tssolver_debugmode=0")
params_list.append("--spgqpsolver_maxit=10000 --spgqpsolver_debugmode=0 --spgqpsolver_stop_difff=false --spgqpsolver_stop_Anormgp=true")
params_list.append("--test_shortinfo=true")
params_list.append("--test_K=2 --test_Theta=1.0 --test_Theta=2.0")
params = ' '.join(params_list);

gpu_problem_name = "weak_G";
gpu_exec_path = "%s/examples/build_gpu/" %(main_folder);
gpu_batch_path = "%s/batch/" %(gpu_exec_path);
gpu_host_string = ["select=1:ncpus=16:mpiprocs=1:host=cn200,walltime=00:20:00",\
                   "select=1:ncpus=16:mpiprocs=1:host=cn200+1:ncpus=16:mpiprocs=1:host=cn201,walltime=00:20:00",\
                   "select=1:ncpus=16:mpiprocs=1:host=cn200+1:ncpus=16:mpiprocs=1:host=cn201+1:ncpus=16:mpiprocs=1:host=cn202,walltime=00:20:00",\
                   "select=1:ncpus=16:mpiprocs=1:host=cn200+1:ncpus=16:mpiprocs=1:host=cn201+1:ncpus=16:mpiprocs=1:host=cn202+1:ncpus=16:mpiprocs=1:host=cn203,walltime=00:20:00"];
gpu_modules_path = "%s/util/module_load_anselm_gpu" %(main_folder);

cpu_problem_name = "weak_C";
cpu_exec_path = "%s/examples/build_cpu/" %(main_folder);
cpu_batch_path = "%s/batch/" %(cpu_exec_path);
cpu_host_string = gpu_host_string;
cpu_modules_path = "%s/util/module_load_anselm_cpu" %(main_folder);


# GPU: generate bash scripts
batchfile_list = [];
for index in range(len(N)):
    print "GPU: Preparing batch scripts: %s/%s" % (index+1,len(N))
    problem_name = "%s%d" %(gpu_problem_name,N[index])
    host_string = gpu_host_string[index]
    exec_path = gpu_exec_path
    params2 = "--test_filename_out=%s --test_shortinfo_header=ngpus, --test_shortinfo_values=%d, --test_shortinfo_filename=shortinfo/%s.txt --spgqpsolver_eps=%s" % (problem_name, N[index], problem_name, spgqpsolver_eps);
    params3 = "--test_filename=data/%s_%d_data.bin --test_filename_solution=data/%s_%d_solution.bin --test_filename_gamma0=data/%s_%d_gamma0.bin" % (inputfile,N[index],inputfile,N[index],inputfile,N[index]);
    
    exec_name_full = "%s -n %d %s %s %s %s > batch_out/%s.log" %(mpiexec, N[index], exec_name, params, params2, params3, problem_name)
    batch_filename = os.path.join(gpu_batch_path, "%s.pbs" % (problem_name))
    write_pbs(problem_name, host_string, batch_filename, exec_path, exec_name_full, gpu_modules_path)
    batchfile_list.append(batch_filename);

# CPU: generate bash scripts
batchfile_list = [];
for index in range(len(N)):
    print "CPU: Preparing batch scripts: %s/%s" % (index+1,len(N))
    problem_name = "%s%d" %(cpu_problem_name,N[index])
    host_string = cpu_host_string[index]
    exec_path = cpu_exec_path
    params2 = "--test_filename_out=%s --test_shortinfo_header=ncpus, --test_shortinfo_values=%d, --test_shortinfo_filename=shortinfo/%s.txt --spgqpsolver_eps=%s" % (problem_name, N[index], problem_name, spgqpsolver_eps)
    exec_name_full = "%s -n %d %s %s %s > batch_out/%s.log" %(mpiexec, N[index], exec_name, params, params2, problem_name)
    batch_filename = os.path.join(cpu_batch_path, "%s.pbs" % (problem_name))
    write_pbs(problem_name, host_string, batch_filename, exec_path, exec_name_full, cpu_modules_path)
    batchfile_list.append(batch_filename);


# run bash scripts
#commit_pbs(batchfile_list)

# show the queue
#show_jobs(username)

