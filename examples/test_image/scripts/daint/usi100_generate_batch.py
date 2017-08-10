#
#
#

# import my stuff
from common_batch import write_batch
from common_batch import commit_batch
from common_batch import show_jobs
import os, shutil, sys, getopt

# path to PASC_inference library
library_path = "~/soft/PASC_inference/";
problem_name="usi";
exec_name = "./test_image";
mpiexec = "srun"
build_path = "%s/build_gpu/" % (os.getenv( "SCRATCH"));
batch_path = "batch/%s/" %(problem_name);
module_name = "module_load_daint_sandbox";

width=1000;
height=600; 
xdim=1;
filename_solution = "data/%s/solution.bin" % (problem_name);
data_type=0; # type of input data [0=TRn, 1=TnR, 2=nTR]
K=2;
mu0="0.498039215686275";
mu1="0.501960784313725";
annealing=1;
fem_reduce=1.0;
fem_type=1;
matrix_type=1;

N = 1;
problem_time = "00:20:00";

# noise of input signal
nmbfilesmax = 2;
Sigma = [0,1,2,3,4,5,6,7,8,9];
#nmbfilesmax = 100;
#Sigma = [1,2,3,4,5,6,7,8,9,10];

# used penalty
#epssqrs = [1e-7, 1e-6, 1e2];
epssqrs = [1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1e0, 5e0, 1e1, 5e1, 1e2, 5e2, 1e3];

# define console parameters
params_list = [];
params_list.append("--test_filename_solution=%s --test_data_type=%s --test_width=%s --test_height=%s --test_xdim=%s" %(filename_solution, data_type, width,height,xdim) );
params_list.append("--test_K=%s --test_annealing=%s --test_printstats=false --test_printinfo=true" % (K, annealing) );
params_list.append("--test_Theta=%s --test_Theta=%s --tssolver_thetasolver_updatebeforesolve=false" % (mu0,mu1) );
params_list.append("--test_cutdata=false --test_scaledata=false" );
params_list.append("--tssolver_maxit=1 --tssolver_eps=1e-5 --tssolver_debugmode=0" );
params_list.append("--log_or_not=false --log_or_not_func_call=true --log_or_not_file_line=true --log_or_not_level=true" );
params_list.append("--spgqpsolver_maxit=5000 --spgqpsolver_debugmode=0 --spgqpsolver_stop_difff=false --spgqpsolver_stop_normgp=true --spgqpsolver_eps=1e-5" );
params_list.append("--spgqpsolver_debug_print_it=false --tssolver_debug_print_gamma=false --tssolver_debug_print_theta=false --tssolver_debug_print_it=false" );
params_list.append("--test_fem_reduce=%s --test_fem_type=%s" % (fem_reduce, fem_type) );
params_list.append("--graphh1femmodel_matrixtype=%s" % (matrix_type) );
params = ' '.join(params_list);

# add penalties to problem parameters
for epssqr in epssqrs:
    params = "--test_epssqr=%.16f %s" % (epssqr,params)

# create directory for batch files
if not os.path.exists(batch_path):
    os.makedirs(batch_path)

# generate bash scripts
batchfile_list = [];
for n in range(0,nmbfilesmax):
    print "Preparing batch scripts: %s" % (n)
    exec_list = [];
    for sigma in Sigma:
        filename_in = "data/%s/%s_id%s_idSigma%s.bin" % (problem_name,problem_name,n,sigma);
        filename_out = "%s_id%s_idSigma%s" % (problem_name,n,sigma);
        shortinfo_header = "n,sigmaid,";
        shortinfo_values = "%s,%s," % (n,sigma);
        shortinfo_filename = "shortinfo/%s_id%s_idSigma%d.txt" % (problem_name,n,sigma);
        params2 = "--test_filename_in=\"%s\" --test_filename_out=\"%s\" --test_shortinfo_header=\"%s\" --test_shortinfo_values=\"%s\" --test_shortinfo_filename=\"%s\"" % (filename_in, filename_out, shortinfo_header, shortinfo_values, shortinfo_filename);
        exec_name_full = "%s -n %d %s %s %s > batch_out/%s.log" %(mpiexec, N, exec_name, params, params2, filename_out)
        exec_list.append(exec_name_full)        
    batch_name = "%s_id%s" % (problem_name, n);
    batch_filename = os.path.join(batch_path, "%s.batch" % (batch_name))
    write_batch(batch_name, N, 1, 1, problem_time, library_path, batch_path, exec_list, module_name)
    batchfile_list.append(batch_filename);

# run bash scripts
commit_batch(batchfile_list, "--account=s747")

# show the queue
#show_jobs("s747")

