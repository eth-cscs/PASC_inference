#
#
#

# import my stuff
from myscripts import write_batchfile
from myscripts import commit_batchfiles
from myscripts import show_jobs

# path to PASC_inference library
library_path = "~/soft/PASC_inference";

# nmbfilesmax
filename_solution = "data/signal1D/signal1D_solution.bin";
problem_name = "test_signal1D";

# noise of input signal
nmbfilesmax = 10;
Sigma = [1,2,3];
#Sigma = [1,2,3,4,5,6,7,8,9,10];


# used penalty
epssqrs = [1e-7, 1e-6, 1e-5, 1e-4];
#epssqrs = [1e-14, 1e-12, 1e-10, 1e-8, 1e-7, 2e-7, 3e-7, 4e-7, 5e-7, 6e-7, 7e-7, 8e-7, 9e-7,
#		   1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 
#		   1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 
#		   1e-4, 5e-4, 1e-3, 1e-2];

# name of problem (name of compiled program)
problem_name = "test_signal1D";

# add parameters
problem_parameters="";
for epssqr in epssqrs:
    problem_parameters = "--epssqr=%.16f %s" % (epssqr,problem_parameters)

problem_parameters = "%s --test_filename_solution=\"%s\" --test_shortinfo=true --test_cutdata=false --test_scaledata=false --test_annealing=1 --tssolver_debugmode=0 --spgqpsolver_maxit=1000 --tssolver_maxit=100 --spgqpsolver_debugmode=0 --test_shortinfo=true --test_K=2 --test_Theta=1.0 --test_Theta=2.0" % (problem_parameters,filename_solution);

# the upper estimation of computing time
problem_time = "00:20:00";

# machine parameters
architecture = "GPU1";
N = 1;
Ngpu = 1;
Nthreads = 1;

# generate bash scripts
batchfile_list = [];
for n in range(1,nmbfilesmax+1):
    print "Preparing batch scripts: %s" % (n)
    for sigma in Sigma:
        filename = "data/signal1D/signal1D_id%s_idSigma%s.bin" % (n,sigma);
        filename_out = "signal1D_id%s_idSigma_%s" % (n,sigma);
        shortinfo_header = "n,sigmaid,";
        shortinfo_values = "%s,%s," % (n,sigma);
        shortinfo_filename = "shortinfo/id%s_idSigma%d.txt" % (n,sigma);
        problem_name_full = "%s_id%s_idSigma%s" % (problem_name,n,sigma);
        problem_parameters_full = "%s --test_filename=\"%s\" --test_filename_out=\"%s\" --test_shortinfo_header=\"%s\" --test_shortinfo_values=\"%s\" --test_shortinfo_filename=\"%s\"" % (problem_parameters, filename, filename_out, shortinfo_header, shortinfo_values, shortinfo_filename);
        batchfile_name = write_batchfile(problem_name, problem_name_full, problem_time, problem_parameters_full, library_path, architecture, N, Nthreads, Ngpu);
        batchfile_list.append(batchfile_name);

# run bash scripts
commit_batchfiles(batchfile_list, "c11", "normal")

# show the queue
show_jobs("c11")

