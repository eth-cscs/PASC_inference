#
#
#

# import my stuff
from ../myscripts import write_batchfile
from ../myscripts import commit_batchfiles
from ../myscripts import show_jobs

# path to PASC_inference library
library_path = "~/soft/PASC_inference";

# file with data
data_path = "data/S001R01.edf";
data_name = "S001R01";

# how many records to take from file (-1 = all)
max_record_nmb = -1;

# file with graph coordinates in PETSc vector format
graph_path = "data/Koordinaten_EEG_P.bin";
graph_coeff = 2.5;

# used penalty
#epssqrs = [1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 1e-4, 5e-4, 1e-3, 1e-2];
epssqrs = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2];

# how many clusters to use (value of K)
Ks = [2];

# on how many nodes to compute
Ns = [1];

# name of problem (name of compiled program)
problem_name = "test_edf";

# common parameters
problem_parameters = "--test_savevtk=true --test_cutdata=true --test_scaledata=false --test_annealing=10 --tssolver_debugmode=2 --spgqpsolver_maxit=1000 --tssolver_maxit=100 --spgqpsolver_debugmode=0 --test_shortinfo=true";

# the upper estimation of computing time
problem_time = "00:20:00"; 

# machine parameters
architecture = "GPU1";
Ngpu = 1;
Nthreads = 1;

# generate bash scripts
print "Preparing batch scripts: %s (Nthreads=%d, Ngpu=%d)" % (architecture,Nthreads,Ngpu)
batchfile_list = []; # in this list store all generated batch files
for epssqr in epssqrs:
    for K in Ks:
        for N in Ns:
            problem_name_full = "%s_%s_epssqr%f_K%s_arch%s_N%s_Nthreads%s_Ngpu%s" % (problem_name,data_name,epssqr,K,architecture,N,Nthreads,Ngpu)
            print " - %s: %s" % (problem_name, problem_name_full);
            problem_parameters_full = "%s --test_data_filename=\"%s\" --test_max_record_nmb=\"%d\" --test_graph_coordinates=\"%s\" --test_graph_coeff=%f --test_data_out=\"%s\" --test_epssqr=%f --test_K=%s --test_shortinfo_header='data_name,epssqr,K,architecture,N,Nthreads,Ngpu,' --test_shortinfo_values='%s,%f,%d,%s,%d,%d,%d,' --test_shortinfo_filename='shortinfo/%s.txt'" % (problem_parameters, data_path, max_record_nmb, graph_path, graph_coeff, problem_name_full, epssqr, K, data_name, epssqr, K, architecture, N, Nthreads, Ngpu, problem_name_full);
            batchfile_name = write_batchfile(problem_name, problem_name_full, problem_time, problem_parameters_full, library_path, architecture, N, Nthreads, Ngpu);
            batchfile_list.append(batchfile_name);

# run bash scripts
commit_batchfiles(batchfile_list, "c11", "normal")

# show the queue
show_jobs("c11")

