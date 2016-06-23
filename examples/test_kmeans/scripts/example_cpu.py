#
#
#

# import my stuff
from myscripts import write_batchfiles
from myscripts import commit_batchfiles
from myscripts import show_jobs

# times to compute
Ts = [1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000];
#Ts = [1000, 5000];

# on how many nodes to compute
Ns = [1, 2, 3, 4, 5, 6];
#Ns = [1, 2];

# name of problem
problem_name = "test_kmeans"
problem_parameters = "--test_K=3 --spgqpsolver_eps=0.00001 --spgqpsolver_stop_Anormgp=false --spgqpsolver_stop_normgp_normb=true --spgqpsolver_debug_mode=3 --tssolver_debug_mode=3 --spgqpsolver_dotfloor=3 --spgqpsolver_m=20 --test_epssqr=10 --test_shortinfo=true --test_savevtk=false --test_savecsv=false";
library_path = "~/soft/PASC_inference"
problem_time = "00:15:00"

# generate bash scripts
write_batchfiles(Ts,Ns,problem_name, problem_time, problem_parameters, library_path, "CPU1", 0, 1)
write_batchfiles(Ts,Ns,problem_name, problem_time, problem_parameters, library_path, "CPU8", 0, 8)

# run bash scripts
commit_batchfiles(Ts,Ns,problem_name, "CPU1", "c11", "normal")
commit_batchfiles(Ts,Ns,problem_name, "CPU8", "c11", "normal")

# show the queue
show_jobs("c11")
