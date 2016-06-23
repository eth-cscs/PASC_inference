#
#
#

# import my stuff
from myscripts import write_batchfiles
from myscripts import commit_batchfiles
from myscripts import show_jobs

# times to compute
Ts = [1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000];
#Ts = [100000];

# on how many nodes to compute
Ns = [1, 2, 3, 4, 5, 6];
#Ns = [1, 2];

# name of problem
problem_name = "test_kmeans"
problem_parameters = "--test_K=3 --spgqpsolver_eps=0.0001 --spgqpsolver_stop_Anormgp_normb=true --spgqpsolver_stop_normgp_normb=true --spgqpsolver_debug_mode=3 --spgqpsolver_maxit=500 --tssolver_debug_mode=3 --spgqpsolver_dotfloor=2 --spgqpsolver_m=20 --test_epssqr=10"
library_path = "~/soft/PASC_inference"
problem_time = "00:10:00"

# generate bash scripts
write_batchfiles(Ts,Ns,problem_name, problem_time, problem_parameters, library_path, "GPU1", 1, 1)

# run bash scripts
commit_batchfiles(Ts,Ns,problem_name, "GPU1", "c11", "normal")

# show the queue
show_jobs("c11")
