#
#
#

# import my stuff
from myscripts import write_batchfiles

# times to compute
Ts = [1000, 5000];

# on how many nodes to compute
Ns = [1, 2];

# name of problem
problem_name = "test_kmeans"
problem_parameters = "--test_K=3 --spgqpsolver_eps=0.00001 --spgqpsolver_stop_Anormgp=false --spgqpsolver_stop_normgp_normb=true --spgqpsolver_debug_mode=3 --tssolver_debug_mode=3 --spgqpsolver_dotfloor=2 --spgqpsolver_m=10 --test_epssqr=1 --test_shortinfo=true";
library_path = "~/soft/PASC_inference"
problem_time = "00:05:00"

# generate bash scripts
write_batchfiles(Ts,Ns,problem_name, problem_time, problem_parameters, library_path, "CPU1", 0, 1)
write_batchfiles(Ts,Ns,problem_name, problem_time, problem_parameters, library_path, "CPU8", 0, 8)
write_batchfiles(Ts,Ns,problem_name, problem_time, problem_parameters, library_path, "GPU1", 1, 1)

# run bash scripts
commit_batchfiles(Ts,Ns,problem_name, "CPU1", "c11", "normal")
commit_batchfiles(Ts,Ns,problem_name, "CPU8", "c11", "normal")
commit_batchfiles(Ts,Ns,problem_name, "GPU1", "c11", "normal")

# show the queue
show_jobs("c11")
