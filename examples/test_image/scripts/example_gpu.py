#
#
#

# import my stuff
from myscripts import write_batchfiles
from myscripts import commit_batchfiles
from myscripts import show_jobs

# path to PASC_inference library
library_path = "~/soft/PASC_inference";

# image_path = [image_dir]/[begin]_[width]_[height]_[noise].bin
image_name = "usi";
image_dir = "data/usi_text";

# dimensions of images to compute [width,height]
dimensions = [[250,150],[500,300],[1000,600]];

# noise of input images (part of the name of image)
noises = ['00','02','04','10'];

# used penalty
epssqrs = [1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 1e-4, 5e-4, 1e-3, 1e-2];

# how many clusters to use (value of K)
Ks = [2];

# on how many nodes to compute
Ns = [1];

# name of problem (name of compiled program)
problem_name = "test_image";

# common parameters
problem_parameters = "--test_cutdata=true --test_scaledata=false --test_annealing=10 --tssolver_debug_mode=2 --spgqpsolver_maxit=1000 --tssolver_maxit=100 --spgqpsolver_debug_mode=0 --test_shortinfo=true --test_Theta=0.4 --test_Theta=0.5";

# the upper estimation of computing time
problem_time = "00:10:00"; 

# generate bash scripts
batchfile_list = write_batchfiles(image_dir, image_name, dimensions, noises, epssqrs, Ks, Ns, problem_name, problem_time, problem_parameters, library_path, "GPU1", 1, 1);

# run bash scripts
commit_batchfiles(batchfile_list, "c11", "normal")

# show the queue
show_jobs("c11")

