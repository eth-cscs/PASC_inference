#
#
#

# import my stuff
from myscripts import write_batchfile
from myscripts import commit_batchfiles
from myscripts import show_jobs

# path to PASC_inference library
library_path = "~/soft/PASC_inference";

# image_path = [image_dir]/[begin]_[width]_[height].bin
image_name = "C_noise_medium";
#image_name = "C_noise_large";

image_dir = "data/illia_image";

# dimensions of images to compute [width,height]
#dimensions = [[640,202]];
dimensions = [[1683,374]];

# used penalty
#epssqrs = [1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 1e-4, 5e-4, 1e-3, 1e-2];
epssqrs = [0, 1e-14, 1e-12, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];

# how many clusters to use (value of K)
Ks = [2];

# on how many nodes to compute
Ns = [1];

# name of problem (name of compiled program)
problem_name = "test_image";

# prepare epssqr parameters
epssqr_parameters = "";
for epssqr in epssqrs:
    epssqr_parameters = "%s --test_epssqr=%.16f" % (epssqr_parameters, epssqr);

# common parameters
problem_parameters = "--test_Theta=0.5 --test_Theta=0.51 --test_cutdata=true --test_scaledata=false --test_annealing=10 --tssolver_debugmode=2 --spgqpsolver_maxit=5000 --tssolver_maxit=100 --spgqpsolver_debugmode=0 --test_shortinfo=true %s" % (epssqr_parameters);

# the upper estimation of computing time
problem_time = "00:40:00"; 

# machine parameters
architecture = "GPU1";
Ngpu = 1;
Nthreads = 1;

# prepare batchscripts
print "Preparing batch scripts: %s (Nthreads=%d, Ngpu=%d)" % (architecture,Nthreads,Ngpu)
batchfile_list = [];
for dimension in dimensions:
    for K in Ks:
        for N in Ns:
            image_path = "%s/%s_%s_%s.bin" % (image_dir,image_name,dimension[0],dimension[1]);
            problem_name_full = "%%j_%s_%s_w%s_h%s_K%s_arch%s_N%s_Nthreads%s_Ngpu%s" % (image_name,dimension[0],dimension[1],architecture,N,Nthreads,Ngpu)
            print " - %s: %s" % (problem_name, problem_name_full);
            problem_parameters_full = "%s --test_image_filename=\"%s\" --test_image_out=\"%s\" --test_width=%s --test_height=%s --test_K=%s --test_shortinfo_header='image_name,architecture,N,Nthreads,Ngpu,' --test_shortinfo_values='%s,%s,%d,%d,%d,' --test_shortinfo_filename='shortinfo/%s.txt'" % (problem_parameters, image_path, problem_name_full, dimension[0], dimension[1], K, image_name, architecture, N, Nthreads, Ngpu, problem_name_full);
            batchfile_name = write_batchfile(problem_name, problem_name_full, problem_time, problem_parameters_full, library_path, architecture, N, Nthreads, Ngpu);
            batchfile_list.append(batchfile_name);

# run bash scripts
commit_batchfiles(batchfile_list, "c11", "normal")

# show the queue
show_jobs("c11")

