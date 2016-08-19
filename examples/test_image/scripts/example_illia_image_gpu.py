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
epssqrs = [1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 1e-4, 5e-4, 1e-3, 1e-2];
#epssqrs = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2];

# how many clusters to use (value of K)
Ks = [2];

# on how many nodes to compute
Ns = [1];

# name of problem (name of compiled program)
problem_name = "test_image";

# common parameters
problem_parameters = "--test_Theta=0.4 --test_Theta=0.6 --test_cutdata=true --test_scaledata=false --test_annealing=10 --tssolver_debugmode=2 --spgqpsolver_maxit=1000 --tssolver_maxit=100 --spgqpsolver_debugmode=0 --test_shortinfo=true";

# the upper estimation of computing time
problem_time = "00:20:00"; 

# machine parameters
architecture = "GPU1";
Ngpu = 1;
Nthreads = 1;

# prepare batchscripts
print "Preparing batch scripts: %s (Nthreads=%d, Ngpu=%d)" % (architecture,Nthreads,Ngpu)
batchfile_list = [];
for dimension in dimensions:
    for epssqr in epssqrs:
        for K in Ks:
            for N in Ns:
                image_path = "%s/%s_%s_%s.bin" % (image_dir,image_name,dimension[0],dimension[1]);
                problem_name_full = "%s_%s_w%s_h%s_epssqr%.10f_K%s_arch%s_N%s_Nthreads%s_Ngpu%s" % (problem_name,image_name,dimension[0],dimension[1],epssqr,K,architecture,N,Nthreads,Ngpu)
                print " - %s: %s" % (problem_name, problem_name_full);
                problem_parameters_full = "%s --test_image_filename=\"%s\" --test_image_out=\"%s\" --test_width=%s --test_height=%s --test_epssqr=%.10f --test_K=%s --test_shortinfo_header='image_name,width,height,epssqr,K,architecture,N,Nthreads,Ngpu,' --test_shortinfo_values='%s,%d,%d,%f,%d,%s,%d,%d,%d,' --test_shortinfo_filename='shortinfo/%s.txt'" % (problem_parameters, image_path, problem_name_full, dimension[0], dimension[1], epssqr, K, image_name, dimension[0], dimension[1], epssqr, K, architecture, N, Nthreads, Ngpu, problem_name_full);
                batchfile_name = write_batchfile(problem_name, problem_name_full, problem_time, problem_parameters_full, library_path, architecture, N, Nthreads, Ngpu);
                batchfile_list.append(batchfile_name);

# run bash scripts
commit_batchfiles(batchfile_list, "c11", "normal")

# show the queue
show_jobs("c11")

