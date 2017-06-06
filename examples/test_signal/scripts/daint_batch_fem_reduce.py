#
#
#

# import my stuff
import os, shutil, sys, getopt

data_properties="--test_cutdata=false --test_scaledata=false --test_saveresult=false"
epssqrs="--test_epssqr=1e0 --test_epssqr=2e0 --test_epssqr=3e0 --test_epssqr=4e0 --test_epssqr=5e0 --test_epssqr=6e0 --test_epssqr=7e0 --test_epssqr=8e0 --test_epssqr=9e0 --test_epssqr=1e1 --test_epssqr=2e1 --test_epssqr=3e1 --test_epssqr=4e1 --test_epssqr=5e1 --test_epssqr=6e1 --test_epssqr=7e1 --test_epssqr=8e1 --test_epssqr=9e1 --test_epssqr=1e2 --test_epssqr=2e2 --test_epssqr=3e2 --test_epssqr=4e2 --test_epssqr=5e2 --test_epssqr=6e2 --test_epssqr=7e2 --test_epssqr=8e2 --test_epssqr=9e2 --test_epssqr=1e3 --test_epssqr=2e3 --test_epssqr=3e3 --test_epssqr=4e3 --test_epssqr=5e3 --test_epssqr=6e3 --test_epssqr=7e3 --test_epssqr=8e3 --test_epssqr=9e3 --test_epssqr=1e4 --test_epssqr=2e4 --test_epssqr=3e4 --test_epssqr=4e4 --test_epssqr=5e4 --test_epssqr=6e4 --test_epssqr=7e4 --test_epssqr=8e4 --test_epssqr=9e4 --test_epssqr=1e5 --test_epssqr=2e5 --test_epssqr=3e5 --test_epssqr=4e5 --test_epssqr=5e5 --test_epssqr=6e5 --test_epssqr=7e5 --test_epssqr=8e5 --test_epssqr=9e5 --test_epssqr=1e6 --test_epssqr=2e6 --test_epssqr=3e6 --test_epssqr=4e6 --test_epssqr=5e6 --test_epssqr=6e6 --test_epssqr=7e6 --test_epssqr=8e6 --test_epssqr=9e6 --test_epssqr=1e7 --test_epssqr=2e7 --test_epssqr=3e7 --test_epssqr=4e7 --test_epssqr=5e7 --test_epssqr=6e7 --test_epssqr=7e7 --test_epssqr=8e7 --test_epssqr=9e7 --test_epssqr=1e8"
solver="--test_annealing=1 --tssolver_maxit=1 --tssolver_debugmode=0 --spgqpsolver_maxit=10000 --spgqpsolver_debugmode=0 --spgqpsolver_stop_difff=false --spgqpsolver_stop_normgp=true --spgqpsolver_eps=1e-6"
general="--test_K=2 --test_Theta=1.0 --test_Theta=2.0"

fem_reduces=["", "--test_fem_type=1 --test_fem_reduce=0.1", "--test_fem_type=1 --test_fem_reduce=0.01"];
fem_reduces_names = ["fem1", "fem01", "fem001"];
fem_reduces_values = [1.0,0.1,0.01];

# path to exec folder
library_path = "~/soft/PASC_inference/";

# GPU: generate bash scripts
print "GPU: Preparing batch scripts:"

batchfile_list = [];
for idfile in range(100):
    outname = "signal1D_id%s" % (idfile+1);
    print "- preparing batch script: %s" % (outname)
    batchfile_name = "batch/%s.batch" % (outname);
    myfile = open(batchfile_name, 'w+');
    myfile.write("#!/bin/bash -l\n")
    myfile.write("\n## sbatch settings\n")
    myfile.write("#SBATCH --nodes=1\n")
    myfile.write("#SBATCH --ntasks-per-node=1\n")
    myfile.write("#SBATCH --ntasks-per-core=1\n")
    myfile.write("#SBATCH --threads-per-core=1\n")
    myfile.write("#SBATCH --time=00:15:00\n")
    myfile.write("#SBATCH --partition=normal\n")
    myfile.write("#SBATCH --output=batch_out/%s.%%j.o\n" % (outname))
    myfile.write("#SBATCH --error=batch_out/%s.%%j.e\n" % (outname))
    myfile.write("\n## load modules\n")
    myfile.write("source %s/util/module_load_daint_sandbox\n" % (library_path))
    myfile.write("\n## set number of threads\n")
    myfile.write("export OMP_NUM_THREADS=1\n")
    myfile.write("\n## run the jobs\n")
    myfile.write("\n\n")
    for idSigma in range(10):
        for fem_id in range(len(fem_reduces)):
            outname2 = "%s_idSigma%s" % (outname, idSigma+1);
            outname_full = "%s_%s" % (outname2, fem_reduces_names[fem_id]);
            params_list = [];
            params_list.append("--test_filename=data/%s.bin" %(outname2))
            params_list.append("--test_filename_solution=data/signal1D_solution.bin")
            params_list.append("--test_filename_out=%s" % (outname_full))
            params_list.append("--test_shortinfo_filename=shortinfo/%s.txt" % (outname_full))
            params_list.append("%s %s %s %s" %(epssqrs,solver,general,data_properties))
            params_list.append("%s --test_shortinfo=true --test_shortinfo_header=fem_reduce,idfile,idSigma, --test_shortinfo_values=%f,%f,%f, " %(fem_reduces[fem_id], fem_reduces_values[fem_id],idfile+1,idSigma+1))
            params = ' '.join(params_list);
            myfile.write("srun -n 1 ./test_signal1D %s > batch_out/%s.txt\n\n\n" %(params,outname_full))
        
    myfile.close()
    batchfile_list.append(batchfile_name);

filename_run = "send_batch.runme"
print "Preparing run script: %s" % (filename_run)
myfile_run = open(filename_run, 'w+');
myfile_run.write("#!/bin/bash\n\n")
for batchfile_id in range(len(batchfile_list)):
    myfile_run.write("echo %s/%s\n"%(batchfile_id,len(batchfile_list)))
    print "- adding file: %s" % (batchfile_list[batchfile_id])
    myfile_run.write("sbatch --account=c11 --constraint=gpu ")
    myfile_run.write(batchfile_list[batchfile_id])
    myfile_run.write("\n\n")

myfile_run.close()
#os.chmod(filename_run, 744)


# run bash scripts
#commit_batch(batchfile_list,"-C gpu --account=c11 --constraint=gpu")

# show the queue
#show_jobs(username)

