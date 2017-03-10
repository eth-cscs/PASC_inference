% load data from provided shortinfo files and plot strong scalability
% graphs
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath('common')

K = 2;
T = 10^7;
cpu_MPI_per_node = 36;

filename_cpu=['results/shortinfo_strong_cpu.txt'];
M_cpu = csvread(filename_cpu,1,0);
    
filename_gpu=['results/shortinfo_strong_gpu.txt'];
M_gpu = csvread(filename_gpu,1,0);
    
    %0col ngpus,K,epssqr_best,abserr_best,energy,
    %5col Theta0,Theta1,it all, t all, t gamma update, 
    %10col t gamma solve, t theta update, t theta solve, SPGQP it, SPGQP hessmult, 
    %15col SPGQP t all,SPGQP t project, SPGQP t matmult, SPGQP t dot, SPGQP t update,  
    %20col SPGQP t stepsize, SPGQP t fs, SPGQP t other, SPGQP fx, SPGQP fx_linear,  
    %25col SPGQP fx_quadratic, SPGQP_sum it, SPGQP_sum hessmult, SPGQP_sum t all, SPGQP_sum t project, 
    %30col SPGQP_sum t matmult, SPGQP_sum t dot, SPGQP_sum t update, SPGQP_sum t stepsize, SPGQP_sum t fs, 
    %31col SPGQP_sum t other,

cpu_nmb = M_cpu(:,1)';
[cpu_nmb,cpu_sortidx] = sort(cpu_nmb);
cpu_it = M_cpu(cpu_sortidx,25+K)';
%all,projection, matmult, dot, update, fs
cpu_times_all = M_cpu(cpu_sortidx,27+K)'./cpu_it;
cpu_times_project = M_cpu(cpu_sortidx,28+K)'./cpu_it;
cpu_times_matmult = M_cpu(cpu_sortidx,29+K)'./cpu_it;
cpu_times_dot = M_cpu(cpu_sortidx,30+K)'./cpu_it;
cpu_times_update = M_cpu(cpu_sortidx,31+K)'./cpu_it;
cpu_times_stepsize = M_cpu(cpu_sortidx,32+K)'./cpu_it;
cpu_times_fs = M_cpu(cpu_sortidx,33+K)'./cpu_it;
cpu_times_other = M_cpu(cpu_sortidx,34+K)'./cpu_it;


gpu_nmb = M_gpu(:,1)';
[gpu_nmb,gpu_sortidx] = sort(gpu_nmb);
gpu_it = M_gpu(gpu_sortidx,25+K)';
%all,projection, matmult, dot, update, fs
gpu_times_all = M_gpu(gpu_sortidx,27+K)'./gpu_it;
gpu_times_project = M_gpu(gpu_sortidx,28+K)'./gpu_it;
gpu_times_matmult = M_gpu(gpu_sortidx,29+K)'./gpu_it;
gpu_times_dot = M_gpu(gpu_sortidx,30+K)'./gpu_it;
gpu_times_update = M_gpu(gpu_sortidx,31+K)'./gpu_it;
gpu_times_stepsize = M_gpu(gpu_sortidx,32+K)'./gpu_it;
gpu_times_fs = M_gpu(gpu_sortidx,33+K)'./gpu_it;
gpu_times_other = M_gpu(gpu_sortidx,34+K)'./gpu_it;


gpu_times = M_gpu(gpu_sortidx,27+K)';
gpu_times_relative = gpu_times./gpu_it;



if true
    figure

    subplot(4,2,1)
    plot_me( 'all', cpu_nmb, cpu_times_all, gpu_nmb, gpu_times_all )
    hold off

    subplot(4,2,2)
    plot_me( 'projection', cpu_nmb, cpu_times_project, gpu_nmb, gpu_times_project )

    subplot(4,2,3)
    plot_me( 'matmult', cpu_nmb, cpu_times_matmult, gpu_nmb, gpu_times_matmult )
    
    subplot(4,2,4)
    plot_me( 'dot', cpu_nmb, cpu_times_dot, gpu_nmb, gpu_times_dot )

    subplot(4,2,5)
    plot_me( 'update', cpu_nmb, cpu_times_update, gpu_nmb, gpu_times_update )

    subplot(4,2,6)
    plot_me( 'stepsize', cpu_nmb, cpu_times_stepsize, gpu_nmb, gpu_times_stepsize )

    subplot(4,2,7)
    plot_me( 'fs', cpu_nmb, cpu_times_fs, gpu_nmb, gpu_times_fs )

    subplot(4,2,8)
    plot_me( 'other', cpu_nmb, cpu_times_other, gpu_nmb, gpu_times_other )
    
    
    hold off
end
    
    
