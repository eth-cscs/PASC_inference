% load data from provided shortinfo files and plot strong scalability
% graphs
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

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
cpu_energy = M_cpu(cpu_sortidx,4); % /cpu_MPI_per_node!!!!
cpu_times = M_cpu(cpu_sortidx,27+K)';

%cpu_times(1) = cpu_times(1)-200;

cpu_times_relative = cpu_times./cpu_it;


gpu_nmb = M_gpu(:,1)';
[gpu_nmb,gpu_sortidx] = sort(gpu_nmb);
gpu_it = M_gpu(gpu_sortidx,25+K)';
%all,projection, matmult, dot, update, fs
gpu_energy = M_gpu(gpu_sortidx,4);
gpu_times = M_gpu(gpu_sortidx,27+K)';
gpu_times_relative = gpu_times./gpu_it;


if true
    figure
    hold on
    set(gca,'fontsize',12);

    cpu_plot = plot(cpu_nmb, cpu_it, 'b', 'LineWidth', 2.0);
    plot(cpu_nmb, cpu_it, 'bo', 'LineWidth', 2.0);

    gpu_plot = plot(gpu_nmb, gpu_it, 'r', 'LineWidth', 2.0);
    plot(gpu_nmb, gpu_it, 'ro', 'LineWidth', 2.0);

    xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('number of iterations', 'Interpreter', 'latex', 'FontSize', 12)

    h = legend([cpu_plot gpu_plot],'CPU','GPU');
    set(h,'Interpreter','latex');
    set(h, 'FontSize', 16);

%    set(gca,'XScale','log');

    hold off
end


if true
    figure
    hold on
    set(gca,'fontsize',12);

    cpu_plot = plot(cpu_nmb, cpu_times_relative(1,:), 'b', 'LineWidth', 2.0);
    plot(cpu_nmb, cpu_times_relative(1,:), 'bo', 'LineWidth', 2.0);

    gpu_plot = plot(gpu_nmb, gpu_times_relative(1,:), 'r', 'LineWidth', 2.0);
    plot(gpu_nmb, gpu_times_relative(1,:), 'ro', 'LineWidth', 2.0);

    xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('one iteration computing time [s]', 'Interpreter', 'latex', 'FontSize', 12)

    h = legend([cpu_plot gpu_plot],'CPU','GPU');
    set(h,'Interpreter','latex');
    set(h, 'FontSize', 16);

    axis([1 max([cpu_nmb, gpu_nmb]) ...
      min([cpu_times_relative(1,:), gpu_times_relative(1,:)]) max([cpu_times_relative(1,:), gpu_times_relative(1,:)])])

    set(gca,'YScale','log');
    set(gca,'XScale','log');

    hold off
end
    
    
if true
    figure
    hold on
    set(gca,'fontsize',12);

%    mybars = zeros(2,max([length(cpu_nmb),length(gpu_nmb)]));
%    mybars(1,1:size(cpu_times,2)) = cpu_times_relative(1)./cpu_times_relative;
%    mybars(2,1:size(gpu_times,2)) = gpu_times_relative(1)./gpu_times_relative;

%    for i=1:size(mybars,2)
%        cpu_plot = fill([(cpu_nmb(i))+0 (cpu_nmb(i))+0.3 (cpu_nmb(i))+0.3 (cpu_nmb(i))+0], [0 0 mybars(1,i) mybars(1,i)], 'b');
%        gpu_plot = fill([(gpu_nmb(i))+0.32 (gpu_nmb(i))+0.62 (cpu_nmb(i))+0.62 (cpu_nmb(i))+0.32], [0 0 mybars(2,i) mybars(2,i)], 'r');
%    end
    
    cpu_speed_up = cpu_times_relative(1)./cpu_times_relative;
    cpu_plot = plot(cpu_nmb, cpu_speed_up, 'b', 'LineWidth', 2.0);
    plot(cpu_nmb, cpu_times_relative(1)./cpu_times_relative, 'bo', 'LineWidth', 2.0);

    gpu_speed_up = gpu_times_relative(1)./gpu_times_relative;
    gpu_plot = plot(gpu_nmb, gpu_speed_up, 'r', 'LineWidth', 2.0);
    plot(gpu_nmb, gpu_speed_up, 'ro', 'LineWidth', 2.0);
    
    opt_plot = plot(cpu_nmb,cpu_nmb,'k--', 'LineWidth', 2.0);

    xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('speed up of one iteration', 'Interpreter', 'latex', 'FontSize', 12)

    h = legend([cpu_plot gpu_plot opt_plot],'CPU','GPU', 'optimal');
    set(h,'Interpreter','latex');
    set(h, 'FontSize', 16);

    hold off
end





if true
    figure
    hold on
    set(gca,'fontsize',12);

    iidx = 1:min(length(cpu_times_relative),length(gpu_times_relative));
    
    cgpu_ratio = cpu_times_relative(iidx)./gpu_times_relative(iidx);
    
    cgpu_plot = plot(cpu_nmb(iidx), cgpu_ratio, 'Color', [0,0.6,0], 'LineWidth', 2.0);
    plot(cpu_nmb(iidx), cgpu_ratio, 'o', 'Color', [0,0.6,0], 'LineWidth', 2.0);
    
    opt_plot = plot(gpu_nmb(iidx),ones(size(iidx,2)),'k--', 'LineWidth', 2.0);

    xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('CPU time / GPU time of one iteration', 'Interpreter', 'latex', 'FontSize', 12)

%    h = legend([cpu_plot gpu_plot opt_plot],'CPU','GPU', 'optimal');
%    set(h,'Interpreter','latex');
%    set(h, 'FontSize', 16);

    axis([min(gpu_nmb(iidx)), max(gpu_nmb(iidx)), 0, max(cgpu_ratio)+0.1])

    hold off
end

    

if true
    figure
    hold on
    set(gca,'fontsize',12);

    cpu_it_per_second = cpu_it./cpu_times;
    gpu_it_per_second = gpu_it./gpu_times;
    
    cpu_plot = plot(cpu_nmb, cpu_it_per_second, 'b', 'LineWidth', 2.0);
    plot(cpu_nmb, cpu_it_per_second, 'bo', 'LineWidth', 2.0);

    gpu_plot = plot(gpu_nmb, gpu_it_per_second, 'r', 'LineWidth', 2.0);
    plot(gpu_nmb, gpu_it_per_second, 'ro', 'LineWidth', 2.0);
    
%    opt_plot = plot(gpu_nmb,ones(size(gpu_times,2)),'k--', 'LineWidth', 2.0);

    xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('number of iterations per second', 'Interpreter', 'latex', 'FontSize', 12)

    h = legend([cpu_plot gpu_plot],'CPU','GPU');
    set(h,'Interpreter','latex');
    set(h, 'FontSize', 16);

    hold off
end



if true
    figure
    hold on
    set(gca,'fontsize',12);

    cpu_node_hours = cpu_nmb.*cpu_times_relative;
    gpu_node_hours = gpu_nmb.*gpu_times_relative;
    
    cpu_plot = plot(cpu_nmb, cpu_node_hours, 'b', 'LineWidth', 2.0);
    plot(cpu_nmb, cpu_node_hours, 'bo', 'LineWidth', 2.0);

    gpu_plot = plot(gpu_nmb, gpu_node_hours, 'r', 'LineWidth', 2.0);
    plot(gpu_nmb, gpu_node_hours, 'ro', 'LineWidth', 2.0);
    
%    opt_plot = plot(gpu_nmb,ones(size(gpu_times,2)),'k--', 'LineWidth', 2.0);

    xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('node seconds for one iteration', 'Interpreter', 'latex', 'FontSize', 12)

    h = legend([cpu_plot gpu_plot],'CPU','GPU');
    set(h,'Interpreter','latex');
    set(h, 'FontSize', 16);

    hold off
end

if true
    figure
    hold on
    set(gca,'fontsize',12);

    cpu_plot = plot(cpu_nmb, cpu_energy./cpu_times', 'b', 'LineWidth', 2.0);
    plot(cpu_nmb, cpu_energy./cpu_times', 'bo', 'LineWidth', 2.0);

    gpu_plot = plot(gpu_nmb, gpu_energy./gpu_times', 'r', 'LineWidth', 2.0);
    plot(gpu_nmb, gpu_energy./gpu_times', 'ro', 'LineWidth', 2.0);
    
%    opt_plot = plot(gpu_nmb,ones(size(gpu_times,2)),'k--', 'LineWidth', 2.0);

    xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
    ylabel('consumed energy [J] per iteration', 'Interpreter', 'latex', 'FontSize', 12)
    

    h = legend([cpu_plot gpu_plot],'CPU','GPU');
    set(h,'Interpreter','latex');
    set(h, 'FontSize', 16);

%    axis([min([cpu_nmb,gpu_nmb]) max([cpu_nmb,gpu_nmb]) 0 100])
    
    hold off
end

