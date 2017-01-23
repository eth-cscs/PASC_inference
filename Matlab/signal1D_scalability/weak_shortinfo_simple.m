% load data from provided shortinfo files and plot strong scalability
% graphs
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

%sampleid = {'10e5', '10e6','10e7'};
sampleid = {'10e6'};

for k=1:length(sampleid)
    filename_cpu=['cpu/weak_shortinfo_final_' sampleid{k} '.txt'];
    M_cpu = csvread(filename_cpu,1,0);

    filename_cput=['cput/weak_shortinfo_final_' sampleid{k} '.txt'];
    M_cput = csvread(filename_cput,1,0);
    
    filename_gpu=['gpu/weak_shortinfo_final_' sampleid{k} '.txt'];
    M_gpu = csvread(filename_gpu,1,0);
    
    %0col ngpus,K,epssqr_best,abserr_best,Theta0,
    %5col Theta1,it all, t all, t gamma update, t gamma solve,
    %10col t theta update, t theta solve, SPGQP it, SPGQP hessmult, SPGQP t all,
    %15col SPGQP t project, SPGQP t matmult, SPGQP t dot, SPGQP t update, SPGQP t stepsize, 
    %20col SPGQP t fs, SPGQP t other, SPGQP fx, SPGQP fx_linear, SPGQP fx_quadratic, 
    %25col SPGQP_sum it, SPGQP_sum hessmult, SPGQP_sum t all, SPGQP_sum t project, SPGQP_sum t matmult, 
    %30col SPGQP_sum t dot, SPGQP_sum t update, SPGQP_sum t stepsize, SPGQP_sum t fs, SPGQP_sum t other, 

    cpu_nmb = M_cpu(:,1)';
    [cpu_nmb,cpu_sortidx] = sort(cpu_nmb);
    cpu_it = M_cpu(cpu_sortidx,26)';
    %all,projection, matmult, dot, update, fs
    cpu_times(k,:) = M_cpu(cpu_sortidx,28)';
    cpu_times_relative(k,:) = cpu_times(k,:)./cpu_it;

    cput_nmb = M_cput(:,1)';
    [cput_nmb,cput_sortidx] = sort(cput_nmb);
    cput_it = M_cput(cput_sortidx,26)';
    %all,projection, matmult, dot, update, fs
    cput_times(k,:) = M_cput(cput_sortidx,28)';
    cput_times_relative(k,:) = cput_times(k,:)./cput_it;
    
    gpu_nmb = M_gpu(:,1)';
    [gpu_nmb,gpu_sortidx] = sort(gpu_nmb);
    gpu_it = M_gpu(gpu_sortidx,26)';
    %all,projection, matmult, dot, update, fs
    gpu_times(k,:) = M_gpu(gpu_sortidx,28)';
    gpu_times_relative(k,:) = gpu_times(k,:)./gpu_it;

end

figure
hold on
set(gca,'fontsize',12);

cpu_plot = plot(1:size(cpu_times,2), cpu_times_relative(1,:), 'b', 'LineWidth', 2.0);
plot(1:size(cpu_times,2), cpu_times_relative(1,:), 'bo', 'LineWidth', 2.0);

cput_plot = plot(1:size(cput_times,2), cput_times_relative(1,:), 'Color', [0.0 0.6 0.0], 'LineWidth', 2.0);
plot(1:size(cput_times,2), cput_times_relative(1,:), 'Color', [0.0 0.6 0.0], 'Marker', 'o', 'LineWidth', 2.0);

gpu_plot = plot(1:size(cpu_times,2), gpu_times_relative(1,:), 'r', 'LineWidth', 2.0);
plot(1:size(cpu_times,2), gpu_times_relative(1,:), 'ro', 'LineWidth', 2.0);

%plot(gpu_nmb,1:size(gpu_times,2),'k--')
xlabel('nmb of nodes', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('computing time [s]', 'Interpreter', 'latex', 'FontSize', 12)

h = legend([cpu_plot cput_plot gpu_plot],'CPU 1core','CPU 24cores','GPU');
set(h,'Interpreter','latex');
set(h, 'FontSize', 16);

set(gca,'YScale','log');

hold off


figure
hold on
set(gca,'fontsize',12);

for k=1:length(sampleid)
    mybars = zeros(2,max([size(gpu_times,2),size(gpu_times,2)]));
    mybars(1,1:size(cpu_times,2)) = cpu_times_relative(k,1)./cpu_times_relative(k,:);
    mybars(2,1:size(cput_times,2)) = cput_times_relative(k,1)./cput_times_relative(k,:);
    mybars(3,1:size(gpu_times,2)) = gpu_times_relative(k,1)./gpu_times_relative(k,:);

%    cpu_plot = plot3(1:size(mybars,2),k*ones(size(mybars(1,:))),mybars(1,:),'b');
%    plot3(1:size(mybars,2),k*ones(size(mybars(1,:))),mybars(1,:),'bo');
    
%    gpu_plot = plot3(1:size(mybars,2),k*ones(size(mybars(2,:))),mybars(2,:),'r');
%    plot3(1:size(mybars,2),k*ones(size(mybars(2,:))),mybars(2,:),'ro');

    for i=1:size(mybars,2)
%        cpu_plot = fill3([i+0 i+0.2 i+0.2 i+0], [k k k k], [0 0 mybars(1,i) mybars(1,i)], 'b');
%        gpu_plot = fill3([i+0.21 i+0.41 i+0.41 i+0.21], [k k k k], [0 0 mybars(2,i) mybars(2,i)], 'r');
        cpu_plot = fill([i+0 i+0.2 i+0.2 i+0], [0 0 mybars(1,i) mybars(1,i)], 'b');
        cput_plot = fill([i+0.22 i+0.42 i+0.42 i+0.22], [0 0 mybars(2,i) mybars(2,i)], 'g');
        gpu_plot = fill([i+0.44 i+0.64 i+0.64 i+0.44], [0 0 mybars(3,i) mybars(3,i)], 'r');
        
        set(cput_plot,'FaceColor', [0.0 0.6 0.0]);
    end
end

opt_plot = plot(gpu_nmb,ones(length(gpu_times),1),'k--', 'LineWidth', 2.0);

%plot(gpu_nmb,1:size(gpu_times,2),'k--')
xlabel('number of nodes', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('efficiency', 'Interpreter', 'latex', 'FontSize', 12)
%zlabel('computing time')


h = legend([cpu_plot cput_plot gpu_plot opt_plot],'CPU 1core','CPU 24cores','GPU','optimal');
set(h,'Interpreter','latex');
set(h, 'FontSize', 16);


hold off


