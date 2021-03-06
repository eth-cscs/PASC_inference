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

    filename_cpu2=['cpu2/weak_shortinfo_final_' sampleid{k} '.txt'];
    M_cpu2 = csvread(filename_cpu2,1,0);
    
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

    cpu2_nmb = M_cpu2(:,1)';
    [cpu2_nmb,cpu2_sortidx] = sort(cpu2_nmb);
    cpu2_it = M_cpu2(cpu2_sortidx,26)';
    %all,projection, matmult, dot, update, fs
    cpu2_times(k,:) = M_cpu2(cpu2_sortidx,28)';
    cpu2_times_relative(k,:) = cpu2_times(k,:)./cpu2_it;
    
    gpu_nmb = M_gpu(:,1)';
    [gpu_nmb,gpu_sortidx] = sort(gpu_nmb);
    gpu_it = M_gpu(gpu_sortidx,26)';
    %all,projection, matmult, dot, update, fs
    gpu_times(k,:) = M_gpu(gpu_sortidx,28)';
    gpu_times_relative(k,:) = gpu_times(k,:)./gpu_it;

end

figure
hold on

cpu_plot = plot(1:size(cpu_times,2), cpu_times_relative(1,:), 'b', 'LineWidth', 2.0);
plot(1:size(cpu_times,2), cpu_times_relative(1,:), 'bo', 'LineWidth', 2.0);

cpu2_plot = plot(1:size(cpu2_times,2), cpu2_times_relative(1,:), 'g', 'LineWidth', 2.0);
plot(1:size(cpu2_times,2), cpu2_times_relative(1,:), 'go', 'LineWidth', 2.0);

gpu_plot = plot(1:size(cpu_times,2), gpu_times_relative(1,:), 'r', 'LineWidth', 2.0);
plot(1:size(cpu_times,2), gpu_times_relative(1,:), 'ro', 'LineWidth', 2.0);

%plot(gpu_nmb,1:size(gpu_times,2),'k--')
xlabel('nmb of nodes')
ylabel('computing time [s]')

ax = legend([cpu_plot cpu2_plot gpu_plot],'CPU 1core','CPU 24cores','GPU');

leg = findobj(ax,'type','text');
set(leg,'FontSize',8)

set(gca,'YScale','log');

hold off


figure
hold on
for k=1:length(sampleid)
    mybars = zeros(2,max([size(gpu_times,2),size(gpu_times,2)]));
    mybars(1,1:size(cpu_times,2)) = cpu_times_relative(k,1)./cpu_times_relative(k,:);
    mybars(2,1:size(cpu2_times,2)) = cpu2_times_relative(k,1)./cpu2_times_relative(k,:);
    mybars(3,1:size(gpu_times,2)) = gpu_times_relative(k,1)./gpu_times_relative(k,:);

%    cpu_plot = plot3(1:size(mybars,2),k*ones(size(mybars(1,:))),mybars(1,:),'b');
%    plot3(1:size(mybars,2),k*ones(size(mybars(1,:))),mybars(1,:),'bo');
    
%    gpu_plot = plot3(1:size(mybars,2),k*ones(size(mybars(2,:))),mybars(2,:),'r');
%    plot3(1:size(mybars,2),k*ones(size(mybars(2,:))),mybars(2,:),'ro');

    for i=1:size(mybars,2)
%        cpu_plot = fill3([i+0 i+0.2 i+0.2 i+0], [k k k k], [0 0 mybars(1,i) mybars(1,i)], 'b');
%        gpu_plot = fill3([i+0.21 i+0.41 i+0.41 i+0.21], [k k k k], [0 0 mybars(2,i) mybars(2,i)], 'r');
        cpu_plot = fill([i+0 i+0.2 i+0.2 i+0], [0 0 mybars(1,i) mybars(1,i)], 'b');
        cpu2_plot = fill([i+0.22 i+0.42 i+0.42 i+0.22], [0 0 mybars(2,i) mybars(2,i)], 'g');
        gpu_plot = fill([i+0.44 i+0.64 i+0.64 i+0.44], [0 0 mybars(3,i) mybars(3,i)], 'r');
    end
end

%plot(gpu_nmb,1:size(gpu_times,2),'k--')
xlabel('nmb of nodes')
ylabel('speed up')
%zlabel('computing time')


ax = legend([cpu_plot cpu2_plot gpu_plot],'CPU 1core','CPU 24cores','GPU');

leg = findobj(ax,'type','text');
set(leg,'FontSize',8)
%title([titles{i} ': relative speed up'])

hold off


