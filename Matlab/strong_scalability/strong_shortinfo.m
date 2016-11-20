% load data from provided shortinfo files and plot strong scalability
% graphs
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

%sampleid = {'10e4', '10e5', '10e6'};
sampleid = {'10e7'};

for k=1:length(sampleid)
    filename_cpu=['cpu/shortinfo_final_' sampleid{k} '.txt'];
    M_cpu = csvread(filename_cpu,1,0);

    filename_gpu=['gpu/shortinfo_final_' sampleid{k} '.txt'];
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
    cpu_times(1,:) = M_cpu(cpu_sortidx,28)';
    cpu_times(2,:) = M_cpu(cpu_sortidx,29)';
    cpu_times(3,:) = M_cpu(cpu_sortidx,30)';
    cpu_times(4,:) = M_cpu(cpu_sortidx,31)';
    cpu_times(5,:) = M_cpu(cpu_sortidx,32)';
    cpu_times(6,:) = M_cpu(cpu_sortidx,34)';
   
    
    gpu_nmb = M_gpu(:,1)';
    [gpu_nmb,gpu_sortidx] = sort(gpu_nmb);
        
    gpu_it = M_gpu(gpu_sortidx,26)';
    %all,projection, matmult, dot, update, fs
    gpu_times(1,:) = M_gpu(gpu_sortidx,28)';
    gpu_times(2,:) = M_gpu(gpu_sortidx,29)';
    gpu_times(3,:) = M_gpu(gpu_sortidx,30)';
    gpu_times(4,:) = M_gpu(gpu_sortidx,31)';
    gpu_times(5,:) = M_gpu(gpu_sortidx,32)';
    gpu_times(6,:) = M_gpu(gpu_sortidx,34)';

        
    titles = cell(6);
    titles{1} = 'all';
    titles{2} = 'projection';
    titles{3} = 'matmult';
    titles{4} = 'dot';
    titles{5} = 'update';
    titles{6} = 'fs';

    for i = 1:6
        cpu_times_relative(i,:) = cpu_times(i,:)./cpu_it;
        gpu_times_relative(i,:) = gpu_times(i,:)./gpu_it;
    end

    figure
    for i = 1:6
        subplot(3,6,i);
        hold on
        plot(cpu_nmb,cpu_times(i,:),'b-o')
        plot(cpu_nmb,cpu_times(i,1)./(1:size(cpu_times,2)),'b--')

        plot(gpu_nmb,gpu_times(i,:),'r-o')
        plot(gpu_nmb,gpu_times(i,1)./(1:size(gpu_times,2)),'r--')
        
        xlabel('nmb of CPUs/GPUs')
        ylabel('computation time [s]')
        ax = legend('CPU','optimal CPU','GPU','optimal GPU');
        leg = findobj(ax,'type','text');
        set(leg,'FontSize',8)
        title([titles{i} ''])
        hold off

        subplot(3,6,6+i);
        hold on
        plot(cpu_nmb,cpu_times_relative(i,:),'b-o')
        plot(cpu_nmb,cpu_times_relative(i,1)./(1:size(cpu_times,2)),'b--')
        
        plot(gpu_nmb,gpu_times_relative(i,:),'r-o')
        plot(gpu_nmb,gpu_times_relative(i,1)./(1:size(gpu_times,2)),'r--')
        xlabel('nmb of CPUs/GPUs')
        ylabel('one iteration computation time [s]')
        ax = legend('CPU','optimal CPU','GPU','optimal GPU');
        leg = findobj(ax,'type','text');
        set(leg,'FontSize',8)
        title([titles{i} ': relative'])
        hold off

        subplot(3,6,12+i);
        hold on
        mybars = zeros(2,max([size(gpu_times,2),size(gpu_times,2)]));
        mybars(1,1:size(cpu_times,2)) = cpu_times_relative(i,1)./cpu_times_relative(i,:);
        mybars(2,1:size(gpu_times,2)) = gpu_times_relative(i,1)./gpu_times_relative(i,:);

        b = bar(mybars');
        plot(gpu_nmb,1:size(gpu_times,2),'k--')
        xlabel('nmb of CPUs/GPUs')
        ylabel('speed up')
        ax = legend('CPU','GPU','optimal');
        leg = findobj(ax,'type','text');
        set(leg,'FontSize',8)
        title([titles{i} ': relative speed up'])
        hold off

    end
    
end


