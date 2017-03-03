% plot the results - original signal, solution and recoved signal

% include functions
addpath(genpath(fullfile(pwd,'../common')));

% filenames
if false
    filename_data = 'tests/small_data.bin';
    filename_rec = 'tests/small_recovered.bin';
    filename_orig = 'tests/small_solution.bin';
    fem_reduce = 0.13;
end
    
if true
    filename_data = 'data/strong_10e8_data.bin';
    filename_orig = 'data/strong_10e8_solution.bin';
    filename_rec = 'data/strong_10e8_recovered.bin';

    fem_reduce = 1e0;
%    filename_rec =  'tests/reduce_1e0_recovered.bin';
%    filename_rec =  'tests/reduce_1e-1_recovered.bin';
%    filename_rec =  'tests/reduce_1e-2_recovered.bin';
%    filename_rec =  'tests/reduce_1e-3_recovered.bin';
%    filename_rec =  'tests/reduce_1e-4_recovered.bin';
end

if false
    filename_data = 'data/strong_10e7_data_small_noise.bin';
    filename_orig = 'data/strong_10e7_solution_small_noise.bin';

    fem_reduce = 1e-4;
%    filename_rec =  ['tests/reduce_1e0_recovered_small_noise.bin'];
%    filename_rec =  'tests/reduce_1e-1_recovered_small_noise.bin';
%    filename_rec =  'tests/reduce_1e-2_recovered_small_noise.bin';
%    filename_rec =  'tests/reduce_1e-3_recovered_small_noise.bin';
    filename_rec =  'tests/reduce_1e-4_recovered_small_noise.bin';
end


% load vectors
vec_data = loadbin(filename_data);
vec_orig = loadbin(filename_orig);
vec_rec = loadbin(filename_rec);

% what to plot
step = 10000;
t_idx1 = step:2*step;%length(vec_data);
t_idx2 = step:2*step;%length(vec_data);

disp(['norm1 = ' num2str(norm(vec_orig-vec_rec,1))])
disp(['norm2 = ' num2str(norm(vec_orig-vec_rec,2))])

if false
   for t=1:length(vec_rec)
      if vec_rec(t) < 1.5
          vec_rec(t) = 1;
      else
          vec_rec(t) = 2;
      end
   end
end

if false
figure
hold on
set(gca,'fontsize',12);

plot(t_idx,vec_data(t_idx),'b')
plot(t_idx,vec_orig(t_idx),'r')
%plot(t_idx,vec_rec(t_idx),'m')

%axis([t_idx(1) t_idx(end) 0 3])

%h = legend('with noise', 'exact', 'recovered');
h = legend('with noise', 'exact');
set(h,'Interpreter','latex');
set(h, 'FontSize', 16);

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$x(t)$', 'Interpreter', 'latex', 'FontSize', 12)

%sizey = 10;
%plot([60000,90000,90000,60000,60000],[1-sizey,1-sizey,2+sizey,2+sizey,1-sizey],'k')

hold off
end


figure
hold on
set(gca,'fontsize',12);

plot(t_idx2,vec_data(t_idx2),'Color', [0.7 0.7 1.0], 'Linewidth',1)
plot(t_idx2,vec_orig(t_idx2),'r', 'Linewidth',2)
plot(t_idx2,vec_rec(t_idx2),'Color', [0 0.4 0], 'Linewidth',2)

axis([t_idx2(1) t_idx2(end) 0 3])

h = legend('data', 'exact', 'recovered');
set(h,'Interpreter','latex');
set(h, 'FontSize', 16);

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$x(t)$', 'Interpreter', 'latex', 'FontSize', 12)

mytitle = {['$T_2 = \lceil ' num2str(fem_reduce) ' T_1 \rceil = ' num2str(ceil(fem_reduce*length(vec_orig))) '$']};
h = title(mytitle);
set(h,'Interpreter','latex');
set(h, 'FontSize', 16);


hold off
