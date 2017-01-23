% plot the results - original signal, solution and recoved signal

% include functions
addpath(genpath(fullfile(pwd,'../common')));

% filenames
filename_data = 'tests/sample_1e5_data.bin';
filename_orig = 'tests/sample_1e5_solution.bin';
filename_rec =  'tests/sample_1e5_recovered.bin';

% load vectors
vec_data = loadbin(filename_data);
vec_orig = loadbin(filename_orig);
vec_rec = loadbin(filename_rec);

if false
   for t=1:length(vec_rec)
      if vec_rec(t) < 1.5
          vec_rec(t) = 1;
      else
          vec_rec(t) = 2;
      end
   end
end

% what to plot
t_idx = 1:500000;

figure
hold on
set(gca,'fontsize',12);

plot(t_idx,vec_data(t_idx),'b')
plot(t_idx,vec_orig(t_idx),'r')
plot(t_idx,vec_rec(t_idx),'m')

%axis([t_idx(1) t_idx(end) 0 3])

h = legend('with noise', 'exact', 'recovered');
set(h,'Interpreter','latex');
set(h, 'FontSize', 16);

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$x(t)$', 'Interpreter', 'latex', 'FontSize', 12)

sizey = 10;
plot([60000,90000,90000,60000,60000],[1-sizey,1-sizey,2+sizey,2+sizey,1-sizey],'k')

hold off



t_idx2 = 60000:90000;

figure
hold on
set(gca,'fontsize',12);

plot(t_idx2,vec_orig(t_idx2),'r')
plot(t_idx2,vec_rec(t_idx2),'m')

axis([t_idx2(1) t_idx2(end) 0 3])

h = legend('exact', 'recovered');
set(h,'Interpreter','latex');
set(h, 'FontSize', 16);

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$x(t)$', 'Interpreter', 'latex', 'FontSize', 12)

hold off
