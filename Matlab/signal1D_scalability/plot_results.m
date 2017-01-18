% which solution to load
epssqr_load = '2e-4';

% include functions
addpath(genpath(fullfile(pwd,'../common')));

% filenames
filename_orig = 'data/sample_original.bin';
filename_rec = ['data/sample_' epssqr_load '_recovered.bin'];

% load vectors
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
t_idx = 5000:8000;

figure
hold on

plot(t_idx,vec_rec(t_idx),'m')
plot(t_idx,vec_orig(t_idx),'r')

axis([t_idx(1) t_idx(end) 0 3])

ax = legend('recovered','original');
leg = findobj(ax,'type','text');
set(leg,'FontSize',8)

xlabel('t')
ylabel('x(t)')

hold off
