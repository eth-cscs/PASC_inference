% plot the results - original signal, solution and recoved signal

% include functions
addpath(genpath(fullfile(pwd,'../common')));
addpath(genpath(fullfile(pwd,'../../common')));

% filenames
filename_data = 'data/strong_10e7_data.bin';
filename_orig = 'data/strong_10e7_solution.bin';
filename_rec = 'data/strong_10e7_recovered.bin';

% load vectors
vec_data = loadbin(filename_data);
vec_orig = loadbin(filename_orig);
vec_rec = loadbin(filename_rec);

% what to plot
step = 20000;
t_idx1 = 1:10*step;%length(vec_data);
t_idx2 = 2*step:3*step;%length(vec_data);

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

if true
figure
hold on
set(gca,'fontsize',12);

plot(t_idx1,vec_data(t_idx1),'Color', [0.5 0.5 0.8], 'Linewidth',1)
plot(t_idx1,vec_orig(t_idx1),'r', 'Linewidth',1)
%plot(t_idx1,vec_rec(t_idx1), 'Color', [0 0.4 0], 'Linewidth',2)
plot([-1],[-1], 'Color', [0 0.4 0], 'Linewidth',1)

axis([t_idx1(1) t_idx1(end) 1.5-15 1.5+15])

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$x(t)$', 'Interpreter', 'latex', 'FontSize', 12)

x_min = min(t_idx2);
x_max = max(t_idx2);
y_min = 0;
y_max = 3;

%plot([60000,90000,90000,60000,60000],[1-sizey,1-sizey,2+sizey,2+sizey,1-sizey],'k')
plot([x_min,x_max,x_max,x_min,x_min],[y_min,y_min,y_max,y_max,y_min],'Color',[0,0,0],'LineWidth',2.0)

h = legend('data', 'exact','recovered');
set(h,'Interpreter','latex');
set(h, 'FontSize', 16);

hold off
end

if true
figure
hold on
set(gca,'fontsize',12);

%plot(t_idx2,vec_data(t_idx2),'Color', [0.7 0.7 1.0], 'Linewidth',1)
plot(t_idx2,vec_orig(t_idx2),'r', 'Linewidth',2)
plot(t_idx2,vec_rec(t_idx2),'Color', [0 0.4 0], 'Linewidth',2)

axis([t_idx2(1) t_idx2(end) 0 3])

%h = legend('exact', 'recovered');
%set(h,'Interpreter','latex');
%set(h, 'FontSize', 16);

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$x(t)$', 'Interpreter', 'latex', 'FontSize', 12)

hold off

end

