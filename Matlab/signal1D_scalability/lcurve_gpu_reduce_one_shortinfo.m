% proceed "shortinfo.txt" file and plot L-curve
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

%filename = 'gpu_lcurve/shortinfo_final_small_noise.txt';
filename = 'gpu_lcurve/shortinfo_final.txt';
    
addpath(genpath(fullfile(pwd,'../common')))

mycolor{1} = [1.0 0.0 0.0];
mycolor{2} = [0.0 0.6 0.0];
mycolor{3} = [0.0 0.0 1.0];
mycolor{4} = [0.6 0.5 0.0];
mycolor{5} = [0.0 0.6 0.5];
mycolor{6} = [0.6 0.0 0.6];
mycolor{7} = [0.1 0.1 0.1];

M = csvread(filename,1,0);
K = 2;
n = 10^7;

fem_reduce = sort(unique(M(:,1)),1,'descend');
for i=1:length(fem_reduce)
    myidx_unsorted = find(M(:,1) == fem_reduce(i));
    [~,sortidx] = sort(M(myidx_unsorted,2));
    myidx = myidx_unsorted(sortidx);
    epssqr{i} = M(myidx,3);
    abserr{i} = M(myidx,4);
    spgit{i} = M(myidx,13);
    spgtime{i} = M(myidx,15);
    lin_final{i} = fem_reduce(i)*M(myidx,22+K);
    qua_final{i} = (1/fem_reduce(i))*M(myidx,23+K)./epssqr{i};
end

lcurvefig = figure;
set(gca,'fontsize',12);

hold on

%title(['L-curve , K=' num2str(K) ' given'])

clear pl mylegend
for i=1:length(fem_reduce)
    plot(lin_final{i}, abs(qua_final{i}), 'Color', mycolor{i}, 'Marker', 'o')
    pl(i) = plot(lin_final{i}, abs(qua_final{i}), 'Color', mycolor{i});
    mylegend(i) = {['$T_2 = \lceil ' num2str(fem_reduce(i)) ' T_1 \rceil = ' num2str(ceil(fem_reduce(i)*n)) '$']};
end

h = legend(pl,mylegend);
set(h,'Interpreter','latex')
set(h, 'FontSize', 16);

%set(gca,'Xscale','log')
set(gca,'Yscale','log')

xlabel('linear term', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('quadratic term $/ \varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)

hold off




abserrfig = figure;
set(gca,'fontsize',12);
hold on

xlabel('$\varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('absolute error', 'Interpreter', 'latex', 'FontSize', 12)

%title(['abserr , K=' num2str(K) ' given'])
for i=1:length(fem_reduce)
    plot(epssqr{i}, abserr{i}, 'Color', mycolor{i}, 'Marker', 'o')
    pl(i) = plot(epssqr{i}, abserr{i}, 'Color', mycolor{i});
end

h = legend(pl,mylegend);
set(h,'Interpreter','latex')
set(h, 'FontSize', 16);

set(gca,'Xscale','log')
%set(gca,'Yscale','log')

%h = legend(pl,mylegend);
%set(h,'Interpreter','latex')
%set(h, 'FontSize', 16);

hold off




figure
set(gca,'fontsize',12);

hold on

xlabel('$\varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('number of iterations', 'Interpreter', 'latex', 'FontSize', 12)

title(['SPG-QP iterations'])

for i=1:length(fem_reduce)
    plot(epssqr{i}, spgit{i}, 'Color', mycolor{i}, 'Marker', 'o')
    pl(i) = plot(epssqr{i}, spgit{i}, 'Color', mycolor{i});
end

h = legend(pl,mylegend);
set(h,'Interpreter','latex')
set(h, 'FontSize', 16);

set(gca,'Xscale','log')

%h = legend(pl,mylegend);
%set(h,'Interpreter','latex')
%set(h, 'FontSize', 16);

set(gca,'Yscale','log')

hold off



figure
set(gca,'fontsize',12);

hold on

xlabel('$\varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('solution time [s]', 'Interpreter', 'latex', 'FontSize', 12)

title(['SPG-QP solution time'])

for i=1:length(fem_reduce)
    plot(epssqr{i}, spgtime{i}, 'Color', mycolor{i}, 'Marker', 'o')
    pl(i) = plot(epssqr{i}, spgtime{i}, 'Color', mycolor{i});
end

h = legend(pl,mylegend);
set(h,'Interpreter','latex')
set(h, 'FontSize', 16);

set(gca,'Xscale','log')

%h = legend(pl,mylegend);
%set(h,'Interpreter','latex')
%set(h, 'FontSize', 16);

set(gca,'Yscale','log')

hold off

