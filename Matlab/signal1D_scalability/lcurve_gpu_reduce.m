% proceed "shortinfo.txt" file and plot L-curve
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath(genpath(fullfile(pwd,'../common')))

K = 2;

mycolor{1} = [1.0 0.0 0.0];
mycolor{2} = [0.0 0.6 0.0];
mycolor{3} = [0.0 0.0 1.0];
mycolor{4} = [0.6 0.5 0.0];
mycolor{5} = [0.0 0.6 0.5];
mycolor{6} = [0.6 0.0 0.6];
mycolor{7} = [0.1 0.1 0.1];

filename{1} = 'gpu_lcurve/shortinfo_final_1e-6.txt';
filename{2} = 'gpu_lcurve/shortinfo_final_1e-6_reduce1e-1.txt';
filename{3} = 'gpu_lcurve/shortinfo_final_1e-6_reduce5e-2.txt';
filename{4} = 'gpu_lcurve/shortinfo_final_1e-6_reduce1e-2.txt';
filename{5} = 'gpu_lcurve/shortinfo_final_1e-6_reduce5e-3.txt';
filename{6} = 'gpu_lcurve/shortinfo_final_1e-6_reduce1e-3.txt';
filename{7} = 'gpu_lcurve/shortinfo_final_1e-6_reduce5e-3.txt';
%filename{5} = 'gpu_lcurve/shortinfo_final_1e-6_reduce1e-4.txt';

mylabel{1} = {'$T_2 = T_1 = 10^7$'};
mylabel{2} = {'$T_2 = 0.1 \cdot T_1 = 10^6$'};
mylabel{3} = {'$T_2 = 0.05 \cdot T_1 = 5 \cdot 10^5$'};
mylabel{4} = {'$T_2 = 0.01 \cdot T_1 = 10^5$'};
mylabel{5} = {'$T_2 = 0.005 \cdot T_1 = 5 \cdot 10^4$'};
mylabel{6} = {'$T_2 = 0.001 \cdot T_1 = 10^4$'};
mylabel{7} = {'$T_2 = 0.0005 \cdot T_1 = 5 \cdot 10^3$'};

for i=1:length(filename)
    M = csvread(filename{i},1,1);
    myidx_unsorted = find(M(:,1)==K);
    [~,sortidx] = sort(M(myidx_unsorted,2));
    myidx = myidx_unsorted(sortidx);
    epssqr{i} = M(myidx,2);
    abserr{i} = M(myidx,3);
    spgit{i} = M(:,12);
    lin_final{i} = M(myidx,21+K);
    qua_final{i} = M(myidx,22+K)./epssqr{i};
end


lcurvefig = figure;
set(gca,'fontsize',12);

hold on

%title(['L-curve , K=' num2str(K) ' given'])

clear pl mylegend
for i=1:length(filename)

    plot(lin_final{i}, abs(qua_final{i}), 'Color', mycolor{i}, 'Marker', 'o')
    pl(i) = plot(lin_final{i}, abs(qua_final{i}), 'Color', mycolor{i});
    mylegend(i) = mylabel{i};

%for i=1:9:length(epssqr3)-1
%   mytext = text(lin_final3(i)-1.5,qua_final3(i),['$\varepsilon^2 = ' num2str(epssqr3(i)) '$']);
%   set(mytext, 'Interpreter', 'latex');
%   set(mytext, 'FontSize', 12);
%   set(mytext, 'HorizontalAlignment', 'right');
%end
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
for i=1:length(filename)
    plot(epssqr{i}, abserr{i}, 'Color', mycolor{i}, 'Marker', 'o')
    pl(i) = plot(epssqr{i}, abserr{i}, 'Color', mycolor{i});

end

h = legend(pl,mylegend);
set(h,'Interpreter','latex')
set(h, 'FontSize', 16);

set(gca,'Xscale','log')

%h = legend(pl,mylegend);
%set(h,'Interpreter','latex')
%set(h, 'FontSize', 16);

hold off




abserrfig = figure;
set(gca,'fontsize',12);

hold on

xlabel('$\varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('number of iterations', 'Interpreter', 'latex', 'FontSize', 12)

%title(['SPG-QP iterations'])

for i=1:length(filename)
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
