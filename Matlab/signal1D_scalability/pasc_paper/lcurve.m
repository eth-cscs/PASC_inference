% proceed "shortinfo.txt" file and plot L-curve
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath(genpath(fullfile(pwd,'../../common')))

K = 2;
T = 10^7;

M = csvread(['results/shortinfo_big.txt'],1,1);
myidx_unsorted = find(M(:,1)==K);
[~,sortidx] = sort(M(myidx_unsorted,2));
myidx = myidx_unsorted(sortidx);
epssqr1 = M(myidx,2);
abserr1 = M(myidx,3);
energy1 = M(myidx,4);
spgit1 = M(:,13);
lin_final1 = M(myidx,22+K);
qua_final1 = M(myidx,23+K)./epssqr1;


lcurvefig = figure;
set(gca,'fontsize',12);

hold on

plot(lin_final1, qua_final1, 'Color', 'b', 'Marker', 'o', 'LineWidth', 2.0)
pl1 = plot(lin_final1, qua_final1, 'b', 'LineWidth', 2.0);

for i=1:9:length(epssqr1)-1
   mytext = text(lin_final1(i)-0.1,qua_final1(i),['$\varepsilon^2 = ' num2str(epssqr1(i)) '$']);
   set(mytext, 'Interpreter', 'latex');
   set(mytext, 'FontSize', 12);
   set(mytext, 'HorizontalAlignment', 'right');
   
end

axis([min(lin_final1)-0.5 max(lin_final1)+0.5 min(qua_final1) max(qua_final1)])

set(gca,'Xscale','log')
set(gca,'Yscale','log')

%h = legend([pl1,pl2,pl3],'$\Vert g^P \Vert < 10^{-4}$','$\Vert g^P \Vert < 10^{-5}$','$\Vert g^P \Vert < 10^{-6}$');
%set(h,'Interpreter','latex')
%set(h, 'FontSize', 16);

xlabel('linear term', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('quadratic term $/ \varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)


hold off




abserrfig = figure;
set(gca,'fontsize',12);
hold on

xlabel('$\varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('absolute error', 'Interpreter', 'latex', 'FontSize', 12)

%title(['abserr , K=' num2str(K) ' given'])

plot(epssqr1, abserr1/T, 'Color', 'r', 'Marker', 'o', 'LineWidth', 2.0)
pl1 = plot(epssqr1, abserr1/T, 'Color', 'r', 'LineWidth', 2.0);

set(gca,'Xscale','log')

%h = legend([pl1,pl2,pl3],'$\Vert g^P \Vert < 10^{-4}$','$\Vert g^P \Vert < 10^{-5}$','$\Vert g^P \Vert < 10^{-6}$');
%set(h,'Interpreter','latex')
%set(h, 'FontSize', 16);

hold off

[~,best_idx] = min(abserr1);
disp(['best epssqr = ' num2str(epssqr1(best_idx))]);
