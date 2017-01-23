% proceed "shortinfo.txt" file and plot L-curve
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath(genpath(fullfile(pwd,'../common')))

K = 2;

color1 = [1.0 0.0 0.0];
color2 = [0.0 0.6 0.0];
color3 = [0.0 0.0 1.0];

M = csvread(['gpu_lcurve/shortinfo_final_1e-4.txt'],1,1);
myidx_unsorted = find(M(:,1)==K);
[~,sortidx] = sort(M(myidx_unsorted,2));
myidx = myidx_unsorted(sortidx);
epssqr1 = M(myidx,2);
abserr1 = M(myidx,3);
spgit1 = M(:,12);
lin_final1 = M(myidx,21+K);
qua_final1 = M(myidx,22+K)./epssqr1;


M = csvread(['gpu_lcurve/shortinfo_final_1e-5.txt'],1,1);
myidx_unsorted = find(M(:,1)==K);
[~,sortidx] = sort(M(myidx_unsorted,2));
myidx = myidx_unsorted(sortidx);
epssqr2 = M(myidx,2);
abserr2 = M(myidx,3);
spgit2 = M(:,12);
lin_final2 = M(myidx,21+K);
qua_final2 = M(myidx,22+K)./epssqr2;


M = csvread(['gpu_lcurve/shortinfo_final_1e-6.txt'],1,1);
myidx_unsorted = find(M(:,1)==K);
[~,sortidx] = sort(M(myidx_unsorted,2));
myidx = myidx_unsorted(sortidx);
epssqr3 = M(myidx,2);
abserr3 = M(myidx,3);
spgit3 = M(:,12);
lin_final3 = M(myidx,21+K);
qua_final3 = M(myidx,22+K)./epssqr3;


lcurvefig = figure;
set(gca,'fontsize',12);

hold on

%title(['L-curve , K=' num2str(K) ' given'])

plot(lin_final1, qua_final1, 'Color', color1, 'Marker', 'o')
pl1 = plot(lin_final1, qua_final1, 'r');

plot(lin_final2, qua_final2, 'Color', color2, 'Marker', 'o')
pl2 = plot(lin_final2, qua_final2, 'Color', color2);

plot(lin_final3, qua_final3, 'Color', color3, 'Marker', 'o')
pl3 = plot(lin_final3, qua_final3, 'Color', color3);


for i=1:9:length(epssqr3)-1
   mytext = text(lin_final3(i)-1.5,qua_final3(i),['$\varepsilon^2 = ' num2str(epssqr3(i)) '$']);
   set(mytext, 'Interpreter', 'latex');
   set(mytext, 'FontSize', 12);
   set(mytext, 'HorizontalAlignment', 'right');
   
end

%set(gca,'Xscale','log')
set(gca,'Yscale','log')

h = legend([pl1,pl2,pl3],'$\Vert g^P \Vert < 10^{-4}$','$\Vert g^P \Vert < 10^{-5}$','$\Vert g^P \Vert < 10^{-6}$');
set(h,'Interpreter','latex')
set(h, 'FontSize', 16);

xlabel('linear term', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('quadratic term $/ \varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)

hold off




abserrfig = figure;
set(gca,'fontsize',12);
hold on

xlabel('$\varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('absolute error', 'Interpreter', 'latex', 'FontSize', 12)

%title(['abserr , K=' num2str(K) ' given'])

plot(epssqr1, abserr1, 'Color', color1, 'Marker', 'o')
pl1 = plot(epssqr1, abserr1, 'Color', color1);

plot(epssqr2, abserr2, 'Color', color2, 'Marker', 'o')
pl2 = plot(epssqr2, abserr2, 'Color', color2);

plot(epssqr3, abserr3, 'Color', color3, 'Marker', 'o')
pl3 = plot(epssqr3, abserr3, 'Color', color3);

set(gca,'Xscale','log')

h = legend([pl1,pl2,pl3],'$\Vert g^P \Vert < 10^{-4}$','$\Vert g^P \Vert < 10^{-5}$','$\Vert g^P \Vert < 10^{-6}$');
set(h,'Interpreter','latex')
set(h, 'FontSize', 16);

hold off




abserrfig = figure;
set(gca,'fontsize',12);

hold on

xlabel('$\varepsilon^2$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('number of iterations', 'Interpreter', 'latex', 'FontSize', 12)

%title(['SPG-QP iterations'])

plot(epssqr1, spgit1, 'Color', color1, 'Marker', 'o')
pl1 = plot(epssqr1, spgit1, 'Color', color1);

plot(epssqr2, spgit2, 'Color', color2, 'Marker', 'o')
pl2 = plot(epssqr2, spgit2, 'Color', color2);

plot(epssqr3, spgit3, 'Color', color3, 'Marker', 'o')
pl3 = plot(epssqr3, spgit3, 'Color', color3);

set(gca,'Xscale','log')

h = legend([pl1,pl2,pl3],'$\Vert g^P \Vert < 10^{-4}$','$\Vert g^P \Vert < 10^{-5}$','$\Vert g^P \Vert < 10^{-6}$');
set(h,'Interpreter','latex')
set(h, 'FontSize', 16);

set(gca,'Yscale','log')

hold off
