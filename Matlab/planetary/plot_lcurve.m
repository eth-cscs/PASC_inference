% proceed "shortinfo.txt" file and plot L-curve
%
% Lukas Pospisil, USI, Lugano 2017
%

% width,height,T,K,depth,
% epssqr,abserr,Theta0,Theta1,it all, 
% t all, t gamma update, t gamma solve, t theta update, t theta solve, 
% SPGQP it, SPGQP hessmult, SPGQP t all, SPGQP t project, SPGQP t matmult, 
% SPGQP t dot, SPGQP t update, SPGQP t stepsize, SPGQP t fs, SPGQP t other,
% 25
% 
% SPGQP fx, SPGQP fx_linear, SPGQP fx_quadratic, SPGQP_sum it, SPGQP_sum hessmult, 
% SPGQP_sum t all, SPGQP_sum t project, SPGQP_sum t matmult, SPGQP_sum t dot, SPGQP_sum t update, 
% SPGQP_sum t stepsize, SPGQP_sum t fs, SPGQP_sum t other,

clear all

addpath(genpath(fullfile(pwd,'../common')))

K = 2;

M = csvread('results/shortinfo.txt',1,0);

myidx_unsorted = find(M(:,4) == K);
[~,sortidx] = sort(M(myidx_unsorted,6));
myidx = myidx_unsorted(sortidx);

epssqr = M(myidx,6);
lin_final = M(myidx,27);
qua_final = M(myidx,28);

%lin_final2 = log(lin_final)/log(10);
qua_final2 = qua_final./epssqr;

lin_final2 = lin_final;
%qua_final2 = qua_final./epssqr;

lcurvefig = figure;
hold on

%xlabel('log(linear term)')
xlabel('linear term')
ylabel('quadratic term / epssqr')
%ylabel('quadratic term / epssqr')
title(['L-curve'])

plot(lin_final2, qua_final2, 'ro')
plot(lin_final2, qua_final2, 'r')

for i=1:length(epssqr)
   text(lin_final2(i),qua_final2(i),['eps = ' num2str(epssqr(i))]) 
end


%set(gca,'Xscale','log')
set(gca,'Yscale','log')

%axis([300 500 10^0 10^9])


hold off


