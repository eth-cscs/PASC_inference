% proceed "shortinfo.txt" file and plot L-curve
%
% Lukas Pospisil, USI, Lugano 2017
%

% width,height,T,K,depth,
% epssqr,abserr,Theta0,Theta1,it all, 
% t all, t gamma update, t gamma solve, t theta update, t theta solve, 
% SPGQP it, SPGQP hessmult, SPGQP t all, SPGQP t project, SPGQP t matmult, 
% SPGQP t dot, SPGQP t update, SPGQP t stepsize, SPGQP t fs, SPGQP t other,
% SPGQP fx, SPGQP fx_linear, SPGQP fx_quadratic, SPGQP_sum it, SPGQP_sum hessmult, 
% SPGQP_sum t all, SPGQP_sum t project, SPGQP_sum t matmult, SPGQP_sum t dot, SPGQP_sum t update, 
% SPGQP_sum t stepsize, SPGQP_sum t fs, SPGQP_sum t other,

clear all

addpath(genpath(fullfile(pwd,'../common')))

problemsize = 256*128*1000;

K = 2;

M = csvread('results/shortinfo.txt',1,0);

myidx_unsorted = find(M(:,4) == K);
[~,sortidx] = sort(M(myidx_unsorted,6));
myidx = myidx_unsorted(sortidx);

epssqr = M(myidx,6);
abserr = M(myidx,7)/problemsize;

errorcurvefig = figure;
hold on

%xlabel('log(linear term)')
xlabel('epssqr')
ylabel('absolute error')
%ylabel('quadratic term / epssqr')
title(['error curve'])

plot(epssqr, abserr, 'ro')
plot(epssqr, abserr, 'r')

set(gca,'Xscale','log')
set(gca,'Yscale','log')

%axis([300 500 10^0 10^9])


hold off


