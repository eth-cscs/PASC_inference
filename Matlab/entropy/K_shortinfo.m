% load data from provided shortinfo files
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath('functions')

Ks = [1 2 3];
k = 2;
T = 3000;

mycolors{1} = 'g';
mycolors{2} = 'b';
mycolors{3} = 'r';

figure
hold on
for i=1:length(Ks)
    K = Ks(i);
    
    filename = ['shortinfo/K' num2str(K) 'k2.txt'];
    C = get_header( filename );
    M = csvread(filename,1,0);

%    epssqrs = M(:,get_id_header(C,'epssqr'));
    epssqrs = 1./M(:,get_id_header(C,'epssqr'));
    LL = M(:,get_id_header(C,'SPGQP fx_linear'));
    nbins = M(:,get_id_header(C,'Gamma nbins'));
    
%    Ls_idx = 2*(K*k + T*(K-1)) + 2 * LL;
%    Ls_idx = 2*(K*k) + 2 * LL;
%    Ls_idx = 2*(K*k + K*nbins - 1) + 2 * T * LL;
%    Ls_idx = 2 * T * LL;
    Ls_idx = 2*(K*k + K*nbins - 1);

    
    ax = plot(epssqrs,Ls_idx,mycolors{i});
    
    
end

legend('K=1','K=2','K=3')

set(gca,'xscale','log');
hold off

