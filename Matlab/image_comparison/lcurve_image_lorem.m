% proceed "shortinfo.txt" file and plot L-curve
%
% Lukas Pospisil, USI, Lugano 2017
%

clear all

addpath(genpath(fullfile(pwd,'../common')))

K = 2;
width = 1000;height = 1000;

M = csvread(['data/image_lorem/shortinfo_final.txt'],1,1);

myidx_unsorted = find(M(:,1)==width & M(:,2)==height & M(:,3) == K);
[~,sortidx] = sort(M(myidx_unsorted,5));
myidx = myidx_unsorted(sortidx);

epssqr = M(myidx,5);
lin_final = M(myidx,25);
qua_final = M(myidx,26);

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
title(['L-curve , dim = ' num2str(width) 'x' num2str(height) ', K=' num2str(K) ' given'])

plot(lin_final2, qua_final2, 'ro')
plot(lin_final2, qua_final2, 'r')

for i=1:length(epssqr)
   text(lin_final2(i),qua_final2(i),['eps = ' num2str(epssqr(i))]) 
end


%set(gca,'Xscale','log')
set(gca,'Yscale','log')

%axis([300 500 10^0 10^9])


hold off


