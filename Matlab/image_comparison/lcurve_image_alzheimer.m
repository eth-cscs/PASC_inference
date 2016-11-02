% proceed "shortinfo.txt" file and plot L-curve
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath('common')

K = 2;
annealing = 10;
width = 1683;height = 374;

M = csvread(['data/image_alzheimer/shortinfo_final.txt'],1,1);

myidx_unsorted = find(M(:,1)==width & M(:,2)==height & M(:,4) == K);
[~,sortidx] = sort(M(myidx_unsorted,3));
myidx = myidx_unsorted(sortidx);

epssqr = M(myidx,3);
lin_final = M(myidx,26);
qua_final = M(myidx,27);

%lin_final2 = log(lin_final)/log(10);
qua_final2 = log(qua_final./epssqr)/log(10);

lin_final2 = lin_final;
%qua_final2 = qua_final./epssqr;

lcurvefig = figure;
hold on

%xlabel('log(linear term)')
xlabel('linear term')
ylabel('log(quadratic term / epssqr)')
%ylabel('quadratic term / epssqr')

title(['L-curve , dim = ' num2str(width) 'x' num2str(height) ', K=' num2str(K) ' given, annealing=' num2str(annealing)])

plot(lin_final2, qua_final2, 'ro')
plot(lin_final2, qua_final2, 'r')

for i=1:length(epssqr)
   text(lin_final2(i),qua_final2(i),['eps = ' num2str(epssqr(i))]) 
end

hold off


