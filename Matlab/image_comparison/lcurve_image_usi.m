% proceed "shortinfo.txt" file and plot L-curve
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath('common')

M = csvread('data/image_usi/shortinfo_final2.txt',1,1);

K = 6;
annealing = 10;
%width = 250;height = 150;
%width = 500;height = 300;
width = 1000;height = 600;

noise_find = 10;
noise_print = 1.0;

myidx_unsorted = find(M(:,1)==width & M(:,2)==height & M(:,3)==noise_find & M(:,5) == K);
[~,sortidx] = sort(M(myidx_unsorted,4));
myidx = myidx_unsorted(sortidx);

epssqr = M(myidx,4);
lin_final = M(myidx,27);
qua_final = M(myidx,28);


%lin_final2 = log(lin_final)/log(10);
%qua_final2 = log(qua_final./epssqr)/log(10);

lin_final2 = lin_final;
qua_final2 = qua_final./epssqr;

lcurvefig = figure;
hold on

%xlabel('log(linear term)')
xlabel('linear term')
ylabel('quadratic term / epssqr')

title(['L-curve, dim = ' num2str(width) 'x' num2str(height) ', noise = ' num2str(noise_print) ', K=' num2str(K) ' given, annealing=' num2str(annealing)])

plot(lin_final2, qua_final2, 'ro')
plot(lin_final2, qua_final2, 'r')

for i=1:length(epssqr)
   text(lin_final2(i),qua_final2(i),['eps = ' num2str(epssqr(i))]) 
end

hold off


