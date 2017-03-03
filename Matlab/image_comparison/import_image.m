% import results of denoised USI text image
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath(genpath(fullfile(pwd,'../common')))

width = 250;height = 150;
K = 2;

name_part = 'usi_test_250_150_epssqr10';
filename_orig = ['results/' name_part '_original.bin'];
filename_rec = ['results/' name_part '_recovered.bin'];
filename_gamma = ['results/' name_part '_gamma.bin'];

fig = figure;


% load vectors
vec_orig = loadbin(filename_orig);
vec_rec = loadbin(filename_rec);
vec_gamma = loadbin(filename_gamma);

% vector to matrices
Mat_orig = vec2mat(vec_orig,width);
Mat_rec = vec2mat(vec_rec,width);
for j = 1:K
   Mat_gamma{j} = vec2mat(vec_gamma(j:K:end),width);
end

% original image
ax(1) = subplot(1,K+2, 1);
hold on
Mat_orig(1:end,:) = Mat_orig(end:-1:1,:);
myimage_show(:,:,1) = Mat_orig;
myimage_show(:,:,2) = Mat_orig;
myimage_show(:,:,3) = Mat_orig;
image(min(max(myimage_show,0),1))
hold off

% recovered image
ax(2) = subplot(1,K+2, 2);
hold on
Mat_rec(1:end,:) = Mat_rec(end:-1:1,:);
myimage_show(:,:,1) = Mat_rec;
myimage_show(:,:,2) = Mat_rec;
myimage_show(:,:,3) = Mat_rec;
image(min(max(myimage_show,0),1))
hold off

% gamma
for j = 1:K
    ax(j+2) = subplot(1,K+2, 2+j);
    hold on
    Mat_gamma{j}(1:end,:) = Mat_gamma{j}(end:-1:1,:);
    myimage_show2(:,:,1) = Mat_gamma{j};
    myimage_show2(:,:,2) = Mat_gamma{j};
    myimage_show2(:,:,3) = Mat_gamma{j};
    image(min(max(myimage_show2,0),1))
    title(['gamma_' num2str(j)]);
    hold off
end

axis(ax,'image')
    

