% import results of denoised USI text image
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all

addpath('common')

width = 250;height = 150;

epssqr{1}  = 0.000007;epssqr_load{1}  = '0.000007';
%epssqr{2}  = 0.000030;epssqr_load{2}  = '0.000030';
%epssqr{3}  = 0.000500;epssqr_load{3}  = '0.000500';

K = 2;
noise = '04';
annealing = 10;
arch = 'GPU1_N1_Nthreads1_Ngpu1';

fig = figure;

for i = 1:length(epssqr)
   
    name_part = ['test_image_usi_w' num2str(width) '_h' num2str(height) '_noise' noise '_epssqr' num2str(epssqr_load{i}) '_K' num2str(K) '_arch' arch];

    filename_orig = ['data/image_usi/results/' name_part '_original.bin'];
    filename_rec = ['data/image_usi/results/' name_part '_recovered.bin'];
    filename_gamma = ['data/image_usi/results/' name_part '_gamma.bin'];

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
    ax(1) = subplot(length(epssqr),K+2, (i-1)*(K+2)+1);
    hold on
    text(0.5,0.5,['eps = ' num2str(epssqr{i})])
    axis off
    hold off
    
    % recovered image
    ax(2) = subplot(length(epssqr),K+2, (i-1)*(K+2)+2);
    hold on
    Mat_rec(1:end,:) = Mat_rec(end:-1:1,:);
    myimage_show(:,:,1) = Mat_rec;
    myimage_show(:,:,2) = Mat_rec;
    myimage_show(:,:,3) = Mat_rec;
    image(min(max(myimage_show,0),1))

    if i == 1
        title('recovered');
    end
    
    hold off

    % gamma
    for j = 1:K
        ax(j+2) = subplot(length(epssqr),K+2, (i-1)*(K+2)+2+j);
        hold on
        Mat_gamma{j}(1:end,:) = Mat_gamma{j}(end:-1:1,:);
        myimage_show2(:,:,1) = Mat_gamma{j};
        myimage_show2(:,:,2) = Mat_gamma{j};
        myimage_show2(:,:,3) = Mat_gamma{j};
        image(min(max(myimage_show2,0),1))
        if i == 1
            title(['gamma_' num2str(j)]);
        end
            
        hold off
    end

    axis(ax,'image')
    
end

