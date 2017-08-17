function [vec_orig, vec_rec, vec_gamma] = openimage( K, width, filename_orig, filename_rec, filename_gamma, image_title, gamma_title, image_filename, gamma_filename )

% load vectors
vec_orig = loadbin(filename_orig);
vec_rec = loadbin(filename_rec);
vec_gamma = loadbin(filename_gamma);

% vector to matrices
Mat_orig = vec2mat(vec_orig,width);
Mat_rec = vec2mat(vec_rec,width);
for i = 1:K
    Mat_gamma{i} = vec2mat(vec_gamma(i:K:end),width);
end

% show images
imagefig = figure;

% original image
ax(1) = subplot(1,2,1);
hold on
Mat_orig(1:end,:) = Mat_orig(end:-1:1,:);
myimage_show(:,:,1) = Mat_orig;
myimage_show(:,:,2) = Mat_orig;
myimage_show(:,:,3) = Mat_orig;
image(min(max(myimage_show,0),1))
title('original');
hold off

% recovered image
ax(2) = subplot(1,2,2);
hold on
Mat_rec(1:end,:) = Mat_rec(end:-1:1,:);
myimage_show(:,:,1) = Mat_rec;
myimage_show(:,:,2) = Mat_rec;
myimage_show(:,:,3) = Mat_rec;
image(min(max(myimage_show,0),1))
title('recovered');
hold off

axis(ax,'image')
mtit(imagefig,image_title);
set(imagefig, 'PaperUnits', 'points');
set(imagefig, 'Position', [100, 100, 2*250, 270]);
set(imagefig, 'PaperPosition', [100, 100, 2*250, 270]);
saveas(imagefig,image_filename)

% gammas
gammafig = figure;

% original image
for i = 1:K
    ax2(i) = subplot(1,K,i);
    hold on
    Mat_gamma{i}(1:end,:) = Mat_gamma{i}(end:-1:1,:);
    myimage_show2(:,:,1) = Mat_gamma{i};
    myimage_show2(:,:,2) = Mat_gamma{i};
    myimage_show2(:,:,3) = Mat_gamma{i};
    image(min(max(myimage_show2,0),1))
%    image(myimage_show2)
    title(['gamma' num2str(i)]);
    hold off
end

axis(ax2,'image')
mtit(gammafig,gamma_title);
set(gammafig, 'PaperUnits', 'points');
set(gammafig, 'Position', [100, 100, K*250, 270]);
set(gammafig, 'PaperPosition', [100, 100, K*250, 270]);
saveas(gammafig,gamma_filename)

end

