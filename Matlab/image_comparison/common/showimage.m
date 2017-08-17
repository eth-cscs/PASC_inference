function showimage(img, mytitle, subx, suby, subnmb, scaleimg)

% scale image to 0=min, 1=max
if scaleimg
    mymax = max(max(img));
    mymin = min(min(img));
    a = 1/(mymax - mymin);
    b = -a*mymin;
    img = a*img + b;
end

ax = subplot(subx,suby,subnmb);
hold on
imgflip = img(end:-1:1, :);
img_show(:,:,1) = imgflip;
img_show(:,:,2) = imgflip;
img_show(:,:,3) = imgflip;
image(min(max(img_show,0),1))
title(mytitle);
axis off
axis(ax,'image')

hold off


end

