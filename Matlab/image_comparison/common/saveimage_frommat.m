function [image_vec] = saveimage_frommat( A, filename_image_out)

% show image
figure

% original image
ax(1) = subplot(1,1,1);
hold on
myimage_flip = A(end:-1:1,:);
myimage_show(:,:,1) = myimage_flip;
myimage_show(:,:,2) = myimage_flip;
myimage_show(:,:,3) = myimage_flip;
image(min(max(myimage_show,0),1))
hold off

axis(ax(1:1),'image')

rows = size(A,1);
cols = size(A,2);

% now save image as vector in petsc format
image_vec = zeros(1,rows*cols);
for y=1:rows
   image_vec((y-1)*cols+1:y*cols) = A(y,:);
end
savebin( filename_image_out, image_vec );



end

