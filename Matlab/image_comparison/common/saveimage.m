function [image_vec, nodes_vec] = saveimage( filename_in, filename_image_out, filename_graph_out, sigma )

% load image
myimage_orig = imread(filename_in);

% to grayscale
myimage(:,:) = (double(myimage_orig(:,:,1))+double(myimage_orig(:,:,2))+double(myimage_orig(:,:,3)))/3;

% to [0,1]
myimage = myimage/255;

% add noise
%mynoise = normrnd(0,sigma,size(myimage));
mynoise = sigma*randn(size(myimage));
%mynoise = sigma*rand(size(myimage));
myimage2 = myimage + mynoise;

% cut out of 0,1
%myimage2 = min(max(myimage2,0),1);

% create graph
rows = size(myimage2,1);
cols = size(myimage2,2);
nodes = zeros(2,rows*cols);
for y = 1:rows
    for x = 1:cols
        nodes(1,(y-1)*cols+x) = x-1;
        nodes(2,(y-1)*cols+x) = y-1;
    end
end

% show images and graph
figure

% original image
ax(1) = subplot(1,2,1);
hold on
myimage_flip = myimage(end:-1:1,:);
myimage_show(:,:,1) = myimage_flip;
myimage_show(:,:,2) = myimage_flip;
myimage_show(:,:,3) = myimage_flip;
image(min(max(myimage_show,0),1))
hold off

% image with noise
ax(2) = subplot(1,2,2);
hold on
myimage2_flip = myimage2(end:-1:1,:);
myimage_show(:,:,1) = myimage2_flip;
myimage_show(:,:,2) = myimage2_flip;
myimage_show(:,:,3) = myimage2_flip;
image(min(max(myimage_show,0),1))
hold off

% graph
%ax(3) = subplot(1,3,3);
%hold on
%plot(nodes(1,:),nodes(2,:),'r.')
%hold off

axis(ax(1:2),'image')


% now save images as vectors in petsc format
nodes_vec = [nodes(1,:) nodes(2,:)];
savebin( filename_graph_out, nodes_vec );

image_vec = zeros(1,rows*cols);
for y=1:rows
   image_vec((y-1)*cols+1:y*cols) = myimage2(y,:);
end
savebin( filename_image_out, image_vec );



end

