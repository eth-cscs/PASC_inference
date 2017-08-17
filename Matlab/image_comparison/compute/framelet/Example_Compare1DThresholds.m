s1 = double(imread('peppers.jpg'));	% load image as a double
s = s1(:,2,3);			% convert to a 2-D image
x = s + 20*randn(size(s));		% add Gaussian noise to the original image
t = 0:5:45;         			% threshold range (0-45), increment by 5
e = double_den1(s,x,t);    		% using double-density method
de = double_dual1(s,x,t);  		% using double-density dual-tree method
figure(1)
plot(t,e,'blue')			% plot double-density (in blue)
hold on
plot(t,de,'red')			% plot double-density dual-tree (in red)
title('RMS Error vs. Threshold Point')
xlabel('Threshold Point');
ylabel('RMS Error');
legend('Double-Density 1-D','Double-Density Dual-Tree 1-D',0);