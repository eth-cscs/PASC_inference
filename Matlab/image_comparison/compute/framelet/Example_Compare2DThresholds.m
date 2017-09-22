s1 = double(imread('st.tif'));
s = s1(:,:,3);
x = s + 20*randn(size(s));
t = 0:5:45;
e = den2(s,x,t);
de = double_den2(s,x,t);
re = real_den2(s,x,t);
ce = complex_den2(s,x,t);
clf
plot(t,e,'m-')
hold on
plot(t,de,'g-.')
hold on
plot(t,re,'b--')
hold on
plot(t,ce,'r:')
hold off
title('RMS Error vs. Threshold Point')
xlabel('Threshold Point');
ylabel('RMS Error');
legend('Standard 2-D', 'Double-Density 2-D', 'Double-Density Dual-Tree Real 2-D', 'Double-Density Dual-Tree Complex 2-D',0);