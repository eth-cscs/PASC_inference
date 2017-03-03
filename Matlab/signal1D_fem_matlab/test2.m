T = 3000;

x_in = sin(20*[0:T-1]/T);

% reduction matrix
T2 = floor(T*0.07);
x1 = reduce2(x_in,T2);

figure

hold on
plot(0:T-1,x_in,'r')
plot(0:T-1,prolongate3(x1,T),'b')
hold off




