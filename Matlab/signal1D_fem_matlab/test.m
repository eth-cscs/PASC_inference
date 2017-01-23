T = 1000;
T2 = 100;

x_in = sin(20*[0:T-1]/T);
x2 = reduce(x_in,T2);
x_out = prolongate(x2,T);

figure

hold on
plot(0:T-1,x_in,'r')
plot(0:T-1,x_out,'b')
%plot(T*[0:T2-1]/(T2-1),x2,'m')
hold off




