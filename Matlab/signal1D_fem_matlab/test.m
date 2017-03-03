T = 1000;
T2 = 17;

x_in = sin(20*[0:T-1]/T);

x_reduced0 = reduce0(x_in,T2);
x_out0 = prolongate0(x_reduced0,T);

x_reduced1 = reduce1(x_in,T2);
x_out1 = prolongate1(x_reduced1,T);

figure

hold on
plot(0:T-1,x_in,'r')
plot(0:T-1,x_out0,'b')
plot(0:T-1,x_out1,'Color',[0,0.6,0])
hold off




