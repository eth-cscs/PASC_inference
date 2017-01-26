T = 100;

x_in = sin(20*[0:T-1]/T);

% reduction matrix
T2 = floor((T-1)/2.1);
%R = zeros(T2,T-1);
%for i=1:T2
%    idx = 2*i - 1;
%    R(i,idx:idx+2) = [1 2 1];
%end
%R = 1/4*R;

x2 = reduce(x_in,T2);
x_out = prolongate(x2,T);

x4 = reduce2(x_in,T2);


%x3 = (R*x_in')';

figure

hold on
plot(0:T-1,x_in,'r')
%plot(0:T-1,x_out,'b')
plot(0:T-1,prolongate(x2,T),'b')
%plot(0:T-1,prolongate(x3,T),'m')
%plot(0:T-1,prolongate(x4,T),'g')
plot(0:T-1,prolongate2(x4,T),'k')
plot(0:T-1,prolongate3(x4,T),'m')
hold off




