figure('name','inner planets');
plot(t./86400./365.25,sqrt(sum((y(:,[1:3]+6.*1)-y(:,1:3)).^2,2)),'k');hold on;
plot(t./86400./365.25,sqrt(sum((y(:,[1:3]+6.*2)-y(:,1:3)).^2,2)),'r')
plot(t./86400./365.25,sqrt(sum((y(:,[1:3]+6.*3)-y(:,1:3)).^2,2)),'g')
plot(t./86400./365.25,sqrt(sum((y(:,[1:3]+6.*4)-y(:,1:3)).^2,2)),'b')
xlabel('time (earth years)');
ylabel('Sun-planet distance (m)')
legend(bodies(2:5),-1);
xlim([0 10])

figure('name','outer planets');
plot(t./86400./365.25,sqrt(sum((y(:,[1:3]+6.*5)-y(:,1:3)).^2,2)),'k');hold on;
plot(t./86400./365.25,sqrt(sum((y(:,[1:3]+6.*6)-y(:,1:3)).^2,2)),'r')
plot(t./86400./365.25,sqrt(sum((y(:,[1:3]+6.*7)-y(:,1:3)).^2,2)),'g')
plot(t./86400./365.25,sqrt(sum((y(:,[1:3]+6.*8)-y(:,1:3)).^2,2)),'b')
plot(t./86400./365.25,sqrt(sum((y(:,[1:3]+6.*9)-y(:,1:3)).^2,2)),'c')
xlabel('time (earth years)');
ylabel('Sun-planet distance (m)')
legend(bodies(6:10),-1);
xlim([0 1000])

