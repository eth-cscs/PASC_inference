% Fourier analysis of the sun-planet distance for all the planets
figure
L=fix(0.6.*length(t));
NFFT=2^nextpow2(L);
Fs=1./abs((t(2)-t(1))./365.25./86400);
f = Fs/2*linspace(0,1,NFFT/2+1); 

for i=2:length(bodies)
    subplot(3,3,i-1);
    
    Y=fft(sqrt(sum((y(:,[1:3]+6.*(i-1))-y(:,[1:3]+6.*(1-1))).^2,2))  ,NFFT)/L;

    plot(1./f,(2*abs(Y(1:NFFT/2+1))),'r')
    set(gca,'yscale','log','xscale','log')
    xlabel('Period (years)');ylabel('Power');
    text(0.1,0.9,bodies{i},'fontsize',15,'units','normalized');
end

