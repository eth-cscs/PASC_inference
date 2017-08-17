% Fourier analysis of the sun's x velocity
figure('name','fourier analysis of sun''s velocity');

% some hacking to prepare the Fourier transform
L=fix(0.6.*length(t1));
NFFT=2^nextpow2(L);
Fs=1./abs((t1(2)-t1(1))./365.25./86400);
f = Fs/2*linspace(0,1,NFFT/2+1); 



% do the fourier transform
Y=fft(y(:,4)  ,NFFT)/L;

% plot result
plot(1./f,(2*abs(Y(1:NFFT/2+1))),'r')

% do the fourier transform
Y1=fft(y1(:,4)  ,NFFT)/L;
hold on;

% plot result
plot(1./f,(2*abs(Y1(1:NFFT/2+1))),'g')

set(gca,'yscale','log','xscale','log')
xlabel('Period (years)');ylabel('Power');
title('Fourier analysis of the Sun''s motion');
