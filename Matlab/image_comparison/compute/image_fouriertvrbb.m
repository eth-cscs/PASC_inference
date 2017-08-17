function [ Fout, it ] = image_fouriertvrbb(F, S, lambda, epsilon, niter)

time_all = tic;
time_fft = 0;
time_grad = 0;
time_div = 0;
time_fro = 0;

[height, width] = size(F);

% Kernel
x = [0:width/2-1, -width/2:-1];
y = [0:height/2-1, -height/2:-1];
[X,Y] = meshgrid(x,y);
h = exp( (-X.^2 - Y.^2)/(2*S^2) );
h = h/sum(h(:));

% Initialization.
fTV = F;%reshape(F,height*width,1);

% Compute the TV norm, usefull to keep track of its decay through iterations.
%tv = sum(d(:));

% define iteration mapping
Phi = @(x,h)real(ifft2(fft2(x).*fft2(h)));

tau = 1.9 / ( 1 + lambda * 8 / epsilon);
alpha = tau;

Gr = grad(fTV);

d = sqrt(sum3(Gr.^2,3)); %sqrt( epsilon^2 + sum(Gr.^2,3) );

time_div2 = tic;
 G = -div( Gr./ repmat( sqrt( epsilon^2 + d.^2 ), [1 1 2]) );
time_div = time_div + toc(time_div2);

%d = sqrt( epsilon^2 + sum3(Gr.^2,3));
%G = -div( Gr./ repmat( d, [1 1 2]) );

time_fft2 = tic;
myPhi = Phi(fTV,h);
time_fft = time_fft + toc(time_fft2);
e = myPhi - F;

time_fft2 = tic;
myPhi = Phi(e,h);
time_fft = time_fft + toc(time_fft2);
dd = myPhi + lambda*G;

% we use E for postprocessing ( to see descent )
clear E
it = 1;
for it = 1:niter
    % step
    fTVold = fTV;
    fTV = fTV - alpha*dd;
    
    % gradient
    ddold = dd;

    time_grad2 = tic;
    Gr = grad(fTV);
    time_grad = time_grad + toc(time_grad2);

    d = sqrt(sum3(Gr.^2,3)); %sqrt( epsilon^2 + sum(Gr.^2,3) );

    time_div2 = tic;
     G = -div( Gr./ repmat( sqrt( epsilon^2 + d.^2 ), [1 1 2]) );
    time_div = time_div + toc(time_div2);

%    d = sqrt( epsilon^2 + sum3(Gr.^2,3));
%    G = -div( Gr./ repmat( d, [1 1 2]) );    

    time_fft2 = tic;
    myPhi = Phi(fTV,h);
    time_fft = time_fft + toc(time_fft2);
    e = myPhi - F;

    time_fft2 = tic;
    myPhi = Phi(e,h);
    time_fft = time_fft + toc(time_fft2);
    dd = myPhi + lambda*G;
    
    s = fTV - fTVold;
%    alpha = dot(s,s)/dot(s,dd-ddold);
    if it > 10
        alpha = 0.01*sum2(s.*s)/sum2(s.*(dd-ddold));
    end
        
    % energy
    time_fro2 = tic;
    E(it) = 1/2*norm(e, 'fro')^2 + lambda*sum(d(:));
    time_fro = time_fro + toc(time_fro2);
    disp([num2str(it) ':' num2str(E(it)) ',' num2str(norm(s)) ',' num2str(norm(dd-ddold)) ',' num2str(norm(dd)) ', alpha=' num2str(alpha)])

    if it>1
        disp([num2str(it) ': norm(E-E_old)=' num2str(norm(E(it)-E(it-1)))])
        if norm(E(it)-E(it-1)) < sqrt(length(fTV))
            break
        end
    end
end

% set output
Fout = fTV;

time_all = toc(time_all);

disp(['number of iterations: ' num2str(it)])
disp([' time all  :' num2str(time_all)])
disp([' time ftt  :' num2str(time_fft)])
disp([' time div  :' num2str(time_div)])
disp([' time grad :' num2str(time_grad)])
disp([' time fro  :' num2str(time_fro)])

if false
    figure
    hold on
    showimage(Fout, ['fourier TVR, s=' num2str(S) ', lambda=' num2str(lambda) ', epsilon=' num2str(epsilon)], 1, 1, 1, true)
    hold off
end

end

function [out] = sum2(A)
    out = sum(sum(A));
end
