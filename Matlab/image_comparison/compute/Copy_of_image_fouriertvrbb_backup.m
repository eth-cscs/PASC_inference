function [ Fout, it ] = image_fouriertvrbb(F, S, lambda, epsilon, niter)

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
G = -div( Gr./ repmat( sqrt( epsilon^2 + d.^2 ), [1 1 2]) );

%d = sqrt( epsilon^2 + sum3(Gr.^2,3));
%G = -div( Gr./ repmat( d, [1 1 2]) );

e = Phi(fTV,h)-F;
dd = Phi(e,h) + lambda*G;

% we use E for postprocessing ( to see descent )
clear E
it = 1;
for it = 1:niter
    % step
    fTVold = fTV;
    fTV = fTV - alpha*dd;
    
    % gradient
    ddold = dd;
    Gr = grad(fTV);

    d = sqrt(sum3(Gr.^2,3)); %sqrt( epsilon^2 + sum(Gr.^2,3) );
    G = -div( Gr./ repmat( sqrt( epsilon^2 + d.^2 ), [1 1 2]) );

%    d = sqrt( epsilon^2 + sum3(Gr.^2,3));
%    G = -div( Gr./ repmat( d, [1 1 2]) );    
    
    e = Phi(fTV,h)-F;
    dd = Phi(e,h) + lambda*G;
    
    s = fTV - fTVold;
    alpha = dot(s,s)/dot(s,dd-ddold);
    
    % energy
    E(it) = 1/2*norm(e, 'fro')^2 + lambda*sum(d(:));
%    disp([num2str(it) ',' num2str(E(it)) ',' num2str(norm(s)) ',' num2str(norm(dd-ddold)) ',' num2str(norm(dd))])

    if it>1
%        disp([num2str(it) ': ' num2str(norm(E(it)-E(it-1)))])
        if norm(E(it)-E(it-1)) < sqrt(length(fTV))
            break
        end
    end
end

% set output
Fout = fTV;

disp(['number of iterations: ' num2str(it)])

if false
    figure
    hold on
    showimage(Fout, ['fourier TVR, s=' num2str(S) ', lambda=' num2str(lambda) ', epsilon=' num2str(epsilon)], 1, 1, 1, true)
    hold off
end

end
