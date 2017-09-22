clear all

width = 4;
height = 3;
T = 5;
xdim=3;
K=2;


addpath(genpath(fullfile(pwd,'../common')))

disp('- preparing movie')

disp(['- width  = ' num2str(width)])
disp(['- height = ' num2str(height)])
disp(['- T      = ' num2str(T)])
disp(['- xdim   = ' num2str(xdim)])
disp(['- K      = ' num2str(K)])

X_exact = zeros(height,width,T,xdim);
% 0.txyn
for t=1:T
    for x=1:width
        for y=1:height
            for n=1:xdim
                X_exact(y,x,t,n) = 0.1*t + 0.01*x + 0.001*y + 0.0001*n;
            end
        end
    end
end

% to vector TxdimR
X_vec = zeros(1,T*xdim*width*height);
for t=1:T
   for n=1:xdim
       for y=1:height
           for x=1:width
               X_vec((t-1)*width*height*xdim + (n-1)*width*height + (y-1)*width + x +1) = X_exact(y,x,t,n);
           end
       end
   end
end

savebin( ['movieseq_' num2str(width) 'w_' num2str(height) 'h_' num2str(T) 'T_' num2str(xdim) 'xdim.bin'], X_vec );

