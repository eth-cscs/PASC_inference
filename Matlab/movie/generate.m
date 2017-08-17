clear all

width = 30;
height = 20;
T = 10;

addpath(genpath(fullfile(pwd,'../common')))

disp('- preparing movie')

disp(['- width  = ' num2str(width)])
disp(['- height = ' num2str(height)])
disp(['- T      = ' num2str(T)])

X_exact = ones(height,width,T);
tic
for t=1:size(X_exact,3)
    c = [size(X_exact,1),size(X_exact,2)].*[0.3+t*0.05,0.5+t*0.01];
%    c = [size(X_exact,1),size(X_exact,2)].*[0.3,0.5];
    
    r = size(X_exact,1)*0.3;
    for i=1:size(X_exact,1)
        for j=1:size(X_exact,2)
            if (i-c(1))^2+(j-c(2))^2 <= r^2
                X_exact(i,j,t) = 0;
            end
        end
    end
end
%X = X_exact + 5*(-0.5 + rand(size(X_exact))); % add noise
%X = min(max(X,0),1); % cut data
mytime = toc;
disp(['  finished in ' num2str(mytime) 's'])

% to vector
X_vec = zeros(1,width*height*T);

tic
for t=1:T
   for y=1:height
        X_vec((t-1)*width*height + (y-1)*width+1:(t-1)*width*height + y*width) = X_exact(y,:,t);
   end
end
mytime = toc;
disp(['  finished in ' num2str(mytime) 's'])

savebin( ['movie_' num2str(width) 'w_' num2str(height) 'h_' num2str(T) 'T.bin'], X_vec );

