function [ X_rec, theta, gamma, it, Lit ] = kmeans_H1_graph( X,K,DG,epssqr, theta_given )

if isempty(theta_given)
    maxit = 1000;
else
    maxit = 1;
end

T = size(X,3);
R = size(X,1)*size(X,2); % the dimension of image
n = 1;

disp('- reshaping X')
tic
x = zeros(R*T,1);
for xi=1:size(X,1)
    for yi=1:size(X,2)
        r = (yi-1)*size(X,1)+xi;
        x((r-1)*T + 1: r*T) = reshape(X(xi,yi,:),T,1);
    end
end
mytime = toc;
disp(['  finished in ' num2str(mytime) 's'])


% generate initial random feasible gamma 
gamma = rand(K*R*T,1); % [gamma^(1);gamma^(2);...;gamma^(K)] 

% initial gamma has to be feasible (i.e sum to one for each t and positive)
gamma_sum = zeros(T*R,1);
for k=1:K
    gamma_sum = gamma_sum+gamma((k-1)*R*T+1:k*R*T);
end
for k=1:K
    gamma((k-1)*R*T+1:k*R*T) = gamma((k-1)*R*T+1:k*R*T)./gamma_sum;
end

disp('- preparing QP objects')
tic
H = sparse(kron(epssqr*eye(K),DG));
B = kron(ones(1,K),eye(R*T));
c = ones(R*T,1);
l = zeros(K*R*T,1);
u = ones(K*R*T,1);
g = zeros(size(gamma));
mytime=toc;
disp(['  finished in ' num2str(mytime) 's'])



% settings of algorithm (quadprog)
options.Algorithm = 'interior-point-convex';
options.Display = 'none'; 

% here will be stored solution of model parameters - one for each cluster
theta = zeros(K,n); % [theta^(1),theta^(2),...,theta(K)]

% initial object function value
L = Inf;

it = 0; % iteration counter
while it < maxit % practical stopping criteria is present after computing new L (see "break")

    % compute Theta
    if isempty(theta_given)
        for k=1:K
           gammak = gamma((k-1)*T*R+1:k*T*R); % gamma^(k)
      
            if sum(gammak) ~= 0 % maybe gamma_k = 0 ? (i.e. this cluster is empty)
                theta(k) = dot(gammak,x)/sum(gammak);
            else
                theta(k) = 0;
            end
        end
    else
        theta = theta_given;
    end;

    % compute new gamma
    % prepare new linear term in QP, 
    % i.e. compute new residuum based on new theta
    for t = 1:R*T
        for k = 1:K
            g((k-1)*R*T+t) = dot(x(t) - theta(k),x(t) - theta(k)); % local model error
        end
    end

    % solve QP problem
    gamma = quadprog(H,g,[],[],B,c,l,u,gamma,options);
    
    % compute new function value
    Lold = L; % store old function value
    L = dot(g,gamma) + dot(H*gamma,gamma);
    deltaL = Lold - L; % the progress of objective function, Lold > L (?)
    
    % display some informations, so we know that something is happening
    disp([num2str(it) '. it: L = ' num2str(L) ', deltaL = ' num2str(deltaL)])
    
    % stopping (breaking) criteria 
    % based on sufficient decrease of objective funtion
    if abs(deltaL) < 10^-4
        break; % stop outer "while" cycle
    end
    
    it = it + 1; % increase number of iterations
    
    Lit(it) = L; % for postprocessing
end

% compute reconstructed data
x_rec = zeros(size(x));
for k=1:K
    gammak = gamma((k-1)*T*R+1:k*T*R);
    for i = 1:n
        x_rec = x_rec + theta(k)*gammak;
    end
end

% reshape to original matrix
disp('- reshape data back to matrix form')
tic
X_rec = zeros(size(X));
for xi=1:size(X,1)
    for yi=1:size(X,2)
        r = (yi-1)*size(X,1)+xi;
        X_rec(xi,yi,:) = reshape(x_rec((r-1)*T + 1: r*T),1,1,T);
    end
end
mytime=toc;
disp(['  finished in ' num2str(mytime) 's'])


end

