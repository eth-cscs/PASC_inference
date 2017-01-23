function [ y, gamma, qp_lin, qp_quad, L ] = compute_femh1( x, H1, Beq, gamma, Theta, epssqr )
%COMPUTE_FEMH1
%
% x        data
% H1       penalisation matrix (= 3diag(...))
% Beq      matrix of equality constraints (necessary for matlab quadprog)
% gamma    initial approximation of gamma
% Theta    given model parameters on each cluster
% epssqr   penalty-regularisation parameter
%

% in this implementation I assume 1D data
T = length(x);
K = length(Theta);

% create H (constant in outer it)
H = H1;
for k=1:K
    gamma_k = gamma((k-1)*T+1:k*T);
%    alpha_k = (1/T)*epssqr*(Theta(k) - dot(gamma_k,x)/sum(gamma))^2;
%    alpha_k = epssqr*Theta(k)^2;
    H((k-1)*T+1:k*T,(k-1)*T+1:k*T) = (epssqr/T)*H((k-1)*T+1:k*T,(k-1)*T+1:k*T);
%    H((k-1)*T+1:k*T,(k-1)*T+1:k*T) = epssqr*H((k-1)*T+1:k*T,(k-1)*T+1:k*T);
end

% create equality constraints (constant in outer it)
B = Beq;
c = ones(T,1);

% lower bound
l = zeros(K*T,1);

% this will be linear term in QP (changing in every outer it)
b = zeros(size(gamma));

% settings of algorithm (quadprog)
%options = optimoptions('quadprog');
options.Algorithm = 'interior-point-convex';
%options.Algorithm = 'trust-region-reflective';
%options.Algorithm = 'active-set';
options.Display = 'none';
options.ConstraintTolerance = 1e-12;


% initial object function value
L = Inf;

it = 0; % iteration counter
while it < 1000 % practical stopping criteria after computing new L (see "break")

    % compute Theta
    % -- Theta is given, therefore we will not compute it
    
    % compute new gamma
    % prepare new linear term in QP
    for t = 1:T
        for k = 1:K
            b((k-1)*T+t) = (x(t) - Theta(k))^2; % local model error
        end
    end
    % scale b
    b = (1/T)*b;
   
    % run QP Solver
    if true
        % use matlab solver
        [gamma, Lnew, exitflag, output] = quadprog(H,b,[],[],B,c,l,[],gamma,options);
        it_qp = output.iterations;
    else
        % use my SPGQP
        [gamma, Lnew, it_qp] = sqpqp(2*H,-b,gamma, K, 1e-9);
    end
    
    % compute function value
    Lold = L;
    L = Lnew;

%    disp(['    it=' num2str(it) ', L=' num2str(L)]);
    disp(['    it=' num2str(it) ', L=' num2str(L) ', it_qp=' num2str(it_qp)]);
    
    if abs(L - Lold) < 1e-9
        break; % stop outer "while" cycle
    end
    
    it = it + 1;
    
end    

% compute reconstructed signal
y = zeros(size(x));
for k=1:K
    y = y + Theta(k)*gamma((k-1)*T+1:k*T)';
end

% compute output values of linear term and quadratic term
qp_lin = dot(gamma,b);
qp_quad = (L - qp_lin)/epssqr;

end

