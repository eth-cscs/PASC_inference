function [ y, gamma, qp_lin, qp_quad, L ] = compute_femh1( x, gamma, Theta, epssqr, fem_reduce )
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

T2 = ceil(T/fem_reduce);
[H1 B2] = get_H1(T2, length(Theta));
H2 = (epssqr/T2)*H1;

gamma2 = zeros(K*T2,1);
for k = 1:K
    gamma2((k-1)*T2+1:k*T2) = reduce2( gamma((k-1)*T+1:k*T), T2 );
end
% create H (constant in outer it)
%for k=1:K
%    gamma_k = gamma((k-1)*T+1:k*T);
%    H((k-1)*T+1:k*T,(k-1)*T+1:k*T) = (epssqr/T)*H((k-1)*T+1:k*T,(k-1)*T+1:k*T);
%end

% create equality constraints (constant in outer it)
c2 = ones(T2,1);

% lower bound
l2 = zeros(K*T2,1);

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
    
    b2 = zeros(K*T2,1);
    for k = 1:K
        b2((k-1)*T2+1:k*T2) = reduce2( b((k-1)*T+1:k*T), T2 );
    end
   
    % run QP Solver
    if true
        % use matlab solver
        [gamma2, Lnew, exitflag, output] = quadprog(H2,b2,[],[],B2,c2,l2,[],gamma,options);
        it_qp = output.iterations;
    else
        % use my SPGQP
        [gamma2, Lnew, it_qp] = sqpqp(2*H,-b,gamma, K, 1e-9);
    end
    
    gamma = zeros(K*T,1);
    for k = 1:K
        gamma((k-1)*T+1:k*T) = prolongate3( gamma2((k-1)*T2+1:k*T2), T );
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

