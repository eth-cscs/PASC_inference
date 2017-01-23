% load data/generate new testing data and show how matlab implementation of
% FEMH1 regularisation works
%
% Lukas Pospisil, USI, Lugano 2016
%

clear all
close all

addpath(genpath(fullfile(pwd,'../common')))
addpath('compute_femh1') % add functions for femh1 computation 

Nreduce = 2;

% load data or create new
C_true=[2*ones(1,60) ones(1,40) 2*ones(1,150) ones(1,50) 2*ones(1,150)];

C_true = [C_true];
sigma=1; 
C_orig=C_true+sigma*randn(size(C_true)); % gaussian

% reduce original problem
C = reduce( C_orig, ceil(length(C_orig)/Nreduce) );


Theta = [1.0, 2.0]; % given Theta, K = length(Theta)
epssqrs = 10.^(0:0.2:7); % used penalties
T = length(C_true);

C_err = Inf*ones(size(epssqrs)); % here store errors (to plot final error curve)
L_value = Inf*ones(size(epssqrs)); % here store values of L (to final plot)
qp_lin = zeros(size(epssqrs)); % here store values of linear term of L (to plot final L-curve)
qp_quad = zeros(size(epssqrs)); % here store values of linear term of L (to plot final L-curve)

best_epssqr = Inf;
best_epssqr_id = 1;

gamma = get_random_gamma0( length(C), length(Theta) ); % random initial gamma for first epssqr 
[H1 B] = get_H1(length(C), length(Theta)); % I want to reuse matrices, it is not neccessary to assemble them for each epssqr
for i = 1:length(epssqrs)
	% get penalty
	epssqr = epssqrs(i);

    % show some info to prove that Matlab is not sleeping
    disp([' - epssqr = ' num2str(epssqr)])
    
	% run the algorithm (reuse previous gamma as initial approximation)
    [C_filtered, gamma, qp_lin(i), qp_quad(i), L ] = compute_femh1( C, H1, B, gamma, Theta, epssqr );
    
	% compute error
%    C_err(i) = mean(abs(C_filtered-C_true));
    C_filtered_orig = prolongate(C_filtered,length(C_orig));
    C_err(i) = norm(C_filtered_orig-C_true,1);
    L_values(i) = L;

	% if this choice of penalty is the best (subject to error), then store the solution
    if C_err(i) < best_epssqr
		best_epssqr_id = i;
        best_C_filtered = C_filtered_orig;
        best_epssqr = C_err(best_epssqr_id);
	end
end


% plot the error curve (epssqr vs. err)
figure
subplot(1,3,1);
hold on
title('Error curve');
xlabel('epssqr');
ylabel('Filtering Error E_t[||x_{true}(t)-x_{filtered}(t)||]')
% plot curve
plot(epssqrs,C_err,'b.-','LineWidth',1.0);
% plot the best choice
plot(epssqrs(best_epssqr_id),C_err(best_epssqr_id),'mo', 'MarkerSize', 15)
set(gca,'xscale','log','yscale','log')
hold off


% plot the L curve (linear_part vs. quadratic_part)
subplot(1,3,2);
hold on
title('L-curve');
xlabel('linear term (model error)');
ylabel('quadratic term (penalty)')
%set(gca,'xscale','log','yscale','log')
% plot curve
plot(qp_lin,qp_quad,'r.-','LineWidth',1.0);
% plot the best choice
plot(qp_lin(best_epssqr_id),qp_quad(best_epssqr_id),'mo', 'MarkerSize', 15)
% plot text
for i=1:length(epssqrs)
   text(qp_lin(i),qp_quad(i),['eps = ' num2str(epssqrs(i))]) 
end
set(gca,'xscale','log','yscale','log')
hold off


% plot values of L
subplot(1,3,3);
hold on
title('values of L');
xlabel('epssqr');
ylabel('L')
% plot curve
plot(epssqrs,L_values,'b.-','LineWidth',1.0);
% plot the best choice
plot(epssqrs(best_epssqr_id),L_values(best_epssqr_id),'mo', 'MarkerSize', 15)
set(gca,'xscale','log')
hold off


% plot solution
figure
hold on
title(['solution, epssqr = ' num2str(epssqrs(best_epssqr_id))]);
xlabel('t');
ylabel('x(t)')
% plot signals
plot(1:length(C_orig),C_orig,'b-','LineWidth',1.0);
plot(1:length(C_true),C_true,'r--','LineWidth',2.0);
plot(1:length(C_true),best_C_filtered,'-','LineWidth',2.0,'Color',[0.0,0.4,0.0]);
legend('with noise','original', 'reconstructed');

% the noise is large, I am interested more in signal, not in noise
axis([1 length(C_true) 0 3])
hold off



