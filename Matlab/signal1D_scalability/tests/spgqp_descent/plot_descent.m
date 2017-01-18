clear all

spgqpsolver_monitor
norm_gp1 = norm_gp;

%spgqpsolver_monitor2
%norm_gp2 = norm_gp;

nmb_it = length(fx);

figure
hold on

plot(1:length(norm_gp1),norm_gp1,'r')
%plot(1:length(norm_gp1),norm_gp2,'b')

set(gca,'Yscale','log')

xlabel('iterations')
ylabel('$\Vert g^P \Vert$')

hold off
