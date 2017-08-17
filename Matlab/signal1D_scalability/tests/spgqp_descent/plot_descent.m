clear all

spgqpsolver_monitor
norm_gp1 = norm_gp;
fx1 = fx;
clear norm_gp fx

spgqpsolver_monitor2
norm_gp2 = norm_gp;
fx2 = fx;
clear norm_gp fx

nmb_it1 = length(fx1);
nmb_it2 = length(fx2);

figure
hold on
set(gca,'fontsize',12);

%brb = 100:2000;
brb = 1:2000;

plot(brb,norm_gp1(brb),'b')
%plot(brb,fx2(brb),'r')

set(gca,'Yscale','log')

xlabel('iterations', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\Vert g^P \Vert$', 'Interpreter', 'latex', 'FontSize', 12)

%h = legend('m = 20', 'm = 1');
%set(h,'Interpreter','latex');
%set(h, 'FontSize', 16);


hold off
