data = ...
[ ...
0.0001,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1,10;...
0.488743,0.474353,0.468427,0.464882,0.46259,0.461054,0.459839,0.458619,0.457734,0.456962,0.456406,0.453983,0.453834,0.454096,0.454501,0.454971,0.455661,0.456196,0.456638,0.45712,0.466243,0.467766;...
];

figure
hold on
plot(data(1,:),data(2,:),'r')
plot(data(1,:),data(2,:),'r.')
xlabel('epssqr')
ylabel('err')
set(gca,'Xscale','log')
hold off