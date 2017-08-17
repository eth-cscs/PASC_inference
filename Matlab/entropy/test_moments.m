mu = 2;
sigma = 3;
x = mu + sqrt(sigma)*randn(1000,1);

T = length(x);
for i = 0:5
    disp(['Mom_' num2str(i) ' = ' num2str(sum(x.^i)/T)])
end
