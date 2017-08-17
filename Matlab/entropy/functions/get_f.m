function f = get_f( LM )

global Mom
k = size(Mom,1);
   
mom_function = @(x,LM,order) power(x,order)*exp(-dot(LM,power(x*ones(size(1,k)),1:k)));
F = integral(@(x)mom_function(x,LM,0),-1,1, 'ArrayValued', true);
F2 = myintegral(LM, 0);

disp(['F = ' num2str(F) ', F2 = ' num2str(F2)])

f = log(F) + dot(Mom,LM);

end

function f = myintegral(lambda, order)

Km = size(lambda,1);

mysum = 0;
N = 1000;
x = -1 + 2*rand(N,1);
y = -1 + 2*rand(N,1);

for i=1:N
   for j=1:N
       if y(j) <= power(x(i),order)*exp(-dot(lambda,power(x(i)*ones(size(1,Km)),1:Km)))
           mysum = mysum + 1;
       end
   end
end

f = mysum/(N*N);

end