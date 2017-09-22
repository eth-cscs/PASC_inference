function mygrad  = get_grad( LM )

global Mom

mygrad = zeros(size(LM));

k = length(LM);

%compute normalization
mom_function = @(x,LM,order) power(x,order)*exp(-dot(LM,power(x*ones(size(1,k)),1:k)));
F = integral(@(x)mom_function(x,LM,0),-1,1, 'ArrayValued', true, 'RelTol',1e-12,'AbsTol',1e-10);

% theoretical moments
for i=1:k
  myI = integral(@(x)mom_function(x,LM,i),-1,1, 'ArrayValued', true, 'RelTol',1e-12,'AbsTol',1e-10);
  mygrad(i) = Mom(i) - myI/F;
end


end

