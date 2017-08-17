function H = get_grad( LM )

global Mom

k = length(LM);
H = zeros(k,k);
    
I = zeros(2*k,1);
    
% compute normalization
mom_function = @(x,LM,order) power(x,order)*exp(-dot(LM,power(x*ones(size(1,k)),1:k)));
F = integral(@(x)mom_function(x,LM,0),-1,1, 'ArrayValued', true, 'RelTol',0,'AbsTol',1e-12);

% theoretical moments
for i = 1:2*k
    I(i) = integral(@(x)mom_function(x,LM,i),-1,1, 'ArrayValued', true, 'RelTol',0,'AbsTol',1e-12);
end

for i = 1:k
    for j = 1:k
        H(i,j) = I(i+j)/F - I(i)*I(j)/(F*F);
    end
end
            
end

