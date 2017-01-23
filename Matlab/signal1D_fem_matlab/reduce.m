function x_out = reduce( x_in, T2 )

T = size(x_in,2);
xdim = size(x_in,1);

x_out = zeros(xdim,T2);

for t2=0:T2-1
    left  = ((t2-1)*(T-1))/(T2-1);
    right = ((t2+1)*(T-1))/(T2-1);
    center = (t2*(T-1))/(T2-1);
    
%    disp([num2str(left) ', ' num2str(right)])

    left_idx = max(floor(left),0);
    right_idx = min(ceil(right),T-1);
    
    nmb = right_idx-left_idx;
    coeffs = zeros(1,nmb);
    
    for i=0:nmb-1
       idx = left_idx + i;
       if idx <= center
            coeffs(i+1) =  (idx - left)/(center - left);
       end
       if idx > center
            coeffs(i+1) =  (idx - right)/(center - right);
       end
   
    end

    coeffs = coeffs/sum(coeffs);
    
%    if t==1 || t==T2
%       coeffs = coeffs*2; 
%    end
    
    mysum = 0;
    for i=0:nmb-1
        mysum = mysum + coeffs(i+1)*x_in(left_idx+i+1);
    end
    
    x_out(t2+1) = mysum;
end

end

