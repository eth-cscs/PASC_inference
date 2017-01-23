function x_out = prolongate( x_in, T )

T2 = size(x_in,2);
xdim = size(x_in,1);

x_out = zeros(xdim,T);

for t=0:T-1
    t2  = (t*(T2-1))/(T-1);
    
    %    center = (t*(T2-1))/(T-1);
    
    left_idx = max(floor(t2),0);
    right_idx = min(ceil(t2),T2-1);
    
    if right_idx-left_idx == 0
        x_out(t+1) = x_in(t2+1);
    else 
        a = (x_in(right_idx+1)-x_in(left_idx+1))/(right_idx-left_idx);
        b = x_in(left_idx+1) - a*left_idx;
    
        x_out(t+1) = a*t2+b;
    end
        
end

end

