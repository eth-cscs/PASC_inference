function x_out = reduce0( x_in, T2 )

T = length(x_in);

%pos_switch = round(1:(T-1)/T2:T);
diff = T/T2;

x_out = zeros(T2,1);
for t2 = 0:T2-1
    left = round(t2*diff);
%    right = min(round((t2+1)*diff), T2-1);
    right = round((t2+1)*diff);
    
    mysum = sum(x_in((left:right-1)+1)) / (right-left);
    
    x_out(t2+1) = mysum;%/(pos_switch(b+1)-pos_switch(b));
end


end

