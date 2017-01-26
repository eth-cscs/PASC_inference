function x_out = prolongate2( x_in, T )

T2 = length(x_in);

pos_switch = round(1:(T-1)/T2:T);
x_out = zeros(T,1);
for b = 1:T2
    x_out(pos_switch(b)+1:pos_switch(b+1)) = x_in(b)* ones(pos_switch(b+1) - pos_switch(b),1);
end

end

