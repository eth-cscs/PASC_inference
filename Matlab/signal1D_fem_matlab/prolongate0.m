function x_out = prolongate0( x_in, T )

T2 = length(x_in);

diff = (T2-1)/(T-1);
%diff = (T2-1)/T;

%pos_switch = round(1:(T-1)/T2:T);
%x_out = zeros(T,1);
%for b = 1:T2
%    x_out(pos_switch(b)+1:pos_switch(b+1)) = x_in(b)* ones(pos_switch(b+1) - pos_switch(b),1);
%end

x_out = zeros(T,1);
for t = 0:T-1
%    x_out(pos_switch(b)+1:pos_switch(b+1)) = x_in(b)* ones(pos_switch(b+1) - pos_switch(b),1);

%    disp([num2str(t) ': ' num2str(round(t*diff)+1)])

    x_out(t+1) = x_in(floor(t*diff)+1);
end


end

