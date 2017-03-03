function x_out = test_ip(A,b,B,c,lb,ub,x0)

disp('############ interior-point-convex ############');

% Matlab solver
options.Algorithm = 'interior-point-convex';
options.Display = 'iter-detailed';
options.ConstraintTolerance = 1e-12;
options.LargeScale = 'on';
%options.Display = 'iter';

%[x, fval, exitflag, output, lambda] = quadprog(A,-b,[],[],B,c,lb,ub,x0,options);
[x, fval, exitflag, output, lambda] = quadprog(A,-b,[],[],B,c,lb,[],x0,options);
disp(['nmb of iterations = ' num2str(output.iterations)]);

% test iterations one after another
if false
    options.MaxIterations = 1;
    options.Display = 'none';
    x2 = x0;
    for i = 1:output.iterations
        [x2, fval2, exitflag2, output2, lambda2] = quadprog(A,-b,[],[],B,c,lb,ub,x2,options);
    end
    disp(['diff (control)    = ' num2str(norm(x - x2))]);
end


% print KKT
ll = (A*x - b + B'*lambda.eqlin);
norm(ll - lambda.lower)
kkt(1) = norm(A*x - b + B'*lambda.eqlin - lambda.lower + lambda.upper);
kkt(2) = norm(B*x-c);
kkt(3) = norm(min(x-lb,0));
kkt(4) = norm(max(x-ub,0));
kkt(5) = norm(min(lambda.lower,0));
kkt(6) = norm(min(lambda.upper,0));
kkt(7) = dot(lambda.lower,x-lb);
kkt(8) = dot(lambda.upper,x-ub);

disp('optimality: ');
for i=1:8
    disp(['- kkt' num2str(i) ' = ' num2str(kkt(i))]);
end

x_out = x;

end

