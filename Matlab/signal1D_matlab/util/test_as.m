function x_out = test_as(A,b,B,c,lb,ub,x0)

disp('############ active-set ############');

% Matlab solver
options.Algorithm = 'active-set';
options.Display = 'iter-detailed';
options.ConstraintTolerance = 1e-12;
%options.Display = 'iter';

A_full = full(A); % :(
B_full = full(B);

[x, fval, exitflag, output, lambda] = quadprog(A_full,-b,[],[],B_full,c,lb,ub,x0,options);
disp(['nmb of iterations = ' num2str(output.iterations)]);

% test iterations one after another
options.MaxIterations = 1;
options.Display = 'none';
x2 = x0;
for i = 1:output.iterations
%    [x2, fval2, exitflag2, output2, lambda2] = quadprog(A_full,-b,[],[],B_full,c,lb,ub,x2,options);
end

disp(['diff (control)    = ' num2str(norm(x - x2))]);



% print KKT
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

