function DG = get_graph_matrix_grid( x,y,T,alpha )

R = x*y;

disp('- preparing AG')
tic
AG = zeros(R,R);
for xi = 1:x
    for yi = 1:y
        if xi - 1 >= 1
            AG(get_r(xi,yi,x,y),get_r(xi-1,yi,x,y)) = 1;
            AG(get_r(xi-1,yi,x,y),get_r(xi,yi,x,y)) = 1;
        end
        if xi + 1 <= x
            AG(get_r(xi,yi,x,y),get_r(xi+1,yi,x,y)) = 1;
            AG(get_r(xi+1,yi,x,y),get_r(xi,yi,x,y)) = 1;
        end
        if yi - 1 >= 1
            AG(get_r(xi,yi,x,y),get_r(xi,yi-1,x,y)) = 1;
            AG(get_r(xi,yi-1,x,y),get_r(xi,yi,x,y)) = 1;
        end
        if yi + 1 <= y
            AG(get_r(xi,yi,x,y),get_r(xi,yi+1,x,y)) = 1;
            AG(get_r(xi,yi+1,x,y),get_r(xi,yi,x,y)) = 1;
        end
    end
end
AG = sparse(AG);
mytime = toc;
disp(['  finished in ' num2str(mytime) 's'])

disp('- preparing AT')
tic
AT = diag(ones(T-1,1),-1) + diag(ones(T-1,1),1);
AT = sparse(AT);
mytime = toc;
disp(['  finished in ' num2str(mytime) 's'])

disp('- preparing A')
tic
A = (1-alpha)*kron(eye(R),AT) + alpha*kron(AG,eye(T));
mytime = toc;
disp(['  finished in ' num2str(mytime) 's'])

disp('- preparing DG')
tic
DG = diag(sum(A,1))-2*A+diag(sum(A,2)); % sparse matrix operations
mytime = toc;
disp(['  finished in ' num2str(mytime) 's'])


end

function r = get_r(xi,yi,x,y)
    r = (yi-1)*x+xi;
end