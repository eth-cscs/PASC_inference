function [A] = loadbin( filename )

indices = 'int32';
precision = 'float64';

% open file
fd = fopen(filename,'r','ieee-be');

header = double(fread(fd,1,indices));

% if header == 1211216 % Petsc Mat Object 
%     header = double(read(fd,3,indices));
%     m      = header(1);
%     n      = header(2);
%     nz     = header(3);
%     if (nz == -1)
%       if arecomplex
%         s     = read(fd,2*m*n,precision);
%         iReal = 1:2:n*m*2-1;
%         iImag = iReal +1 ;
%         A     = complex(reshape(s(iReal),n,m)',reshape(s(iImag),n,m)') ;
%       else
%         s   = read(fd,m*n,precision);
%         A   = reshape(s,n,m)';
%       end
%     else
%       nnz = double(read(fd,m,indices));  %nonzeros per row
%       sum_nz = sum(nnz);
%       if(sum_nz ~=nz)
%         str = sprintf('No-Nonzeros sum-rowlengths do not match %d %d',nz,sum_nz);
%         error(str);
%       end
%       j   = double(read(fd,nz,indices)) + 1;
%       if arecomplex
%         s   = read(fd,2*nz,precision);
%       else 
%         s   = read(fd,nz,precision);
%       end
%       i   = ones(nz,1);
%       cnt = 1;
%       for k=1:m
%         next = cnt+nnz(k)-1;
%         i(cnt:next,1) = (double(k))*ones(nnz(k),1);
%         cnt = next+1;
%       end
%       if arecomplex
%         A = sparse(i,j,complex(s(1:2:2*nz),s(2:2:2*nz)),m,n,nz);
%       else
%         A = sparse(i,j,s,m,n,nz);
%       end
%     end
%     if arecell
%       result{l} = A;
%     else 
%       varargout(l) = {A};
%     end
% end 
    
if header == 1211214 % Petsc Vec Object
  m = double(fread(fd,1,indices)); 

  A = fread(fd,m,precision);
end

% if  header == 1211218 % Petsc IS Object
%     m = double(read(fd,1,indices));
%     v = read(fd,m,'int') + 1; % Indexing in Matlab starts at 1, 0 in PETSc
%     if arecell
%       result{l} = v;
%     else 
%       varargout(l) = {v};
%     end
% end

%   if header == 1211219 % Petsc Bag Object
%     b = PetscBagRead(fd);
%     if arecell
%       result{l} = b;
%     else 
%       varargout(l) = {b};
%     end

%   if header == 1211221 % Petsc DM Object
%     m  = double(read(fd,7,indices));
%     me = double(read(fd,5,indices));
%     b = [' dm ' int2str(m(3)) ' by ' int2str(m(4)) ' by ' int2str(m(5))];
%     if arecell
%       result{l} = b;
%     else 
%       varargout(l) = {b};
%     end

% close file
fclose(fd);


end



