function savebin( filename, A )

indices = 'int32';
precision = 'float64';

if ~issparse(A)
   A = sparse(A);
end

% get size of matrix
[m,n] = size(A);

% compute number of zeros in each row
n_nz = full(sum(A' ~= 0));

% compute number of all nonzeros
nz = sum(n_nz);

% open file
fd = fopen(filename,'w','ieee-be');

% write header of file
fwrite(fd,[1211216,m,n,nz,n_nz],'int32');

% prepare CRS format
[i,j,s] = find(A');

% write nonzero cols indexes
fwrite(fd,i-1,'int32');

% write values
fwrite(fd,s,'float64');

% close file
fclose(fd);




end



