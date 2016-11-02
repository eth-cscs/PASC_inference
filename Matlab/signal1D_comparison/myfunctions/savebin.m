function savebin( filename, A )

indices = 'int32';
precision = 'float64';

% get size of matrix
[m,n] = size(A);

% open file
fd = fopen(filename,'w','ieee-be');

if issparse(A) || min(size(A)) > 1
    disp('write matrix (1211216)')

    % compute number of zeros in each row
    n_nz = full(sum(A' ~= 0));

    % compute number of all nonzeros
    nz = sum(n_nz);

    % write header of file
    fwrite(fd,[1211216,m,n,nz,n_nz],'int32');

    % prepare CRS format
    [i,j,s] = find(A');

    % write nonzero cols indexes
    fwrite(fd,i-1,'int32');

    % write values
    fwrite(fd,s,'float64');
else
    disp('write vector (1211214)')
    
    % write header of file
    fwrite(fd,[1211214,m*n],'int32');
    fwrite(fd,A,'float64');

end

% close file
fclose(fd);

end



