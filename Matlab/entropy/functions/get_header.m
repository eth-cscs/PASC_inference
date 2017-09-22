function [ C ] = get_header( filename )

% Open the file
fid = fopen(filename);

% read first line
ln = fgetl(fid);

% split line
remain = ln;
C = [];

i = 1;
while ~isempty(remain)
   [token,remain] = strtok(remain, ',');
   
   % remove space at the begining
   C{i} = strtrim(token);
   i = i + 1;
end

end

