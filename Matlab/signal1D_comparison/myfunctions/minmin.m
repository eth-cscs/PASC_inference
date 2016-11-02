function [ value idx1 idx2 ] = minmin( M )

[temp_value, temp_idx] = min(M);
[value, idx2] = min(temp_value);

idx1 = temp_idx(idx2);

end

