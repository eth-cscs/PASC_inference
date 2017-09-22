function [ value idx1 idx2 idx3] = minminmin( M )

[temp_value1, temp_idx1] = min(M);
[temp_value2, temp_idx2] = min(temp_value1);
[value, idx3] = min(temp_value2);

idx2 = temp_idx2(idx3);
idx1 = temp_idx1(:,idx2,idx3);

end

