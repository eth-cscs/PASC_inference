function [ id ] = get_id_header( C, name )

id = -1;
for i=1:length(C)
    if strcmp(C{i},name) 
        id = i;
    end
end

end

