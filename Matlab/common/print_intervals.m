function print_intervals(x, value)
% print intervals where value \in x in pretty form

last = value-1; % ~= value
mybegin = 0;
myend = 0;

for i = 1:length(x)
    if x(i) == value && last ~= value
        mybegin = i;
        last = value;
    else
        if x(i) ~= value && last == value
            myend = i-1;
            disp(['interval: [' num2str(mybegin) ', ' num2str(myend) ']'])
        end
        last = x(i);
    end
    
end

% last value is not checked using previous algorithm
if last == value
   myend = length(x);
   disp(['interval: [' num2str(mybegin) ', ' num2str(myend) ']'])
end

end

