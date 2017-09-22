function [Fout, Ferr, Sbest] = optimal_fourier(image_data,image_solution, S)

% f0 is a signal with noise which has to be filtered
%f0 = )';
%if mod(size(f0,2),2) == 1
%   f0 = f0(:,1:end-1);
%end
        
% compute fourier for different sizes of window "s".
%S=[5 20 30 40 60 80];

[height, width] = size(image_data);

Ferr = Inf;
Fout = zeros(size(image_data));
Sbest = Inf;

for s=1:length(S)
    Fout_s = image_fourier(image_data, S(s));
    Ferr_s = mean(reshape(abs(Fout_s - image_solution),width*height,1));

    if Ferr_s < Ferr
       % this solution is better than previous one
       Ferr = Ferr_s;
       Fout = Fout_s;
       Sbest = S(s);
    end
end


end