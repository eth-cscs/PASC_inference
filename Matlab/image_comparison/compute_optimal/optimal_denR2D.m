function [Fout, Ferr, Tbest] = optimal_denR2D(image_data,image_solution, T)

[height, width] = size(image_data);

Ferr = Inf;
Fout = zeros(size(image_data));
Tbest = Inf;

for t=1:length(T)
    Fout_s = denR2D(image_data, T(t));
    Ferr_s = mean(reshape(abs(Fout_s - image_solution),width*height,1));

    if Ferr_s < Ferr
       % this solution is better than previous one
       Ferr = Ferr_s;
       Fout = Fout_s;
       Tbest = T(t);
    end
end


end