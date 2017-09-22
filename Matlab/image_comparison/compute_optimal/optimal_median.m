function [Fout, Ferr, Nbest] = optimal_median(image_data,image_solution, N)

[height, width] = size(image_data);

Ferr = Inf;
Fout = zeros(size(image_data));
Nbest = Inf;

for n=1:length(N)
    Fout_s = image_median(image_data, N(n));
    Ferr_s = mean(reshape(abs(Fout_s - image_solution),width*height,1));

    if Ferr_s < Ferr
       % this solution is better than previous one
       Ferr = Ferr_s;
       Fout = Fout_s;
       Nbest = N(n);
    end
end


end