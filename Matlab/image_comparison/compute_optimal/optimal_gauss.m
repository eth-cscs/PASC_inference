function [Fout, Ferr, Nbest, Sbest] = optimal_gauss(image_data,image_solution, N, S)

[height, width] = size(image_data);

Ferr = Inf;
Fout = zeros(size(image_data));
Sbest = Inf;
Nbest = Inf;

for s=1:length(S)
    for n=1:length(N)
        Fout_s = image_gauss(image_data, N(n), S(s));
        Ferr_s = mean(reshape(abs(Fout_s - image_solution),width*height,1));

        if Ferr_s < Ferr
            % this solution is better than previous one
            Ferr = Ferr_s;
            Fout = Fout_s;
            Sbest = S(s);
            Nbest = N(n);
        end
    end
end


end