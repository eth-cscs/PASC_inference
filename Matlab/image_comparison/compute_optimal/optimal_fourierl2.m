function [Fout, Ferr, Sbest, Lambdabest] = optimal_fourierl2(image_data,image_solution, S, Lambda)

[height, width] = size(image_data);

Ferr = Inf;
Fout = zeros(size(image_data));
Sbest = Inf;
Lambdabest = Inf;

for s=1:length(S)
    for l=1:length(Lambda)
        Fout_s = image_fourierl2(image_data, S(s), Lambda(l));
        Ferr_s = mean(reshape(abs(Fout_s - image_solution),width*height,1));

        if Ferr_s < Ferr
            % this solution is better than previous one
            Ferr = Ferr_s;
            Fout = Fout_s;
            Sbest = S(s);
            Lambdabest = Lambda(l);
        end
    end
end


end