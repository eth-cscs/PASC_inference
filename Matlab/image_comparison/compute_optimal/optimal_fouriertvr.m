function [Fout, Ferr, Sbest, Lambdabest, Epsilonbest] = optimal_fouriertvr(image_data,image_solution, S, Lambda, Epsilon)

[height, width] = size(image_data);

Ferr = Inf;
Fout = zeros(size(image_data));
Sbest = Inf;
Lambdabest = Inf;
Epsilonbest = Inf;

for s=1:length(S)
    for l=1:length(Lambda)
        for ee = 1:length(Epsilon)
            Fout_s = image_fouriertvrbb(image_data, S(s), Lambda(l), Epsilon(ee), 100);
            Ferr_s = mean(reshape(abs(Fout_s - image_solution),width*height,1));

            if Ferr_s < Ferr
                % this solution is better than previous one
                Ferr = Ferr_s;
                Fout = Fout_s;
                Sbest = S(s);
                Lambdabest = Lambda(l);
                Epsilonbest = Epsilon(ee);
            end
        end
    end
end


end