function [F] = load_shortinfo( filename )

M = csvread(filename,1,0);

ns = unique(M(:,1));
sigmas = unique(M(:,2));

Ftemp = zeros(length(sigmas),2);
for i=1:length(sigmas)
    sigma = sigmas(i);

    myabserr_best = M(M(:,2) == sigma,5);
        
    Ftemp(i,1) = sigma;
    Ftemp(i,2) = sum(myabserr_best)/length(myabserr_best);
end

[~,sortidx] = sort(Ftemp(:,1));
F = Ftemp(sortidx,2);

end