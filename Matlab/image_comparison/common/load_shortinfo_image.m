function [F] = load_shortinfo_image( filename, K )

M = csvread(filename,1,0);

ns = unique(M(:,1));
sigmas = unique(M(:,2));

Ftemp = zeros(length(sigmas),2);
for i=1:length(sigmas)
    sigma = sigmas(i);

    myabserr_best = zeros(length(ns),1);
    for j=1:length(ns)
        n = ns(j);
        
        myabserr_best(j) = min(M(M(:,2) == sigma & M(:,1) == n, 8));
    end
    
    Ftemp(i,1) = sigma;
    Ftemp(i,2) = sum(myabserr_best)/length(myabserr_best);
end

[~,sortidx] = sort(Ftemp(:,1));
F = Ftemp(sortidx,2);

end