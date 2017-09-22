function [F, Freduces] = load_shortinfo_image_fem( filename )

M = csvread(filename,1,0);

ns = unique(M(:,1));
fem_reduces = unique(M(:,2));

Ftemp = zeros(length(fem_reduces),2);
for i=1:length(fem_reduces)
    fem_reduce = fem_reduces(i);

    myabserr_best = zeros(length(ns),1);
    for j=1:length(ns)
        n = ns(j);

        sum(M(:,2) == fem_reduce & M(:,1) == n)
        
        myabserr_best(j) = min(M(M(:,2) == fem_reduce & M(:,1) == n, 8));
    end
    
    Ftemp(i,1) = fem_reduce;
    Ftemp(i,2) = sum(myabserr_best)/length(myabserr_best);
end

[~,sortidx] = sort(Ftemp(:,1));
F = Ftemp(sortidx,2);
Freduces = fem_reduces(sortidx);

end