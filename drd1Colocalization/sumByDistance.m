function [centers, sums] = sumByDistance(dist, vals, binWidth)
edges = 0:binWidth:100000;
centers = edges(1:end-1) + binWidth/2;
[~,~,bin] = histcounts(dist, edges);
centers = centers(unique(bin))';
sums = accumarray(bin, vals, [], @sum);


end