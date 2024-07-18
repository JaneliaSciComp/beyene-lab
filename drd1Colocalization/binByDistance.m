function [centers, means, cis] = binByDistance(dist, vals, binWidth, confidenceLevel)

confidenceLevelFunc = @(x) tinv(0.5+confidenceLevel/2, length(x)-1) ...
              * std(x, 0, "omitmissing") / sqrt(sum(~isnan(x)));
stdErrorFunc = @(x) std(x, 0, "omitmissing") / sqrt(sum(~isnan(x)));  % Standard error
edges = 0:binWidth:100000;
centers = edges(1:end-1) + binWidth/2;
[~,~,bin] = histcounts(dist, edges);
centers = centers(unique(bin))';
meanOmitmissing = @(x) mean(x, "omitmissing");
means = accumarray(bin, vals, [], meanOmitmissing);
cis = accumarray(bin, double(vals), [], confidenceLevelFunc);


end