function [z_bin] = customHeatmap(x,y,z, binWidth,sumOrMean, XLabel,YLabel,ZLabel,filename)
%% chatgpt

% Define the bin edges for x and y
row_edges = 0:binWidth:1;
column_edges = 0:binWidth:1;

% Use the discretize function to determine which bin each x and y value belongs to
row_bin = discretize(y, row_edges);
column_bin = discretize(x, column_edges);
sz = [length(row_edges)-1,length(column_edges)-1];

% Use accumarray to average the z values for each bin
if sumOrMean=="mean"
    z_bin = accumarray([row_bin', column_bin'], z, sz, @mean, NaN); %  empty values will be NaN
    %z_sum = accumarray([x_bin', y_bin'], ones(size(z)), sz, @sum);
    %z_bin(z_sum==0) = NaN;
else
    z_bin = accumarray([row_bin', column_bin'], z, sz, @sum, .1);
    z_bin = log10(z_bin);
end

% Generate a 2D plot of the binned z data
[X, Y] = meshgrid(row_edges(1:end-1)+binWidth/2, column_edges(1:end-1)+binWidth/2);
%Z = reshape(z_bin, numel(x_edges)-1, numel(y_edges)-1);
figure();
surf(X, Y, z_bin);
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
axis equal
grid off
view(2);
shading interp
colormap(jet);
c = colorbar;
ylabel(c, ZLabel);

xlabel(XLabel);
ylabel(YLabel);
box on;
% saveas(gcf,"images/"+filename,"epsc")
print('-painters','-depsc2',"images/"+filename)

end

