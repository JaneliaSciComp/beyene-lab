function plotWithShadedError(x,y,e, color)
if size(x,1)~=1
    x = x';
end
if size(y,1)~=1
    y = y';
end
if size(e,1)~=1
    e=e';
end
xs = [x'; flip(x')];
ys = [y'-e'; flip(y'+e')];
invalid = isnan(ys);
xs(invalid) = [];
ys(invalid) = [];


%dontFill = ys(:,1) ==0 & ys(:,2)==0;
%xs(dontFill,:) = [];
%ys(dontFill,:) = [];
%hold on;
fill(xs, ys,color, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;

plot(x, y, 'color',color, 'LineWidth', .75);
set(gca,'Color','k')

end