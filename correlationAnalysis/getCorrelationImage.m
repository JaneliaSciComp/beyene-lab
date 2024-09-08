function getCorrelationImage(TDTmask,BSN,DF,correlationWindowRadius)
%CORRELATIONImage Plot correlation between BSN and DF using TDTmask as mask

[sizeX,sizeY] = size(TDTmask);
correlationImage = NaN*ones(sizeX,sizeY);
pImage = NaN*ones(sizeX,sizeY);

TDTmaskPadded = padarray(TDTmask,[correlationWindowRadius,correlationWindowRadius],NaN,'both');
BSNpadded = padarray(BSN,[correlationWindowRadius,correlationWindowRadius],NaN,'both');
DFpadded = padarray(DF,[correlationWindowRadius,correlationWindowRadius],NaN,'both');

% generate circular crop
circleX = [];
circleY = [];
for x = 0:2*correlationWindowRadius
    dx = correlationWindowRadius - x;
    for y= 0:2*correlationWindowRadius
        dy = correlationWindowRadius - y;
        if dx^2 + dy^2 <=correlationWindowRadius^2
            circleX = [circleX; x];
            circleY = [circleY; y];
        end
    end
end

parfor windowCenterX = 1:sizeX
    for windowCenterY = 1:sizeY
        if ~isnan(TDTmaskPadded(correlationWindowRadius+windowCenterX, correlationWindowRadius+windowCenterY))
            ind = sub2ind(size(TDTmaskPadded),circleX + windowCenterX, circleY + windowCenterY);
            validIndices = ind(~isnan(TDTmaskPadded(ind)));
    
            if length(validIndices)>1%=length(ind)*0.5
                windowedBSNpadded = BSNpadded(validIndices);
                windowedDFpadded = DFpadded(validIndices);
                [correlation_coefficient,P] = corrcoef(windowedBSNpadded, windowedDFpadded);
                correlationImage(windowCenterX,windowCenterY)= correlation_coefficient(1,2);
                pImage(windowCenterX, windowCenterY) = P(1,2);
            end
        end
    end
end
figure();
imagesc(correlationImage)
colormap([0,0,0;flip(redgreencmap,1)])
c = colorbar;
ylabel(c,'Correlation');
title(['Correlation Image Masked (Window Radius = ', num2str(correlationWindowRadius),')'])
end

