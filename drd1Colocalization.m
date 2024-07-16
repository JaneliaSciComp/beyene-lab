%% read in files
close all;
drd1 = imread('/groups/beyene/beyenelab/Imaging Data/Ackerman/06-25 Drd1_tdt_Dapi/DRD1.tif');
dapi = imread('/groups/beyene/beyenelab/Imaging Data/Ackerman/06-25 Drd1_tdt_Dapi/DAPI.tif');
tdt = imread('/groups/beyene/beyenelab/Imaging Data/Ackerman/06-25 Drd1_tdt_Dapi/TDT.tif');

% for now hardcode pixel size
pixelSize = 0.2071607;

% convert to binary using automatic thresholding
drd1bw = imbinarize(imgaussfilt(drd1,2));
dapibw = imbinarize(imgaussfilt(dapi,2));
tdtbw = imbinarize(imgaussfilt(tdt,2), 'adaptive');

imwrite(dapibw, ['tifs' filesep 'DAPIbw.tif'])
imwrite(drd1bw, ['tifs' filesep 'DRD1bw.tif'])
imwrite(tdtbw, ['tifs' filesep 'TDTbw.tif'])

%% Option to do masking of DAPI with DRD1
doMask = true;
% dilate for mask
prefix = '';
if doMask
    prefix = 'Masked_';
    se = strel('disk',10);
    drd1Dilated = imdilate(drd1bw, se);
    imwrite(drd1Dilated,['tifs' filesep 'DRD1mask.tif'])
    dapibw(drd1Dilated) = 0;
    imwrite(dapibw,['tifs' filesep 'DAPIbwmasked.tif'])
    dapi = double(dapi);
    dapi(drd1Dilated) = NaN;
end

%% Preprocess files
% calculate distance transform
drd1dist = bwdist(drd1bw);
dapidist = bwdist(dapibw);
tdtdist = bwdist(tdtbw);

% vectorize
drd1bw = drd1bw(:);
dapibw = dapibw(:);

tdt = tdt(:);
drd1 = drd1(:);
dapi = dapi(:);

drd1dist = drd1dist(:);
dapidist = dapidist(:);
tdtdist = tdtdist(:);



%% distance from non-tdt channels
% group tdt values by bins and take averages
% bin by every 2 pixles
[tdtCentersDapi, tdtMeanDapi, tdtConfidenceLevelDapi] = binByDistance(dapidist, tdt, 2,0.99);
[tdtCentersDrd1, tdtMeanDrd1, tdtConfidenceLevelDrd1 ] = binByDistance(drd1dist, tdt, 2,0.99);

% plot
figure()
hold on;
plotWithShadedError(tdtCentersDapi*pixelSize,tdtMeanDapi, tdtConfidenceLevelDapi, 'b')
plotWithShadedError(tdtCentersDrd1*pixelSize,tdtMeanDrd1, tdtConfidenceLevelDrd1, 'r')
xlabel('Distance from DAPI or DRD1 (\mum)')
ylabel('TDT Intensity')
legend('','DAPI','','DRD1')
print(gcf, '-dtiff', ['tifs' filesep prefix 'TDT_Intensity_vs_Distance.tiff']);

%% distance from tdt
[centersDapi, meanDapi, confidenceLevelDapi] = binByDistance(tdtdist, normalize(double(dapi)), 2,0.99);
[centersDrd1, meanDrd1, confidenceLevelDrd1 ] = binByDistance(tdtdist, normalize(double(drd1)), 2,0.99);
figure()
hold on;
plotWithShadedError(centersDapi*pixelSize,meanDapi, confidenceLevelDapi, 'b')
plotWithShadedError(centersDrd1*pixelSize,meanDrd1, confidenceLevelDrd1, 'r')
xlabel('Distance from TDT (\mum)')
ylabel('Normalized Intensity')
legend('','DAPI','','DRD1')
print(gcf, '-dtiff', ['tifs' filesep prefix 'Intensity_vs_TDT_Distance.tiff']);


%% fraction of pixels
[centersDapi, sumsDapi] = sumByDistance(tdtdist, dapibw/sum(dapibw(:)), 2);
[centersDrd1, sumsDrd1] = sumByDistance(tdtdist, drd1bw/sum(drd1bw(:)), 2);
figure()
hold on;
plot(centersDapi*pixelSize,sumsDapi*100, 'b', 'LineWidth', 2)
plot(centersDrd1*pixelSize,sumsDrd1*100, 'r', 'LineWidth', 2)
xlabel('Distance from TDT (\mum)')
ylabel('Percent Of Segmented Pixels')
legend('DAPI','DRD1')
print(gcf, '-dtiff', ['tifs' filesep prefix 'Percent_Segmented_vs_TDT_Distance.tiff']);

