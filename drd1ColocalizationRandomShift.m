% read in file
close all;
drd1 = imread('/groups/beyene/beyenelab/Imaging Data/Ackerman/06-25 Drd1_tdt_Dapi/DRD1.tif');
dapi = imread('/groups/beyene/beyenelab/Imaging Data/Ackerman/06-25 Drd1_tdt_Dapi/DAPI.tif');
tdtOriginal = imread('/groups/beyene/beyenelab/Imaging Data/Ackerman/06-25 Drd1_tdt_Dapi/TDT.tif');




% convert to binary using automatic thresholding
drd1bw = imbinarize(imgaussfilt(drd1,2));
dapibw = imbinarize(imgaussfilt(dapi,2));

% dilate for mask
doMask = false;
prefix = '';
if doMask
    prefix = 'Masked_';
    se = strel('disk',10);
    drd1Dilated = imdilate(drd1bw, se);
    dapibw(drd1Dilated) = 0;
    dapi = double(dapi);
    dapi(drd1Dilated) = NaN;
end

drd1dist = bwdist(drd1bw);
dapidist = bwdist(dapibw);


% vectorize
drd1bw = drd1bw(:);
dapibw = dapibw(:);

drd1 = drd1(:);
dapi = dapi(:);

drd1dist = drd1dist(:);
dapidist = dapidist(:);

numTrials = 100;
binWidth=2;
allCenters = (binWidth/2:binWidth:1000);

allTdtMeanDapi = NaN(numTrials, length(allCenters));
allTdtMeanDrd1 = NaN(numTrials, length(allCenters));

allMeanDapi = NaN(numTrials, length(allCenters));
allMeanDrd1 = NaN(numTrials, length(allCenters));

allSumDapi = NaN(numTrials, length(allCenters));
allSumDrd1 = NaN(numTrials, length(allCenters));


for i=1:numTrials
    tdt = circshift(tdtOriginal, [randi(size(tdtOriginal,1),1,1), randi(size(tdtOriginal,2),1,1)]);
    
    % for now hardcode pixel size
    pixelSize = 0.2071607;
    
    % convert to binary using automatic thresholding
    tdtbw = imbinarize(imgaussfilt(tdt,2), 'adaptive');
    
    % calculate distance transform
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
    [tdtCentersDapi, tdtMeanDapi, tdtConfidenceLevelDapi] = binByDistance(dapidist, tdt, binWidth,0.99);
    [tdtCentersDrd1, tdtMeanDrd1, tdtConfidenceLevelDrd1 ] = binByDistance(drd1dist, tdt, binWidth,0.99);
    % Set the value of allTdtMeanDapi where tdtCentersDapi equals allCenters
    allTdtMeanDapi(i,ismember(allCenters, tdtCentersDapi)) = tdtMeanDapi;
    allTdtMeanDrd1(i,ismember(allCenters, tdtCentersDrd1)) = tdtMeanDrd1;

    
    %% distance from tdt
    [centersDapi, meanDapi, confidenceLevelDapi] = binByDistance(tdtdist, normalize(double(dapi)), binWidth,0.99);
    [centersDrd1, meanDrd1, confidenceLevelDrd1 ] = binByDistance(tdtdist, normalize(double(drd1)), binWidth,0.99);
    allMeanDapi(i,ismember(allCenters, centersDapi)) = meanDapi;
    allMeanDrd1(i,ismember(allCenters, centersDrd1)) = meanDrd1;
    %% fraction of pixels
    [centersDapi, sumDapi] = sumByDistance(tdtdist, dapibw/sum(dapibw(:)), binWidth);
    [centersDrd1, sumDrd1] = sumByDistance(tdtdist, drd1bw/sum(drd1bw(:)), binWidth);
    allSumDapi(i,ismember(allCenters, centersDapi)) = sumDapi;
    allSumDrd1(i,ismember(allCenters, centersDrd1)) = sumDrd1;
end

    % plot
    figure()
    hold on;
    plotWithShadedError(allCenters*pixelSize,mean(allTdtMeanDapi, "omitmissing"), tinv(0.99,numTrials)*std(allTdtMeanDapi), 'b')
    plotWithShadedError(allCenters*pixelSize,mean(allTdtMeanDrd1, "omitmissing"), tinv(0.99,numTrials)*std(allTdtMeanDrd1), 'r')
    legend('','DAPI','','DRD1')
    xlabel('Distance from DAPI or DRD1 (\mum)')
    ylabel('TDT Intensity')
    print(gcf, '-dtiff', ['tifs' filesep prefix 'RandomShifts_TDT_Intensity_vs_Distance.tiff']);

    figure()
    hold on;
    plotWithShadedError(allCenters*pixelSize,mean(allMeanDapi, "omitmissing"), tinv(0.99,numTrials)*std(allMeanDapi), 'b')
    plotWithShadedError(allCenters*pixelSize,mean(allMeanDrd1, "omitmissing"), tinv(0.99,numTrials)*std(allMeanDrd1), 'r')
    xlabel('Distance from TDT (\mum)')
    ylabel('Normalized Intensity')
    legend('','DAPI','','DRD1')
    print(gcf, '-dtiff', ['tifs' filesep prefix 'RandomShifts_Intensity_vs_TDT_Distance.tiff']);

    figure()
    hold on;
    plotWithShadedError(allCenters*pixelSize,mean(allSumDapi, "omitmissing"), tinv(0.99,numTrials)*std(allSumDapi), 'b')
    plotWithShadedError(allCenters*pixelSize,mean(allSumDrd1, "omitmissing"), tinv(0.99,numTrials)*std(allSumDrd1), 'r')
    xlabel('Distance from TDT (\mum)')
    ylabel('Percent Of Segmented Pixels')
    legend('','DAPI','','DRD1')
    print(gcf, '-dtiff', ['tifs' filesep prefix 'RandomShifts_Percent_Segmented_vs_TDT_Distance.tiff']);

    