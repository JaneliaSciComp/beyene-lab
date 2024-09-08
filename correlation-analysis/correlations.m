close all
clear all
%% Distance transform of MAP to either "DF masked by TDT" or "TDT"
%path = '../../Colocalization-Analysis-Data/type-2-data/CB-2/Composite_MAP2_TDT_GAD1_DAPI_DFF.tif';    
path = '../../Colocalization-Analysis-Data/fromAbrahamEmail_20240105/TDT-MAP2-DFF-BSN-1.tif'
%path = '../../Colocalization-Analysis-Data/type-1-data/AB-1/Composite_MAP2_TDT_DFF-Example1.tif';    
%path = '../../Colocalization-Analysis-Data/type-1-data/AB-2/Whole_FOV_Example/TDT-MAP2-DFF-Composite.tif' ;

% As for the standardization of the file names, here is what I propose your codes to assume:
% All channels will be saved as separate files and placed in the root folder.
% Channel 1 will have the following name: C1_DFF
% Channel 2 will have the following name: C2_TDT
% Channel 3 will have the following name: C3_"P3"
% Channel 4 will have the following name: C4_"P4"
% Channel 5 will have the following name: C4_"P5"
% Proteins in channels 3, 4 and 5 will vary from one specimen to another (they may not even be proteins, it could be DAPI, for example). In this case, a Bassoon signal in channel 3 will be named as such: C3_BSN. Similarly, MUNC13 in channel 4 will be called: C4_MUNC13 etc

% desired order:  TDT,MAP2,DF
if contains(path, 'fromAbrahamEmail_20240105')
    channelOrder = [1,2,3];
elseif contains(path,'TDT-MAP2-DFF')
    channelOrder = [1,3,2];
elseif contains(path, 'MAP2_TDT_DFF')
    channelOrder = [2,3,1];
elseif contains(path, 'MAP2_TDT_GAD1_DAPI_DFF')
    channelOrder = [2,3,5];
end
TDT = mat2gray(imread(path,channelOrder(1)));
MAP2 = mat2gray(imread(path,channelOrder(2)));
DF = mat2gray(imread(path,channelOrder(3)));

corr_single = corr2(DF,TDT);
corr_combined = corr2(DF,TDT.*MAP2);
fprintf('%f %f\n',corr_single,corr_combined);

corr_single_dt = corr2(DF,1-mat2gray(bwdist(imbinarize(TDT))));
corr_combined_dt = corr2(DF,1-mat2gray(bwdist(imbinarize(TDT).*imbinarize(MAP2))));
fprintf('%f %f\n',corr_single_dt,corr_combined_dt);

corr_single_dt = corr2(DF,1-mat2gray(bwdist(imbinarize(TDT))));
corr_combined_dt = corr2(DF,1-mat2gray(bwdist(imbinarize(TDT).*imbinarize(MAP2))));
fprintf('%f %f\n',corr_single_dt,corr_combined_dt);

figure();

%DF_TDT = mat2gray(imbinarize(DF.*(imbinarize(TDT)>0)*1.0));

%TODO: maybe try thresholding

DF_TDT = mat2gray(imregionalmax(DF).*(imbinarize(TDT,'adaptive')>0)*1.0);
%DF_TDT = mat2gray(imextendedmax(DF,0.01).*(imbinarize(TDT)>0)*1.0);
fontsize = 20;
subplot(1,2,1)
% TODO: maybe more stringent for binarization
image=zeros(size(TDT,1),size(TDT,2),3); 
image(:,:,1) = imbinarize(TDT,'adaptive');
image(:,:,2) = DF_TDT;
MAP2_dist = bwdist(imbinarize(MAP2));
image(:,:,3)=mat2gray(MAP2_dist);
imshow(image)
title('\color{red}TDT, \color{green}DF masked by TDT, \color{blue}MAP2 Distance')
set(gca,'FontSize',fontsize)


subplot(1,2,2)
axis square
% set(gca,'YScale','log')

DF_TDT(DF_TDT==0) = NaN;

TDT_binarized = mat2gray(imbinarize(TDT));
TDT_binarized(TDT_binarized==0) = NaN;

MAP2_dist_to_DF_TDT = MAP2_dist.*DF_TDT;

MAP2_dist_to_TDT = MAP2_dist.*TDT_binarized;
binWidth = 2;
binEdges = 0:binWidth:max([MAP2_dist_to_DF_TDT(:); MAP2_dist_to_TDT(:)])+binWidth;
hold on

MAP2_dist_to_DF_TDT = MAP2_dist_to_DF_TDT(~isnan(MAP2_dist_to_DF_TDT));
MAP2_dist_to_TDT = MAP2_dist_to_TDT(~isnan(MAP2_dist_to_TDT));
histogram(MAP2_dist_to_TDT,'BinEdges',binEdges,'Normalization','probability','FaceColor','k')
histogram(MAP2_dist_to_DF_TDT,'BinEdges',binEdges,'Normalization','probability','FaceColor','red');
legend(sprintf('TDT, mean = %f',mean(MAP2_dist_to_TDT)),sprintf('DF masked by TDT, mean = %f',mean(MAP2_dist_to_DF_TDT)))
ylabel('Fraction')
xlabel('Distance To MAP2')
set(gca,'FontSize',fontsize)

% Test if the values are from the same distribution
[h,p,k]=kstest2(MAP2_dist_to_TDT,MAP2_dist_to_DF_TDT,'tail','smaller');

fprintf('Mean distance of MAP2 to binarized TDT: %f \n', mean(MAP2_dist_to_TDT));
fprintf('Mean distance of MAP2 to binarized DF masked by TDT: %f \n', mean(MAP2_dist_to_DF_TDT));

%% 1. Spatial “clustering” of hotspots:
% There is clustering of release sites in each field of view of imaging. We believe this cannot happen by chance. If you take the DF+ TDT boutons, they appear to be "congregated" and not "spread apart”. DGA: TDT dopamine neurons, coculture them with other nerves (the nucleii, mostly non dopamanergic).
% How can we show this?
% Compute “inter-pixel distances” of DF+ TDT boutons and compare that with “inter-pixel distances” of randomly selected TDT boutons. Choose the TDT boutons such that they are “boutons” (intensity > 2x average could be a good metric, or use thresholding). DGA: calculate if they are clustered more in a 2D sense, NOT along the process itself. If boutons are stochastic, then why we see clustering. POSTSYNAPTIC PARTNERS. Some nuclei may come from neurons, some may be from non-neuronal which may be why some nuclei aren’t displaying DF. Could be different z or may not be neuron

% TODO: Good way of identifying boutons. May use computational approaches, or ML: DeepBouton: Automated
% Identification of Single-Neuron Axonal Boutons at the Brain-Wide Scale,


% try without finding boutons first, just binarize imregionalmax and use those
% pixels as before. create a 2d image heatmap where each axis is the DF
% for the voxel and the heatmap is the distance.

% really gotta do maximum within a region

close all;

TDTbinarized = imbinarize(TDT);
boutonPoints = imregionalmax(imgaussfilt(TDT,5)).*TDTbinarized;
%why this
[boutonRow,boutonColumn,~] = find(imregionalmax(boutonPoints)~=0);

figure();
imshow(TDTbinarized);
hold on;
plot(boutonColumn,boutonRow, 'rx', 'MarkerSize', 5, 'LineWidth', 1);
print('-painters','-depsc2',"images/boutons")

numPoints = length(boutonRow);
distances = pdist([boutonRow,boutonColumn]);
iDF = zeros(size(distances));
jDF = zeros(size(distances));
numPairs = length(distances);
for idx=1:numPairs
    [i,j] = triind2sub(numPoints,idx);
    iDF(idx) = DF(boutonRow(i),boutonColumn(i));
    jDF(idx) = DF(boutonRow(j),boutonColumn(j));
end
% for symmetry
iDF_ = [iDF,jDF];
jDF = [jDF,iDF];
iDF = iDF_;

distances = [distances, distances];

% twicecuz it fails to colormap first time
Z = customHeatmap(iDF,jDF,distances,0.1,'mean',"Bouton {\Delta}F (Normalized)","Bouton {\Delta}F (Normalized)","Average Distance (Pixels)","bouton-bouton-distances");
Z = customHeatmap(iDF,jDF,distances,0.1,'mean',"Bouton {\Delta}F (Normalized)","Bouton {\Delta}F (Normalized)","Average Distance (Pixels)","bouton-bouton-distances");


% h = figure();
% imshow(mat2gray(Z));
% % Apply the colormap and show the colorbar
% colormap(h, jet);
% colorbar;

% Restrict to only distance to nearest,will be asymmetric
iDFs = [];
jDFs = [];
distances = [];
for i = 1:numPoints
    iDF = DF(boutonRow(i),boutonColumn(i));
    jBins = ones(10,1)*Inf;
    currentDistances = sqrt((boutonRow(i) - boutonRow).^2 + (boutonColumn(i)-boutonColumn).^2);
    validIDs = find(currentDistances>0); % 0 is current point
    for idx = 1:length(validIDs)
        validID = validIDs(idx);
        jDF = DF(boutonRow(validID),boutonColumn(validID));
        jBin = min(floor(jDF/0.1)+1,10);
        if currentDistances(validID)<jBins(jBin)
            jBins(jBin) = currentDistances(validID);
        end
    end
    for jBin = 1:length(jBins)
        currentDistance = jBins(jBin);
        if currentDistance<Inf
            iDFs = [iDFs;iDF];
            jDFs = [jDFs;(jBin-1)/10+0.05]; %put it in terms of DF again since will be binning during plotting anyway
            distances = [distances; jBins(jBin)];
        end
    end
end
Z = customHeatmap(iDFs',jDFs',distances',0.1,'mean',"Bouton 1 {\Delta}F (Normalized)","Bouton 2 {\Delta}F (Normalized)","Min Distance (Pixels)","bouton-bouton-min-distances");

% Just care about if the hotspots are clumped together
TDTbinarizedNaN = 1.0*TDTbinarized;
TDTbinarizedNaN(TDTbinarized==0)=NaN;

DFmasked = DF.*TDTbinarizedNaN;
DFthreshold = 2*median(DFmasked(:),"omitnan");
DFhotspots = DFmasked>DFthreshold;
DFnonhotspots = DFmasked<=DFthreshold;

figure();
rgb_image = zeros(size(DFhotspots,1),size(DFhotspots,2),3);
rgb_image(:,:,2) = DFmasked.*DFhotspots;
rgb_image(:,:,1) = DFmasked.*DFnonhotspots;
imshow(rgb_image)
print('-painters','-depsc2',"images/hotspots-and-nonhotspots")

figure();
%TODO: also do distances
hold on;
binSize = 20;
maxBinEdge = sqrt(size(TDT,1)^2 + size(TDT,2)^2)+10;
binEdges = (0:binSize:maxBinEdge);
for t = {TDTbinarized,DFhotspots,DFnonhotspots}
    [r,c]=find(t{:}==1);
    fprintf("num points %d \n",length(r));
    d = pdist([r,c]);
    [ct,edges] = histcounts(d,binEdges,'Normalization', 'probability');
    binCenters = (edges(1:end-1)+edges(2:end))/2;
    %[f,xi] = ksdensity(pdist([r,c]));
    plot(binCenters,ct,'linewidth',3);
end
xlabel("Pairwise Pixel Distance (pixels)")
ylabel("Probability")
legend({'All TDT Masked Voxels';'TDT Masked Voxels at DF Hotspots';'TDT Masked Voxels at Non hotspots'})
print('-painters','-depsc2',"images/hotspot-clustering")


%% 2. Clustered hotspots are in close proximity to a postsynaptic marker:
% We believe that these clusters happen to be in close proximity to a MAP2, DAPI, GFAP , GAD1 or VGLUT1 signal or some other postsynaptic process marker 
% 
% How can we show this?
% Compute distance transforms of DF+ TDT boutons vs. randomly selected boutons relative to said postsynaptic process marker (MAP2, DAPI or GFAP or GAD1 or VGLUT1)

TDTbinarized = imbinarize(TDT);
boutonPoints = imregionalmax(imgaussfilt(TDT,5)).*imbinarize(TDT);
figure();
imshow(cat(3,TDTbinarized,boutonPoints,zeros(size(TDT))))
[boutonRow,boutonColumn,~] = find(imregionalmax(boutonPoints)~=0);
MAP2_dist = mat2gray(bwdist(imbinarize(MAP2)));
DFs = zeros([1,length(boutonRow)]);
distances = zeros([1, length(boutonRow)]);

for idx=1:length(boutonRow)
    distances(idx) = MAP2_dist(boutonRow(idx),boutonColumn(idx));
    DFs(idx) = DF(boutonRow(idx),boutonColumn(idx));
end

% Issue here is in part that most of the stuff is close to MAP2
Z = customHeatmap(DFs,distances,ones([1,length(boutonRow)]),0.1,'sum','Bouton {\Delta}F (Normalized)','Bouton Distance to MAP2 (Normalized)','log(count)','bouton-MAP2-distances');

% try overall TDT
[TDTRow,TDTColumn] = find(TDTbinarized==1);
distances = zeros([1, length(TDTRow)]);
DFs = zeros([1,length(TDTRow)]);
for idx=1:length(TDTRow)
    distances(idx) = MAP2_dist(TDTRow(idx),TDTColumn(idx));
    DFs(idx) = DF(TDTRow(idx),TDTColumn(idx));
end

Z = customHeatmap(DFs,distances,ones([1,length(TDTRow)]),0.1,'sum','TDT {\Delta}F (Normalized)','TDT Distance to MAP2 (Normalized)','log(count)','TDT-MAP2-distances');

% try again, but don't do heatmap use kernel density
%h=figure();
%hold on;
binWidth = 0.1;
binCenters = binWidth/2:binWidth:1;
binEdges = 0:binWidth/2:1;

for wholeOrBouton = ["TDT {\Delta}F (Normalized)","Bouton {\Delta}F (Normalized)"]
    currentDF = DF;
    %set to -1 outside of mask
    if wholeOrBouton=="TDT {\Delta}F (Normalized)"
        currentDF(TDTbinarized==0)=-1;
        outputName = "TDT-MAP2-distances-errorbars";
    else
        mask = ones(size(currentDF));
        for i=1:length(boutonRow)
            mask(boutonRow(i),boutonColumn(i)) = 0;
        end
        currentDF(mask==1)=-1;
        outputName = "bouton-MAP2-distances-errorbars";
    end
    means = [];
    CIs = [];
    for startBinEdge=0:binWidth:1-binWidth
        endBinEdge = startBinEdge+binWidth;
        if startBinEdge==1-binWidth
            distances = MAP2_dist(currentDF>=startBinEdge & currentDF<=endBinEdge);
        else
            distances = MAP2_dist(currentDF>=startBinEdge & currentDF<endBinEdge);
        end
        [m,CI] = myConfidenceInterval(distances,0.95);
        means = [means;m];
        CIs = [CIs;CI];
    end
    figure();
    errorbar(binCenters,means,CIs)
    xlabel(wholeOrBouton);
    ylabel("MAP2 Distance (Normalized)");
    print('-painters','-depsc2',"images/"+outputName)
end

%% 3. Compute correlation coefficients:

% Boutons are where releases of dopamine occur, producing hotspots of activity in the DFF image. We would like to examine the molecular machinery (protein) content of each bouton and compute how well a give protein’s expression correlates with DFF. 
% 
% How can we show this?
% 
% There are many proteins we could stain for. Let’s take the following as an example:
% 
% 1. Stain for Bassoon and MUNC13 (two proteins involved in release)
% 2. Produce a spatial map of correlation coefficients between Bassoon and DFF signal
% 3. Produce a spatial map of correlation coefficients between MUNC13 and DFF signal
% 4. Produce a spatial map of correlation coefficients between Bassoon and MUNC13.
% 5. Show that as Bassoon and MUNC13 overlap, their correlation with DFF improves
% Repeat such analysis with various protein stains.
% 
% The ideas here is to show that as DFF signal correlates better with multiprotein signal vs. just one protein signal.
close all;
path = '../../Colocalization-Analysis-Data/type-2-data/CB-1/1201Dish#9/';    

TDT = mat2gray(imread(path+"TDT.tif"));
DF = mat2gray(imread(path+"DFF.tif"));
BSN = mat2gray(imread(path+"BSN.tif"));
VGLUT2 = mat2gray(imread(path+"VGLUT2.tif"));

%CROP
TDT = TDT(1:1683,1:2000);
DF = DF(1:1683,1:2000);
BSN = BSN(1:1683,1:2000);
VGLUT2 = VGLUT2(1:1683,1:2000);

%do local corr2
TDTmasked = NaN*ones(size(TDT));
DFmasked = NaN*ones(size(DF));
TDTmasked(imbinarize(TDT)) = TDT(imbinarize(TDT));
DFmasked(imbinarize(DF)) = DF(imbinarize(DF));

% figure();
% localCorr2Image = localCorr2(BSN,DF,10);
% imshow(localCorr2Image.*(imbinarize(TDT)+imbinarize(DF)));

% try structural similarity
figure();
[ssimval,ssimmap] = ssim(DF,TDT);
imshow(ssimmap.*imbinarize(TDT));

figure();
imshow(cat(3,DF,imbinarize(BSN),zeros(size(TDTmasked))))

% histograms of DF values
% this shows that if there is for example BSN, is there high DF
figure()
hold on
edges = 0:.01:1;
binarized = {imbinarize(TDT),imbinarize(BSN),imbinarize(VGLUT2),imbinarize(BSN).*imbinarize(VGLUT2)==1};
for idx=1:length(binarized)
    histogram(DF(binarized{idx}),'BinEdges',edges,'Normalization','Probability');
end
legend(["TDT","BSN","VGLUT2","BSN+VGLUT2"])
%close all
figure();
hold on;
for idx=1:length(binarized)
    [f,xi] = ksdensity(DF(binarized{idx})); 
    plot(xi,f,'linewidth',3);
end
legend(["TDT","BSN","VGLUT2","BSN+VGLUT2"])
xlabel("{\Delta}F")
ylabel("Probability")
print('-painters','-depsc2',"images/combined")
pause(1)

%%for =[imbinarize(TDT),imbinarize(BSN),imbinarize(VGLUT2),imbinarize(BSN).*imbinarize(VGLUT2)==1]
%[f,xi] = ksdensity(x); 
%plot(xi,f);
%end

% try again, but we want to see if there is a high DF is there for example
% BSN
% one axis will be DF,the other axis will be BSN
% or can do one axis will be DF, the other is the mean BSN intensity
% or can do "what is the distance from DF values 0-.1 to BSN. if there is correlation"

% see what the distance is from eg. BSN to DF values of different ranges
figure();
hold on
binWidth = 0.1;
binCenters = binWidth/2 : binWidth : 1;
binarized = {imbinarize(TDT),imbinarize(BSN),imbinarize(VGLUT2),imbinarize(BSN).*imbinarize(VGLUT2)==1};
for idx=1:length(binarized)
    means = [];
    CIs = [];
    currentDist = bwdist(binarized{idx});
    for startBinEdge=0:binWidth:1-binWidth
        if startBinEdge==1-binWidth
            distances = currentDist(DF>=startBinEdge & DF<=startBinEdge+binWidth);
        else
            distances = currentDist(DF>=startBinEdge & DF<startBinEdge+binWidth);
        end
        [m,CI] = myConfidenceInterval(distances,.95);
        means = [means; m];
        CIs = [CIs; CI];
    end
    errorbar(binCenters:.1:1,means,CIs,'linewidth',3);
end
legend(["TDT","BSN","VGLUT2","BSN+VGLUT2"])
xlabel('{\Delta}F (Normalized)')
ylabel('Distance (Pixels)')
print('-painters','-depsc2',"images/combined-errorbars")

% use regional maxs
regionalMax = imregionalmax(DF);

%try correlation as well


