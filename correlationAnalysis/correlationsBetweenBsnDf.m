close all
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3074624/
% The most obvious advantage of MCC is that it is a more intuitive measure of colocalization than PCC. MCC is also more useful for data that are poorly suited to the simple, linear model that underlies PCC. For example, in the image shown in Fig. 4E, nearly all of the immunolocalized EEA1 appears to occur in compartments associated with GFP-RhoB. This appearance is quantitatively supported by MCC analysis indicating that 85% of the EEA1 localizes to compartments associated with GFP-RhoB. However, PCC measurement indicates a relatively poor association between the probes (PCC = 0.73) due to the fact that, while colocalized, the ratio of the probes varies widely (Fig. 4H).
% Hi David, BSN is from the neurons that express TDT. TDT is the tdtomato channel (these are the dopamine neuron axons) 
% 
%  The ∆F channel is the dopamine activity released from the dopamine neuron axons.
% 
% MAP2 stands for microtube associated protein 2. It is broadly expressed in dendrites and cell bodies of all neurons (but not in axons-- the TDT is dopamine neuron axon and hence doesn't have any MAP2 signal).
%
% Since axons synapse onto soma/dendrites, we use TDT and MAP2 images to visualize those synapse (or contact) sites.
% 
% We want to show (quantitatively) that BSN is enriched in TDT axons that release a lot of dopamine, and expression is low in TDT axons that do not release dopamine.
% 
% Hope this makes sense. Let me know if it doesn't.
% 
close all;

%path = '../../Colocalization-Analysis-Data/fromAbrahamEmail_20240105/TDT-MAP2-DFF-BSN-1.tif';
path = '../../Jan11-2024/';
if contains(path, "20240105")
    TDT = imread(path,1);
    MAP2 = imread(path,2);
    DF = double(imread(path,3));
    BSN = double(imread(path,4));
else
    TDT = imread(path + "TDT.tif");
    MAP2 = imread(path + "MAP2.tif");
    DF = double(imread(path + "DFF_Ninox_32Bit.tif"));
    BSN = double(imread(path + "BSN-Original.tif"));
end
%BSNmat2gray = mat2gray(BSN);

TDTbinarized = imbinarize(TDT);
TDTbinarizedNaN = 1.0*TDTbinarized;
TDTbinarizedNaN(TDTbinarized==0)=NaN;

DFmasked = DF.*TDTbinarizedNaN;
BSNmasked = double(BSN).*TDTbinarizedNaN;

% basic figures
figure(); imagesc(BSNmasked);colormap([hot;0,0,0]);
title('BSN Masked By TDT')
c=colorbar;
ylabel(c,'Intensity')


figure(); imagesc(DFmasked);colormap([0,0,0; jet]);
title('DF Masked By TDT')
c=colorbar;
ylabel(c,'Intensity')

getBSNtimesDFimage(TDTbinarizedNaN,BSN,DF);
%% loglogplot and line of best fit
getCorrelationPlot(TDTbinarizedNaN,BSN,DF);
% moderate if r=0.3
% https://psychology.emory.edu/clinical/bliwise/Tutorials/SCATTER/scatterplots/effect.htm
%% correlation image
correlationWindowRadius = 25;
getCorrelationImage(TDTbinarizedNaN, BSN, DF, correlationWindowRadius)

%% histogram
%figure();
%hist3([x, y],'CdataMode','auto','nbins',[25,25])
%xlabel('BSN')
%ylabel('DF')
%colorbar;
%view(2)