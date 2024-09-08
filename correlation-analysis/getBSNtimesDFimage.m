function getBSNtimesDFimage(TDTmask, BSN, DF)
%SHOWBSNTIMESDFIMAGE Multiply and show BSN times DF with TDTmask as mask
BSN(isnan(TDTmask))=NaN;
DF(isnan(TDTmask))=NaN;
figure(); imagesc(mat2gray(BSN).*mat2gray(DF));colormap([jet;0,0,0]);
title('Normalized BSN Masked x Normalized DF Masked')
colorbar

end