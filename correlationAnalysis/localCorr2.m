function corrImage = localCorr2(A,B, width)
%LOCALCORR2 Summary of this function goes here
%   Detailed explanation goes here
assert(isequal(size(A), size(B)), "Matrix sizes do not agree")
[rows,cols] = size(A);
corrImage = zeros(size(A));
for row=width+1:rows-width
    for col=width+1:cols-width
        subregionA = A(row-width:row+width,col-width:col+width);
        subregionB = B(row-width:row+width,col-width:col+width);
        corrImage(row,col) = corr2(subregionA,subregionB);
    end
end

