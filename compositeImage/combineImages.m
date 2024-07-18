function [rgbImageF, rgbImageDeltaFoverF] = combineImages(F,deltaFoverF)
    numColors = 256;
    redColormap = [linspace(0, 1, numColors)', zeros(numColors, 1), zeros(numColors, 1)];

    fireLUT = readmatrix('fireLUT.csv');
    fireLUT = fireLUT(:,2:end)/255.0;

    indexedImage = gray2ind(double(F)/16464, size(redColormap, 1));
    rgbImageF = ind2rgb(indexedImage, redColormap);

    indexedImage = gray2ind((double(deltaFoverF)-100)/(600-100),size(fireLUT, 1));
    rgbImageDeltaFoverF = ind2rgb(indexedImage, fireLUT);

end