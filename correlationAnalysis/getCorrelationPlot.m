function getCorrelationPlot(TDTmask,BSN,DF)
%CORRELATIONPLOT Plot correlation between BSN and DF using TDTmask as mask

validIndices = find(~isnan(TDTmask));
    validBSN = double(BSN(validIndices));
    validDF = DF(validIndices);
    
    figure()
    loglog(validBSN, validDF,'.')
    hold on
    xlabel("BSN Masked by TDT")
    ylabel("DF Masked by TDT")
    
    %[coefficients, S] = polyfit(validBSN, validDF, 1);
    %validBSNfit = linspace(min(validBSN), max(validBSN), 1000);
    %[validDFfit, delta] = polyval(coefficients , validBSNfit, S);
    %loglog(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
    [correlation_coefficient,P] = corrcoef(validBSN, validDF);%y(randperm(length(y))));
    disp(['Correlation Coefficient: ', num2str(correlation_coefficient(1,2))]);
    disp(['P-value: ', num2str(P(1,2))]);
    title(['Correlation Coefficient: ', num2str(correlation_coefficient(1,2))])
    % moderate if r=0.3
    % https://psychology.emory.edu/clinical/bliwise/Tutorials/SCATTER/scatterplots/effect.htm
end

