function [m, CI] = myConfidenceInterval(x, lvl)
    m = mean(x);
    SEM = std(x)/sqrt(length(x));               % Standard Error
    delta = 1-lvl;
    ts = tinv(1-delta/2,length(x)-1);      % T-Score
    CI = ts*SEM;                      % Confidence Intervals
end

