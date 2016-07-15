function [statData] = spectroNorm(TFR1,TFR2,normType)
%% Normalize data
%% Average time-frequency data across trials
meanValsEvent = squeeze(nanmean(TFR1.powspctrm,1));
meanValsBase = squeeze(nanmean(TFR2.powspctrm,1));
%% Normalize
if (strcmp(normType,'absolute'))
    statData = meanValsEvent - meanValsBase;
elseif (strcmp(normType,'relative'))
    statData = meanValsEvent ./ meanValsBase;
elseif (strcmp(normType,'relchange'))
    statData = (meanValsEvent - meanValsBase) ./ meanValsBase;
end

