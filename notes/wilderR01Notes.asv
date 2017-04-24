%% Binge vs. Rest vs. Binge-Rest Prediction Comparison
load('C:\Users\Lucas\Desktop\GreenLab\data\wilderR01\compData.mat')
% Variables
% data = 2x2 cell where first row is data from baseline condition and
%   second is from food deprivation 24 hours condition; first column is
%   binge data and second is rest
% bingeSize = 12x2 matrix where row = animal and column = condition (1:
%   base; 2: food dep24)
% subIndex = subset indices becuase of double recordings at baseline.
%% Get binge change; both absolute and relative (normalized)
% Absolute (grams)
bingeSizeDiffAbs = bingeSize(:,2) - bingeSize(:,1);
% Relative (percent change)
bingeSizeDiffRel = bingeSizeDiffAbs./bingeSize(:,1);
%% Get ephys change from baseline tables [dep24-base]; both absolute and relative
% Subset base data
data{1,1} = data{1,1}(subInd,:);
data{1,2} = data{1,2}(subInd,:);
% Absolute differences
% depb-baseb
bingeDiffAbs = data{2,1}-data{1,1};
%depr-baser
restDiffAbs = data{2,2}-data{1,2};
%[depb-depr] - [baseb-baser]
bingeNormDiffAbs = (data{2,1}-data{2,2}) - (data{1,1}-data{1,2});
% Stack for cycling
absoluteDat = cat(3,bingeDiffAbs,restDiffAbs,bingeNormDiffAbs);
% Relative
bingeDiffRel = bingeDiffAbs./data{1,1};
restDiffRel = restDiffAbs./data{2,2};
bingeNormDiffRel = bingeNormDiffAbs./(data{1,1}-data{1,2});
% Stack for cycling
relativeDat = cat(3,bingeDiffRel,restDiffRel,bingeNormDiffRel);
%% Run data through glmnet and get error distributions
% Set up config for lassoNet
cfg = lassoNetCfg([],'n','y','n',100,'1se');
% Cycle through data sets
for ii = 1:size(absoluteDat,3)
    % Absolute data
    [~,allAbsLambda{ii},allAbsBeta{ii},~,~,~] = lassoNet(absoluteDat(:,:,ii),bingeSizeDiffAbs,'gaussian','mse',1,4,1,cfg);
    allAbsErr(ii,:) = allAbsLambda{ii}{1,1}.allErr';
    % Relative data
%     [~,allRelLambda{ii},~,~,~,~] = lassoNet(relativeDat(:,:,ii),bingeSizeDiffRel,'gaussian','mse',0:0.1:1,4,1,cfg);
%     allRelErr(ii,:) = allRelLambda{ii}{1,1}.allErr';
end
%% Plot distributions
figure
hold on
for ii = 1:size(allAbsErr)
    histogram(allAbsErr(ii,:),'Normalization','probability','BinWidth',2)
end
legend({'Binge','Rest','Binge-Rest'},'Location','northeast')
title('Error Distributions')
xlabel('Error (mse)')
ylabel('Probability')
figure
hold on
for ii = 1:size(allAbsErr)
    ecdf(allAbsErr(ii,:))
end
title('Error Cumulative Density')
xlabel('Error (mse)')
ylabel('Cumulative Density')
legend({'Binge','Rest','Binge-Rest'},'Location','southeast')
% figure
% hold on
% for ii = 1:size(allRelErr)
%    histogram(allRelErr(ii,:),'Normalization','probability','BinWidth',0.5) 
% end
% Plot ecdf

