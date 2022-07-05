for ii = 1:216
    [~,p(ii),~,stats] = ttest2(fedData(y==1,ii),fedData(y==0,ii));
    t(ii) = stats.tstat;
end
sfc = 0.8*max(abs(t))*sqrt(1/24+1/22);
%%
load('F:\irdmRound2\interData.mat','effect')
% Load last baseline data
load('F:\irdmRound2\lastBaseNorm2.mat')
% Remove empty data from both
theseInd = logicFind(1,~cellfun(@isempty,lastBase(:,1)),'==');
thisLastBase = lastBase(theseInd,:);
% Round effect to get single value
thisEffect = round(cellfun(@mean,effect));
thisEffect = thisEffect(theseInd,:);
% Average across time for each baseline
mBase = cellfun(@(x) mean(x,'omitnan'),thisLastBase,'UniformOutput',false);
% Core
thisX = cell(1);
thisY = cell(1);
c = 1;
for ii = 1:size(thisEffect,1)
    if ~isnan(thisEffect(ii,2))
        thisX{c} = cat(1,mBase{ii,1:3});
        thisY{c} = repmat(thisEffect(ii,2),3,1);
        c = c+1;
    end
end
theseX = cat(1,thisX{:});
theseY = cat(1,thisY{:});
for ii = 1:216
    [~,p(ii),~,stats] = ttest2(theseX(theseY==1,ii),theseX(theseY==0,ii));
    t(ii) = stats.tstat;
end
sfc = 0.8*max(abs(t))*sqrt(1/sum(theseY)+1/sum(theseY==0));
% sfc = 0.8*max(abs(t))*sqrt(1/7+1/9);