%% NeruoDay 2017 Poster Figures
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\4conditionBingeSize.mat')
%% Voracity between groups
figure
hold on
plot(1,voracity(:,2),'.k','MarkerSize',10)
plot(1,nanmean(voracity(:,2)),'rs','MarkerSize',5,'MarkerFaceColor','r')
plot(1.5,voracity(:,4),'.','Color',[0.5 0.5 0.5],'MarkerSize',10)
plot(1.5,nanmean(voracity(:,4)),'rs','MarkerSize',5,'MarkerFaceColor','r')
xlim([0.5 2])
set(gca,'XTick',[1,1.5],'XTickLabel',{'PF','Chow'})
title('Voracity between Groups')
ylabel('Voracity (gm/sec)')
[h,p] = ttest2(voracity(:,2),voracity(:,4));
%% Voracity between base and deps
figure
hold on
plot(1,voracity(:,1),'.k','MarkerSize',10)
plot(1,nanmean(voracity(:,1)),'rs','MarkerSize',5,'MarkerFaceColor','r')
plot(1.5,voracity(:,2),'.k','MarkerSize',10)
plot(1.5,nanmean(voracity(:,2)),'rs','MarkerSize',5,'MarkerFaceColor','r')
plot(2,voracity(:,3),'.k','MarkerSize',10)
plot(2,nanmean(voracity(:,3)),'rs','MarkerSize',5,'MarkerFaceColor','r')
xlim([0.5 2.5])
set(gca,'XTick',[1,1.5,2],'XTickLabel',{'Base','Dep24','Dep48'})
title('Voracity between Groups')
ylabel('Voracity (gm/sec)')
%% Calories eaten
figure
hold on
plot(1,bingeCal(:,2),'.k','MarkerSize',10)
plot(1,nanmean(bingeCal(:,2)),'rs','MarkerSize',5,'MarkerFaceColor','r')
plot(1.5,bingeCal(:,4),'.','Color',[0.5 0.5 0.5],'MarkerSize',10)
plot(1.5,nanmean(bingeCal(:,4)),'rs','MarkerSize',5,'MarkerFaceColor','r')
xlim([0.5 2])
set(gca,'XTick',[1,1.5],'XTickLabel',{'PF','Chow'})
title('Calories Consumed between Groups')
ylabel('kCalories')
[h,p] = ttest2(bingeCal(:,2),bingeCal(:,4));
sigstar([1,1.5],p)
%% Representative Binge size change with dep
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\4conditionBingeSize.mat')
% Normalize calories by baseline (column 1)
normCal = (bingeCal-repmat(bingeCal(:,1),1,4))./repmat(bingeCal(:,1),1,4); 
figure
hold on
plot([1:3],normCal(5,1:3).*100)
plot([1:3],normCal(1,1:3).*100)
plot([1:3],normCal(7,1:3).*100)
set(gca,'XTick',(1:3),'XTickLabel',{'Base','Dep24','Dep48'})
ylabel('% \Delta in Binge Size')
title('Representative Binge Size Change Groups')
%% All animal t-tests
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\palatableData.mat')
depDiff = (data{1,1}-data{1,2});
chowDiff = (data{2,1}-data{2,2});

[h,p,pAdj] = bulkT(depDiff,0,0,'fdr');
% [h,p,pAdj] = bulkT(data{1,1},data{1,2},1,'fdr');
[hChow,pChow,pAdjChow] = bulkT(chowDiff,0,0,'fdr');
% [h2,p2,pAdj2] = bulkT(data{2,1},data{2,2},1,'fdr');
stripPlot(depDiff,[],pAdj,[],'Palatable')
stripPlot(chowDiff,[],pAdjChow,[],'Non-Palatable')
%% Split data using median split of voracity data
% Find difference for nine animals across conditions (take advantage of
% NaNs)
% Dep - Chow; thus positive = faster chewing during Dep (palatable)
vorDiff = voracity(:,2)-voracity(:,4);
% Get non-NaN indices
inds = logicFind(1,~isnan(vorDiff),'==');
% Get median split
medSplit = median(voracity(inds,2));
bot = logicFind(medSplit,voracity(inds,2),'<=');
top = logicFind(medSplit,voracity(inds,2),'>');
% Do t-tests with only those nine animals
depDiff2 = depDiff(inds,:);
[h,p,pAdj] = bulkT(depDiff2,0,0,'fdr');
[hChow,pChow,pAdjChow] = bulkT(chowDiff,0,0,'fdr');
% Plot similarly to stripPlot but with colored means for medSplit
pInds = logicFind(0.05,pAdj,'<=');
pIndsChow = logicFind(0.05,pAdjChow,'<=');
% Get means of top and bot
c = 1;
for iP = pInds
   botMean(c) = mean(depDiff2(bot,iP));
   topMean(c) = mean(depDiff2(top,iP));
   c = c + 1;
end
c = 1;
for iP = pIndsChow
   botMeanChow(c) = mean(chowDiff(bot,iP));
   topMeanChow(c) = mean(chowDiff(top,iP));
   c = c + 1;
end
%% Plot voracity with split line
figure
hold on
for vi = inds
    plot([1 1.5],[voracity(vi,2) voracity(vi,4)],'-ok','MarkerFaceColor','k','Color',[0.5 0.5 0.5])
end
plot([1 1.5],[medSplit medSplit],'--')
xlim([0.9 1.6])
set(gca,'XTick',[1,1.5],'XTickLabel',{'PF','Chow'})
title('Voracity between Groups')
ylabel('Voracity (gm/sec)')
%% Plot - full all Binge-Rest
% Palatable
figure
hold on
c = 1;
for iP = pInds
    plot(repmat(iP,1,numel(top)),depDiff2(top,iP),'.k')
    plot([iP-0.25 iP+0.25],[topMean(c) topMean(c)],'-r')
    plot(repmat(iP,1,numel(bot)),depDiff2(bot,iP),'.','Color',[0.5 0.5 0.5])
    plot([iP-0.25 iP+0.25],[botMean(c) botMean(c)],'-b')
    c = c + 1;
end
plot([0:size(depDiff2,2)],zeros(size(depDiff2,2)+1),'-k')
title('Palatable: 9 subset, full')
legH(1) = plot(NaN,'.k');
legH(2) = plot(NaN,'-r');
legH(3) = plot(NaN,'.','Color',[0.5 0.5 0.5]);
legH(4) = plot(NaN,'-b');
legend(legH,{'Top','Top Mean','Bottom','Bottom Mean'},'Location','southwest')
% Chow
figure
hold on
c = 1;
for iP = pIndsChow
    plot(repmat(iP,1,numel(top)),chowDiff(top,iP),'.k')
    plot([iP-0.25 iP+0.25],[topMeanChow(c) topMeanChow(c)],'-r')
    plot(repmat(iP,1,numel(bot)),chowDiff(bot,iP),'.','Color',[0.5 0.5 0.5])
    plot([iP-0.25 iP+0.25],[botMeanChow(c) botMeanChow(c)],'-b')
    c = c + 1;
end
plot([0:size(chowDiff,2)],zeros(size(chowDiff,2)+1),'-k')
title('Non-Palatable: 9 subset, full')
legH(1) = plot(NaN,'.k');
legH(2) = plot(NaN,'-r');
legH(3) = plot(NaN,'.','Color',[0.5 0.5 0.5]);
legH(4) = plot(NaN,'-b');
legend(legH,{'Top','Top Mean','Bottom','Bottom Mean'},'Location','southwest')
%% Plot - condensed all
% Palatable
figure
hold on
c = 1;
for iP = pInds
    plot(repmat(c,1,numel(top)),depDiff2(top,iP),'.k')
    plot([c-0.25 c+0.25],[topMean(c) topMean(c)],'-r')
    plot(repmat(c,1,numel(bot)),depDiff2(bot,iP),'.','Color',[0.5 0.5 0.5])
    plot([c-0.25 c+0.25],[botMean(c) botMean(c)],'-b')
    c = c + 1;
end
plot([0:c],zeros(c+1),'-k')
title('Palatable: subset 9, condensed')
legH(1) = plot(NaN,'.k');
legH(2) = plot(NaN,'-r');
legH(3) = plot(NaN,'.','Color',[0.5 0.5 0.5]);
legH(4) = plot(NaN,'-b');
legend(legH,{'Top','Top Mean','Bottom','Bottom Mean'},'Location','southwest')
% Chow
figure
hold on
c = 1;
for iP = pIndsChow
    plot(repmat(c,1,numel(top)),chowDiff(top,iP),'.k')
    plot([c-0.25 c+0.25],[topMeanChow(c) topMeanChow(c)],'-r')
    plot(repmat(c,1,numel(bot)),chowDiff(bot,iP),'.','Color',[0.5 0.5 0.5])
    plot([c-0.25 c+0.25],[botMeanChow(c) botMeanChow(c)],'-b')
    c = c + 1;
end
plot([0:c],zeros(c+1),'-k')
title('Non-Palatable: subset 9, condensed')
legH(1) = plot(NaN,'.k');
legH(2) = plot(NaN,'-r');
legH(3) = plot(NaN,'.','Color',[0.5 0.5 0.5]);
legH(4) = plot(NaN,'-b');
legend(legH,{'Top','Top Mean','Bottom','Bottom Mean'},'Location','southwest')
%% Get number of times  bot is closer to zero
% Palatable
absBotMean = abs(botMean);
absTopMean = abs(topMean);
% Get percentage of bot being smaller and top being smaller
botSmall = sum(absTopMean>absBotMean)/length(absBotMean);
topSmall = sum(absTopMean<absBotMean)/length(absBotMean);
% Plot using log y scale
figure
for iM = 1:length(botMean)
    if absTopMean(iM) > absBotMean(iM)
        semilogy(iM,absTopMean(iM),'sr');
        hold on
        semilogy(iM,absBotMean(iM),'sb','MarkerFaceColor','b');
    else
        semilogy(iM,absTopMean(iM),'sr','MarkerFaceColor','r');
        hold on
        semilogy(iM,absBotMean(iM),'sb');
    end
end
title('Palatable: Means Split')
legH(1) = plot(NaN,'rs','MarkerFaceColor','r');
legH(2) = plot(NaN,'bs','MarkerFaceColor','b');
legend(legH,{['Top: ',num2str(round(topSmall*100,1)),'%'],['Bottom: ',num2str(round(botSmall*100,1)),'%']},'Location','southeast')
% Palatable
absBotMeanChow = abs(botMeanChow);
absTopMeanChow = abs(topMeanChow);
% Get percentage of bot being smaller and top being smaller
botSmallChow = sum(absTopMeanChow>absBotMeanChow)/length(absBotMeanChow);
topSmallChow = sum(absTopMeanChow<absBotMeanChow)/length(absBotMeanChow);
% Plot using log y scale
figure
for iM = 1:length(botMeanChow)
    if absTopMeanChow(iM) > absBotMeanChow(iM)
        semilogy(iM,absTopMeanChow(iM),'sr');
        hold on
        semilogy(iM,absBotMeanChow(iM),'sb','MarkerFaceColor','b');
    else
        semilogy(iM,absTopMeanChow(iM),'sr','MarkerFaceColor','r');
        hold on
        semilogy(iM,absBotMeanChow(iM),'sb');
    end
end
title('Non-Palatable: Means Split')
legH(1) = plot(NaN,'rs','MarkerFaceColor','r');
legH(2) = plot(NaN,'bs','MarkerFaceColor','b');
legend(legH,{['Top: ',num2str(round(topSmallChow*100,1)),'%'],['Bottom: ',num2str(round(botSmallChow*100,1)),'%']},'Location','southeast')
%% Elastic Net
% Prep response variable - 0 = chow, 1 = palatable
% y = [zeros(12,1);ones(9,1)];
% %% Binge
% % Stack binge data and add voracity measure
% xBinge = [data{1,1},voracity(:,2);data{2,1},voracity([1,2,3,5,6,7,8,10,11],4)];
% cfg = lassoNetCfg([],'n','y','n',100,'1se');
% [allAlpha,bingeLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(xBinge,y,'binomial','class',0:0.01:1,4,1,cfg);
% bingeBeta = allBeta{1,1}.survBeta;
% bingeErr = bingeLambda{1,1}.allErr;
% % Rest
% xRest = [data{1,2},voracity(:,2);data{2,2},voracity([1,2,3,5,6,7,8,10,11],4)];
% cfg = lassoNetCfg([],'n','y','n',100,'1se');
% [allAlpha,restLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(xRest,y,'binomial','class',0:0.01:1,4,1,cfg);
% restBeta = allBeta{1,1}.survBeta;
% restErr = restLambda{1,1}.allErr;
%% Binge-Rest 
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\palatableData.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\4conditionBingeSize.mat')
y = [zeros(12,1);ones(9,1)];
xNorm = [(data{1,1}-data{1,2}),voracity(:,2);(data{2,1}-data{2,2}),voracity([1,2,3,5,6,7,8,10,11],4)];
%%
cfg = lassoNetCfg([],'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(xNorm,y,'binomial','class',1,5,1,cfg);
%% Randomize Binge-Rest: use same alpha as above
cfg = lassoNetCfg([],'y','y','n',100,'1se');
[randAlpha,randLambda,randBeta,randCvFitsArray,randAccArray,randHist] = lassoNet(xNorm,y,'binomial','class',1,5,100,cfg);
randBeta = randBeta{1,1}.survBeta;
randErr = randLambda{1,1}.allErr;
%% Plot real Binge-Rest vs. Randomized
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\normLasso.mat')
realAcc = 1-allLambda{1,1}.allErr;
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\normLassoRand.mat')

randAcc = [];
for ii = 1:100
   randAcc = [randAcc;1-randLambda{1,ii}.allErr];
end

%% Histogram
figure
histogram(realAcc.*100,15,'Normalization','probability')
hold on
histogram(randAcc.*100,'Normalization','probability')
ylim([0 0.2])
xlabel('Accuracy')
ylabel('Proportion of Models')
title('Average Accuracy in Palatability Prediction: Binge-Rest','FontName','Arial')
legend({'Real','Permuted'},'location','northwest')
%% ECDF
figure
ecdf(realAcc.*100)
hold on
ecdf(randAcc.*100)
xlabel('Accuracy')
ylabel('Cumulative Density')
title('Palatability Prediction: Binge-Rest')
legend({'Real','Permuted'},'location','northwest')

[h,p] = kstest2(realAcc,randAcc);
%% Predict baseline binge size from baseline (b-r)
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\4conditionData.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\4conditionBingeSize.mat')
x = data{1,1}-data{1,2};
y = bingeSizes(:,1);
cfg = lassoNetCfg([],'y','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x,y,'gaussian','mae',1,4,100,cfg);
%% Plot baseline binge size prediction
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\baseBingeSizeLasso.mat')
realErr = allLambda{1,1}.allErr;
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\baseBingeSizeLassoRand.mat')

mReal = mean(realErr);
sReal = std(realErr);
mRand = mean(randErr);
sRand = std(randErr);
%%
figure
histogram(realErr,'BinWidth',0.1,'Normalization','probability')
hold on
histogram(randErr,'Normalization','probability')
title('Predicting Baseline Binge Size')
xlabel('Mean Absolute Error: Grams')
ylabel('Proportion of Models')
legend({['Real: \mu = ',num2str(round(mReal),2),'\pm',num2str(round(sReal,2)),' gm'],['Permuted: \mu = ',num2str(round(mRand,2)),'\pm',num2str(round(sRand,2)),' gm']})
%% Try to predict change in binge size from baseline using both deps (24 and 48)
% Percent change in binge size from base
percBinge24 = (bingeSizes(:,2)-bingeSizes(:,1))./bingeSizes(:,1);
percBinge48 = (bingeSizes(:,3)-bingeSizes(:,1))./bingeSizes(:,1);
percBinge = [percBinge24;percBinge48];
% Remove NaNs
percBinge = percBinge(~isnan(percBinge));
% Subtract baseline data from deps (24 and 48); then stack
dep24 = (data{2,1}-data{2,2}) - (data{1,1}-data{1,2});
% Get indices using bingeSize NaNs
dInds = logicFind(1,~isnan(bingeSizes(:,3)),'==');
dep48 = (data{3,1}-data{3,2}) - (data{1,1}(dInds,:)-data{1,2}(dInds,:));
bingeRest = [dep24;dep48];
% Predict change in binge size
cfg = lassoNetCfg([],'y','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(bingeRest,percBinge,'gaussian','mae',1,5,100,cfg);
%% Plot binge size change prediction
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\percBingeSizeLassoRand.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\percBingeSizeLasso.mat')
figure
histogram(allLambda{1,1}.allErr.*100,'BinWidth',2,'Normalization','probability')
hold on
histogram(randErr.*100,'Normalization','probability')
title('Predicting % Change in Binge Size from Baseline to Deprivation')
xlabel('% \Delta in Binge Size')
ylabel('Proportion of Models')
legend({['Real: \mu = ',num2str(round(mean(allLambda{1,1}.allErr.*100),2)),'\pm',num2str(round(std(allLambda{1,1}.allErr.*100),2)),'%'],['Permuted: \mu = ',num2str(round(mean(randErr.*100),2)),'\pm',num2str(round(std(randErr.*100),2)),'%']})
%% Try to predict amount eaten (calories)
normX = [data{1,1}-data{1,2};data{2,1}-data{2,2};data{3,1}-data{3,2};data{4,1}-data{4,2}];
dInds = logicFind(1,~isnan(bingeCal(:,3)),'=='); 
cInds = logicFind(1,~isnan(bingeCal(:,4)),'=='); 
calEaten = [bingeCal(:,1);bingeCal(:,2);bingeCal(dInds,3);bingeCal(cInds,4)];
%%
cfg = lassoNetCfg([],'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(normX,calEaten,'gaussian','mse',0:0.01:1,5,1,cfg);
%% Try to predict state of animal 
%1 = base, 0 = 48dep
normX = [(data{1,1}-data{1,2})./data{1,2};(data{3,1}-data{3,2})./data{3,2}];%;data{3,1}-data{3,2};data{4,1}-data{4,2}];
% normX = [(data{1,1}-data{1,2});(data{2,1}-data{2,2})];
state = [repmat(0,12,1);repmat(1,9,1)];%;repmat(1,9,1);repmat(4,9,1)];

cfg = lassoNetCfg([],'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(normX,state,'binomial','class',0:0.01:1,5,1,cfg);
%% Plot percent change in feature against absolute difference in voracity
% dep48 - base; thus positive = faster chewing during Dep48
vorDiff = [voracity(:,3)-voracity(:,1);voracity(:,2)-voracity(:,1)];
% Get non-NaN indices
inds24 = logicFind(1,~isnan(voracity(:,2)),'==');
inds48 = logicFind(1,~isnan(voracity(:,3)),'==');

dep48Feats = data{3,1}-data{3,2};
dep24Feats = data{2,1}-data{2,2};
baseFeats = data{1,1}-data{1,2};

feats = [(dep48Feats-baseFeats(inds48,:))./baseFeats(inds48,:);(dep24Feats-baseFeats(inds24,:))./baseFeats(inds24,:)];
%%
for fi = 1:size(feats,2)
    [thisR,thisP] = corrcoef(vorDiff(~isnan(vorDiff)),feats(:,fi));
    r(fi) = thisR(1,2)^2;
    p(fi) = thisP(1,2);
end
% Get indices of significant p-values to exclude those features
pInds = logicFind(0.05,p,'<=');
%% Baseline Binge vs. Rest T-tests
[h,p,pAdj] = bulkT(data{1,1}-data{1,2},0,1,'fdr');
% Effectively remove data with high correlations by setting their pAdj to 1
pAdj(pInds) = 1;
%% Plot power features t-test
stripPlot((data{1,1}(:,1:24)-data{1,2}(:,1:24)).*100,[],pAdj(1:24),names(1:24),{'Baseline Binge - Rest: Power Features'})
ylabel('Difference in Binge from Rest (%)')
%% Plot coherence features t-test
stripPlot((data{1,1}(:,25:60)-data{1,2}(:,25:60)).*100,[],pAdj(25:60),names(25:60),{'Baseline Binge - Rest: Coherence Features'})
ylabel('Difference in Binge from Rest (%)')
%% Plot power-corr features t-test
stripPlot(data{1,1}(:,61:end)-data{1,2}(:,61:end),[],pAdj(61:end),names(61:end),{'Baseline Binge - Rest: Power Correlation Features'})
ylabel('Difference in Binge from Rest (R)')
%% Plot 'best' and 'worst' correlation (lowest and highest p value)
minInd = logicFind(min(p),p,'==');
maxInd = logicFind(max(p),p,'==');
%%
figure
plot(vorDiff(~isnan(vorDiff)),feats(:,minInd),'.k','MarkerSize',10)
lsline
title('Best Correlation: Core Left to Core Right \theta Coherence')
xlabel('\Delta Voracity (gm/sec)')
ylabel('\Delta in % Normalized Coherence')
text(0,-2,['R^2 = ',num2str(round(r(minInd),2))])
text(0,-2.25,['p = ',num2str(p(minInd),'%.1d')])
%% Hand picked index 18 for fewer 'outliers'
figure
plot(vorDiff(~isnan(vorDiff)),feats(:,18),'.k','MarkerSize',10)
lsline
title('Best Correlation: Core Left High \gamma Power')
xlabel('\Delta Voracity (gm/sec)')
ylabel('\Delta in % Normalized Power')
text(0,-5,['R^2 = ',num2str(r(18),'%.1d')])
text(0,-5.5,['p = ',num2str(round(p(18),2))])
%%
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\H10BaseOct15_binge_vs_rest.mat')
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\H10FoodDep24Sep25_binge_vs_rest.mat')
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\H10FoodDep48Oct14_binge_vs_rest.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\processed\H13BaseNov12_binge_vs_rest.mat')
for j = 1:size(hist.bands,1)
    bandInd(j,1) = find(psdTrls.F>=hist.bands{j,2}(1),1);
    bandInd(j,2) = find(psdTrls.F<=hist.bands{j,2}(2),1,'last');
end
% Grab all power and normalize
for ii = 1:size(psdTrls.event1.Pow,2)
    for c = 1:4
        bingeTot(ii,c) = trapz(psdTrls.event1.Pow{1,ii}(c,1:45));
    end
    bingeRel(:,:,ii) = psdTrls.event1.Pow{2,ii}./repmat(bingeTot(ii,:),6,1);
    bingePow(ii,:) = reshape(bingeRel(:,:,ii),1,24);
end
for ii = 1:size(psdTrls.event2.Pow,2)
    for c = 1:4
        restTot(ii,c) = trapz(psdTrls.event2.Pow{1,ii}(c,1:45));
    end
    restRel(:,:,ii) = psdTrls.event2.Pow{2,ii}./repmat(restTot(ii,:),6,1);
    restPow(ii,:) = reshape(restRel(:,:,ii),1,24);
end
% Grab all coherence and normalize
% Band x cmb x trial
for bi = 1:size(bandInd,1)
    bingeAvgCoh(bi,:,:) = mean(coh{1,1}.Cxy(:,bandInd(bi,1):bandInd(bi,2),:),2);
    restAvgCoh(bi,:,:) = mean(coh{1,2}.Cxy(:,bandInd(bi,1):bandInd(bi,2),:),2);
end
% cmb x trial
allBingeCoh = squeeze(mean(coh{1,1}.Cxy(:,1:45,:),2));
allRestCoh = squeeze(mean(coh{1,2}.Cxy(:,1:45,:),2));
for ii = 1:size(allBingeCoh,2)
    bingeRelCoh(:,:,ii) = bingeAvgCoh(:,:,ii)./repmat(allBingeCoh(:,ii)',6,1);
    bingeCoh(ii,:) = reshape(bingeRelCoh(:,:,ii),1,36);
end
for ii = 1:size(allRestCoh,2)
    restRelCoh(:,:,ii) = restAvgCoh(:,:,ii)./repmat(allRestCoh(:,ii)',6,1);
    restCoh(ii,:) = reshape(restRelCoh(:,:,ii),1,36);
end
% Combine power and coherence
allData = [bingePow,bingeCoh;restPow,restCoh];
resp = [ones(size(bingePow,1),1);zeros(size(restPow,1),1)];
%% Predict binge vs. rest - within animal, same condition
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H10Base1Dat.mat')
cfg = lassoNetCfg(0.2,'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(base1Dat,base1Resp,'binomial','auc',1,5,1,cfg);
% Construct AUC
[baseX,baseY,~,baseA] = perfcurve(base1Resp(hist.testInd),accArray{1,1}.pred,1);
%% Predict binge vs. rest - within animal, base to dep24
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H10Base1Dat.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H10Dep24Dat.mat')
% Uses all data from base to train, and all from dep to test
cfg = lassoNetCfg({dep24Dat,dep24Resp},'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(base1Dat,base1Resp,'binomial','auc',1,5,1,cfg);
[baseDepX,baseDepY,~,baseDepA] = perfcurve(dep24Resp,accArray{1,1}.pred,1);
%% Predict binge vs. rest - within animal, base to dep48
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H10Base1Dat.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H10FoodDep48Dat.mat')
% Uses all data from base to train, and all from dep to test
cfg = lassoNetCfg({dep48Dat,dep48Resp},'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(base1Dat,base1Resp,'binomial','auc',1,5,1,cfg);
[baseDep48X,baseDep48Y,~,baseDep48A] = perfcurve(dep48Resp,accArray{1,1}.pred,1);
%% Predict binge vs. rest - across animals (H10-H13) base to base
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H10Base1Dat.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H13BaseNov12Dat.mat')
% Uses all data from base to train, and all from dep to test
cfg = lassoNetCfg({h13BaseDat,h13BaseResp},'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(base1Dat,base1Resp,'binomial','auc',1,5,1,cfg);
[h10Baseh13BaseX,h10Baseh13BaseY,~,h10Baseh13BaseA] = perfcurve(h13BaseResp,accArray{1,1}.pred,1);
%%
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H10BaseOct15_Dep24Sep25Lasso.mat','baseDepA','baseDepX','baseDepY')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H10BaseOct15_Dep48Oct14Lasso.mat','baseDep48A','baseDep48X','baseDep48Y')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\H10BaseOct15Lasso.mat','baseA','baseX','baseY')
load
figure
plot(baseX,baseY,'LineWidth',3)
hold on
plot(baseDepX,baseDepY,'LineWidth',3)
plot(baseDep48X,baseDep48Y,'LineWidth',3)
plot(h10Baseh13BaseX,h10Baseh13BaseY,'LineWidth',3)
legend({'Base','Dep 24','Dep 48'})
%% Predict binge vs. rest - within animal, 24dep


