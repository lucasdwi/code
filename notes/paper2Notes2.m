% % Load ephys data for all trials of binge and rest
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\bingeRestAllTrial.mat')
% for ii = 1:size(data,2)
%     for k = 1:size(data{1,ii},2)
%         trialDat{k,ii} = cat(1,data{1,ii}{:,k});
%     end
% end
% Load ephys average data of binge and rest
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeRestAllAvg.mat')
for ii = 1:size(data,2)
    for k = 1:size(data{1,ii},2)
        avgDat{k,ii} = cat(1,data{1,ii}{:,k});
    end
end
% Load voracity data
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\4conditionBingeSize.mat')
% Correlate voracity with all features
vorDiff = [voracity(:,3)-voracity(:,1);voracity(:,2)-voracity(:,1)];
% Get non-NaN indices
inds24 = logicFind(1,~isnan(voracity(:,2)),'==');
inds48 = logicFind(1,~isnan(voracity(:,3)),'==');
% Do same subtractions for features
dep48Feats = avgDat{1,3}-avgDat{2,3};
dep24Feats = avgDat{1,2}-avgDat{2,2};
baseFeats = avgDat{1,1}-avgDat{2,1};

feats = [(dep48Feats-baseFeats(inds48,:))./baseFeats(inds48,:);(dep24Feats-baseFeats(inds24,:))./baseFeats(inds24,:)];
for fi = 1:size(feats,2)
    [thisR,thisP] = corrcoef(vorDiff(~isnan(vorDiff)),feats(:,fi));
    r(fi) = thisR(1,2)^2;
    p(fi) = thisP(1,2);
end
% Get indices of significant p-values to exclude those features
pInds = logicFind(0.05,p,'<=');
% Save
% save('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\featCorr.mat','avgDat','pInds','feats','vorDiff')
%% Predicting baseline binge size using baseline data
x = avgDat{1,1}-avgDat{2,1};
y = bingeSizes(:,1);
% Remove above indices which may be contaminated by chewing noise
x = x(:,~ismember(1:size(avgDat{1,1},2),pInds));
% Real data
cfg = lassoNetCfg([],[],'n','y','n',100,'1se');
[~,allLambda,allBeta,cvFitsArray,~,hist] = lassoNet(x,y,'gaussian','mae',1,4,1,cfg);
real.lambda = allLambda; real.beta = allBeta; real.fits = cvFitsArray; real.hist = hist; 
real.err = real.lambda{1,1}.allErr;
% Permuted data
cfg = lassoNetCfg([],[],'y','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x,y,'gaussian','mae',1,4,100,cfg);
rand.lambda = allLambda; rand.beta = allBeta; rand.fits = cvFitsArray; rand.hist = hist;
rand.err = [];
for ii = 1:100
   rand.err = [rand.err;rand.lambda{ii}.allErr];
end
% Get Mann-Whitney U statistic
[p,h,stats] = ranksum(real.err,rand.err);
% Get r-family effect size
r = abs(stats.zval/sqrt(numel(real.err)+numel(rand.err)));
% Get Cohen's d from r
d = (2*r)/sqrt(1-r.^2);
% save('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\predictBaseBinge.mat','real','rand','d')
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\predictBaseBinge.mat')

figure
histogram(real.err,'Normalization','probability','BinWidth',.1)
hold on
histogram(rand.err,'Normalization','probability','BinWidth',.1)
title('Predicting Baseline Binge Size: Base Binge - Base Rest')
xlabel('Mean Absolute Error (gm)')
ylabel('Model Frequency')
legend({['Real: \mu = ',num2str(round(mean(real.err)),2),'\pm',num2str(round(std(real.err),2)),' gm'],['Permuted: \mu = ',num2str(round(mean(rand.err),2)),'\pm',num2str(round(std(rand.err),2)),' gm']})
%% Subset predicting baseline binge size
for ii = 1:14
   load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\baseBingeSize\baseBingeSize',num2str(ii),'.mat']) 
   realAll{ii} = real{ii};
   permAll{ii} = perm{ii};
end
%% Predicting percent change in binge size from baseline to deprived states 
% using difference 
% Percent change in binge size from base
percBinge24 = (bingeSizes(:,2)-bingeSizes(:,1))./bingeSizes(:,1);
percBinge48 = (bingeSizes(:,3)-bingeSizes(:,1))./bingeSizes(:,1);
percBinge = [percBinge24;percBinge48];
% Remove NaNs
percBinge = percBinge(~isnan(percBinge));
% Normalize dep24 by baseline: (binge-rest)-(binge-rest)
dep24 = (avgDat{1,2}-avgDat{2,2}) - (avgDat{1,1}-avgDat{2,1});
% Get indices using bingeSize NaNs
dInds = logicFind(1,~isnan(bingeSizes(:,3)),'==');
% Normalize dep48 by baseline using matched indices: (binge-rest)-(binge-rest)
dep48 = (avgDat{1,3}-avgDat{2,3}) - (avgDat{1,1}(dInds,:)-avgDat{2,1}(dInds,:));
x = [dep24;dep48];
x = x(:,~ismember(1:size(avgDat{1,1},2),pInds));
% Real data
cfg = lassoNetCfg([],[],'n','y','n',100,'1se');
[~,allLambda,allBeta,cvFitsArray,~,hist] = lassoNet(x,percBinge,'gaussian','mae',1,5,1,cfg);
real.lambda = allLambda; real.beta = allBeta; real.fits = cvFitsArray; real.hist = hist; 
real.err = real.lambda{1,1}.allErr;
% Permuted data
cfg = lassoNetCfg([],[],'y','y','n',100,'1se');
[~,allLambda,allBeta,cvFitsArray,~,hist] = lassoNet(x,percBinge,'gaussian','mae',1,5,100,cfg);
rand.lambda = allLambda; rand.beta = allBeta; rand.fits = cvFitsArray; rand.hist = hist;
rand.err = [];
for ii = 1:100
   rand.err = [rand.err;rand.lambda{ii}.allErr];
end
% Get Mann-Whitney U statistic
[p,h,stats] = ranksum(real.err,rand.err);
% Get r-family effect size
r = abs(stats.zval/sqrt(numel(real.err)+numel(rand.err)));
% Get Cohen's d from r
d = (2*r)/sqrt(1-r.^2);
%save('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\predictBingeChange.mat','real','rand','d')
%load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\predictBingeChange.mat')
figure
histogram(real.err.*100,'Normalization','probability','BinWidth',5)
hold on
histogram(rand.err.*100,'Normalization','probability','BinWidth',5)
title('Predicting %\Delta in Binge Size: Base(Binge-Rest)-Dep(Binge-Rest)')
xlabel('Mean Absolute Error (%)')
ylabel('Model Frequency')
legend({['Real: \mu = ',num2str(round(mean(real.err)*100),2),'\pm',num2str(round(std(real.err)*100,2)),' %'],['Permuted: \mu = ',num2str(round(mean(rand.err)*100,2)),'\pm',num2str(round(std(rand.err)*100,2)),' %']})
%% Subset predicting change in binge from baseline
for ii = 1:14
   load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeSizeChange\bingeSizeChange',num2str(ii),'.mat']) 
   realAll{ii} = real{ii};
   if ii == 1
       permAll = perm{ii};
   end
end
%% Predicting palatability of food (chow vs. dep24)
x = [avgDat{1,4}-avgDat{2,4};avgDat{1,2}-avgDat{2,2}];
x = x(:,~ismember(1:size(avgDat{1,1},2),pInds));
y = [ones(size(avgDat{1,4},1),1);zeros(size(avgDat{1,2},1),1)];
% Real data
cfg = lassoNetCfg([],[],'n','y','n',100,'1se');
[~,allLambda,allBeta,cvFitsArray,~,hist] = lassoNet(x,y,'binomial','class',1,5,1,cfg);
real.lambda = allLambda; real.beta = allBeta; real.fits = cvFitsArray; real.hist = hist; 
real.err = real.lambda{1,1}.allErr;
% Permuted data
cfg = lassoNetCfg([],[],'y','y','n',100,'1se');
[~,allLambda,allBeta,cvFitsArray,~,hist] = lassoNet(x,y,'binomial','class',1,5,100,cfg);
rand.lambda = allLambda; rand.beta = allBeta; rand.fits = cvFitsArray; rand.hist = hist;
rand.err = [];
for ii = 1:100
   rand.err = [rand.err;rand.lambda{ii}.allErr];
end
% Get Mann-Whitney U statistic
[p,h,stats] = ranksum(real.err,rand.err);
% Get r-family effect size
r = abs(stats.zval/sqrt(numel(real.err)+numel(rand.err)));
% Get Cohen's d from r
d = (2*r)/sqrt(1-r.^2);
% save('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\predictPalatability.mat','rand','real','d')
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\predictPalatability.mat')
figure
histogram((1-real.err).*100,'Normalization','probability','BinWidth',2)
hold on
histogram((1-rand.err).*100,'Normalization','probability','BinWidth',2)
title('Predicting Palatability: Chow(Binge-Rest) vs. Dep(Binge-Rest)')
xlabel('Accuracy (%)')
ylabel('Model Frequency')
legend({['Real: \mu = ',num2str(round(mean(1-real.err)*100),2),'\pm',num2str(round(std(1-real.err)*100,2)),' %'],['Permuted: \mu = ',num2str(round(mean(1-rand.err)*100,2)),'\pm',num2str(round(std(1-rand.err)*100,2)),' %']},'Location','northwest')
%% Predicting palatability of food (chow vs. dep24)
% Set up real and rand arrays
real = cell(1,size(bingeRestInd,2));
rand = cell(size(real));
for ii = 1:size(bingeRestInd,2)
    % Set up data set
    x = [avgDat{1,4}-avgDat{2,4};avgDat{1,2}-avgDat{2,2}];
    % Remove pInds, if they exist in this subset
    inds = bingeRestInd{1,ii}(~ismember(bingeRestInd{1,ii},pInds));
    % Subset
    x = x(:,inds);
    y = [ones(size(avgDat{1,4},1),1);zeros(size(avgDat{1,2},1),1)];
    % Real data
    cfg = lassoNetCfg([],[],'n','y','n',100,'1se');
    [~,allLambda,allBeta,cvFitsArray,~,hist] = lassoNet(x,y,'binomial','class',1,5,1,cfg);
    real{ii}.lambda = allLambda; real{ii}.beta = allBeta; real{ii}.fits = cvFitsArray; real{ii}.hist = hist;
    real{ii}.err = real{ii}.lambda{1,1}.allErr;
    % Permuted data
    cfg = lassoNetCfg([],[],'y','y','n',100,'1se');
    [~,allLambda,allBeta,cvFitsArray,~,hist] = lassoNet(x,y,'binomial','class',1,5,100,cfg);
    rand{ii}.lambda = allLambda; rand{ii}.beta = allBeta; rand{ii}.fits = cvFitsArray; rand{ii}.hist = hist;
    rand{ii}.err = [];
    for k = 1:100
        rand{ii}.err = [rand{ii}.err;rand{ii}.lambda{k}.allErr];
    end
end
%% Convert trial data into sets
for ii = 1:size(data,2)
    for k = 1:size(data{ii},1)
        thisX = [data{1,ii}{k,1};data{1,ii}{k,2}];
        thisY = [ones(size(data{1,ii}{k,1},1),1);zeros(size(data{1,ii}{k,2},1),1)];
        trlDat{ii,k} = [thisY,thisX]; 
    end
end
save('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeRestAllTrial.mat','data','trlDat')
%%
for ii = 1%:size(trlDat,2)
    % First build model and test on self
    cfg = lassoNetCfg(0.20,'balanced','n','y','n',100,'1se');
    [~,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(trlDat{1,ii}(:,2:end),trlDat{1,ii}(:,1),'binomial','auc',1,5,1,cfg);    
    % Find internal accuracy
    [allX{ii},allY{ii},~,allA(ii,ii)] = perfcurve(trlDat{1,ii}(hist.testInd,1),accArray{1,1}.pred,1);
     % Find accuracy when applied to other animals
    for oi = 1:size(trlDat,2)
        if oi ~= ii
            thisDat = trlDat{1,oi}(:,2:end);
            [predY] = cvglmnetPredict(cvFitsArray{1,1}{allLambda{1,1}.bestLambdaInds},thisDat,'lambda_1se','response');
            [allX{ii,oi},allY{ii,oi},~,allA(ii,oi)] = perfcurve(trlDat{1,oi}(:,1),predY,1);
        end
    end
end
%% Binge vs. not-binge 
% Convert trial data into sets
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeNotBingeTrial.mat')
for ii = 1:size(data,2)
    for k = 1:size(data{ii},1)
        thisX = [data{1,ii}{k,1};data{1,ii}{k,2}];
        thisY = [ones(size(data{1,ii}{k,1},1),1);zeros(size(data{1,ii}{k,2},1),1)];
        trlDat{ii,k} = [thisY,thisX];
        sampDat{ii,k} = [samp{1,ii}{k,1};samp{1,ii}{k,2}];
    end
end
%% Use each animal's baseline to predict all other animals
for ii = 1:size(trlDat,2)
    cfg = lassoNetCfg(0.20,'balanced','n','y','n',100,'1se');
    [~,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(trlDat{1,ii}(:,2:end),trlDat{1,ii}(:,1),'binomial','auc',1,5,1,cfg);
    % Find internal accuracy
    [allX{ii},allY{ii},~,allA(ii,ii)] = perfcurve(trlDat{1,ii}(hist.testInd,1),accArray{1,1}.pred,1);
    % Apply to other files
    for oi = 1:size(trlDat,2)
        if oi ~= ii
            thisDat = trlDat{1,oi}(:,2:end);
            [predY{oi}] = cvglmnetPredict(cvFitsArray{1,1}{allLambda{1,1}.bestLambdaInds},thisDat,'lambda_1se','response');
            [allX{ii,oi},allY{ii,oi},~,allA(ii,oi)] = perfcurve(trlDat{1,oi}(:,1),predY{oi},1);
        end
    end
end
% save('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\predictBingeNot.mat')
%% Binge vs. Not Subsets
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeNot\')
for ii = 1:13
   load([num2str(ii),'.mat'],'rocX','rocY','auc')
   x{ii} = rocX; y{ii} = rocY;
   a(ii) = auc;
end
save('subsetBingeNot.mat','x','y','a')
%% Subsetting indices for trial data (i.e. 60 dimensions; 24 from power and
% 36 from coherence)
% Power indices (channel,band)
powInd = reshape(1:24,6,4)';
% Coherence indices (cmb,band)
cohInd = reshape(25:60,6,6)';
% Delta indices
dInd = [powInd(:,1);cohInd(:,1)]';
% Theta indices
tInd = [powInd(:,2);cohInd(:,2)]';
% Alpha indices
aInd = [powInd(:,3);cohInd(:,3)]';
% Beta indices
bInd = [powInd(:,4);cohInd(:,4)]';
% Low gamma indices
lgInd = [powInd(:,5);cohInd(:,5)]';
% High gamma indices
hgInd = [powInd(:,6);cohInd(:,6)]';
% Left indices
lInd = [powInd([1,3],:);cohInd(2,:)];
% Right indices
rInd = [powInd([2,4],:);cohInd(5,:)];
% Shell indices
sInd = [powInd([1,2],:);cohInd(1,:)];
% Core indices
cInd = [powInd([3,4],:);cohInd(6,:)];
% Group together all indices of groups for binge vs. not-binge
bingeNotInd = {reshape(powInd',1,24),reshape(cohInd',1,36),dInd,tInd,aInd,bInd,lgInd,hgInd,lInd,rInd,sInd,cInd};
%% Subsetting indices for average data 
% Power indices (channel,band)
powInd = reshape(1:24,6,4)';
% Coherence indices (cmb,band)
cohInd = reshape(25:60,6,6)';
% Power correlation indices
corrInd = reshape(61:96,6,6)';
% Delta indices
dInd = [powInd(:,1);cohInd(:,1);corrInd(:,1)]';
% Theta indices
tInd = [powInd(:,2);cohInd(:,2);corrInd(:,2)]';
% Alpha indices
aInd = [powInd(:,3);cohInd(:,3);corrInd(:,3)]';
% Beta indices
bInd = [powInd(:,4);cohInd(:,4);corrInd(:,4)]';
% Low gamma indices
lgInd = [powInd(:,5);cohInd(:,5);corrInd(:,5)]';
% High gamma indices
hgInd = [powInd(:,6);cohInd(:,6);corrInd(:,6)]';
% Left indices
lInd = [powInd([1,3],:);cohInd(2,:);corrInd(2,:)];
% Right indices
rInd = [powInd([2,4],:);cohInd(5,:);corrInd(5,:)];
% Shell indices
sInd = [powInd([1,2],:);cohInd(1,:);corrInd(1,:)];
% Core indices
cInd = [powInd([3,4],:);cohInd(6,:);corrInd(6,:)];
% Group together all indices of groups for binge vs. rest
bingeRestInd = {1:96,reshape(powInd',1,24),reshape(cohInd',1,36),reshape(corrInd',1,36),dInd,tInd,aInd,bInd,lgInd,hgInd,lInd,rInd,sInd,cInd};
%% Compare subsetting
load('subsetBaseBingeSize.mat')
for ii = 1:14
    mErr(1,ii) = mean(real{ii}.err);
    sErr(1,ii) = std(real{ii}.err);
end
pMErr(1) = mean(perm.err);
pSErr(1) = std(perm.err);
load('subsetBingeSizeChange.mat')
for ii = 1:14
    mErr(2,ii) = mean(real{ii}.err);
    sErr(2,ii) = std(real{ii}.err);
end
pMErr(2) = mean(perm.err);
pSErr(2) = std(perm.err);
load('subsetPalatability.mat')
for ii = 1:14
    mErr(3,ii) = mean(real{ii}.err);
    sErr(3,ii) = std(real{ii}.err);
end
pMErr(3) = mean(perm.err);
pSErr(3) = std(perm.err);
% Rank all but first column
[~,sorted] = sort(mErr(:,2:end),2);
for iR = 1:size(sorted,1)
   for iC = 1:size(sorted,2)
      ranks(iR,iC) = logicFind(iC,sorted(iR,:),'=='); 
   end
end
[~,mRanks] = sort(mean(ranks,1),2);
for iC = 1:size(mRanks,2)
   ranks(4,iC) = logicFind(iC,mRanks,'==');
end
load('subsetBingeNot.mat')
% Add NaN for power corr
aNaN = [a(2:3),NaN,a(4:end)];
[~,sorted] = sort(aNaN,2,'descend');
for iC = 1:13
    if iC == 3
        ranks(5,iC) = NaN;
    else
        ranks(5,iC) = logicFind(iC,sorted(1,:),'==')-1;
    end
end
% Set all .50s to the same value
inds = logicFind(0.5,round(aNaN,2),'==');
ranks(5,inds) = ranks(5,inds(1));
figure
imagesc(ranks)
set(gca,'XTick',1:13,'XTickLabel',{'pow','coh','corr','\Delta','\theta','\alpha','\beta','l\gamma','h\gamma','left','right','shell','core'},'YTick',1:4,'YTickLabel',{'Base Binge Size','\Delta Binge Size','Palatability','Ranked Mean'})
title('Subset Ranks (from bingeRestAllData avgData)')
save('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\subsetRanks.mat','mErr','pMErr','pSErr','sErr','ranks')
%% Collate bingeNot data
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeNot\')
for ii = 1:13
    load([num2str(ii),'.mat'])
    testXAll(:,:,ii) = testX;
    testYAll(:,ii) = testY;
    trainXAll(:,:,ii) = trainX;
    trainYAll(:,ii) = trainY;
    alpha{ii} = allAlpha{1,1};
    lambda{ii} = allLambda{1,1};
    beta{ii} = allBeta{1,1};
    acc{ii} = accArray{1,1};
    histAll{ii} = hist;
    if auc == 0.5
        allRocX(ii,:) = 0:1/2515:1;
        allRocY(ii,:) = 0:1/2515:1;
    else
        allRocX(ii,:) = rocX;
        allRocY(ii,:) = rocY;
    end
    allAuc(ii) = auc;
end
%%
figure
for ii = 1:20
   hold on
   plot(allRocX(ii,:),allRocY(ii,:))
end
%%
figure
for ii = 1:size(real,2)
   % Create normal distribution from statistics of permuted data
   normPerm = normrnd(mean(perm{ii}.err),std(perm{ii}.err),10000,1); 
   [~,pCrit(ii)] = kstest2(normPerm,perm{ii}.err);
   mRErr(ii) = mean(real{ii}.err);
   mPErr(ii) = mean(perm{ii}.err);
   [~,p(1,ii)] = kstest2(real{ii}.err,perm{ii}.err);
   [~,p(2,ii)] = ttest2(real{ii}.err,perm{ii}.err);
   hold on
   ecdf(real{ii}.err);
   if ii == 1
      h = get(gca,'children');
      set(h,'LineWidth',2) 
   end
end
%% Other conds
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\otherConds.mat')
% dep24 = 1,4,7,10,13,16,22,25,28,29,30
% dep48 = 2,5,8,11,14,17,20,23,26
% chow = 3,6,9,12,15,18,21,24,27
%%
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\otherCond\')
for ii = 1:30
   load(['otherCond',num2str(ii),'.mat'])
   allX{ii} = x;
   allY{ii} = y;
   allA{ii} = a;
end
%%
dep24 = [1,4,7,10,13,16,19,22,25,28,29,30];
dep48 = [2,5,8,11,14,17,20,23,26];
chow = [3,6,9,12,15,18,21,24,27];
figure
for ii = 1:30
    hold on
    h = plot(allX{ii},allY{ii});
    if ismember(ii,dep24)
       set(h,'Color','g')
    elseif ismember(ii,dep48)
        set(h,'Color','b')
    elseif ismember(ii,chow)
        set(h,'Color','r')
    end
end
plot([0,1],[0,1],'--k','LineWidth',2)
%% Binge Not One to All Other Animals - Baseline
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeNotBingeTrial.mat')
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeNotOneToAll\')
for ii = 1:12
    load([num2str(ii),'.mat'])
    [~,~,~,auc(1,ii)] = perfcurve(trlDat{1,ii}(hist.testInd,1),accArray{1,1}.pred,1);
    auc(2,ii) = a;
end
%% Binge Not one animal baseline to all other animals and conditions
for ii = 1:12
   load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\oneToRest\',num2str(ii),'.mat'])
   x{ii,:} = allX{ii,:}; y{ii,:} = allY{ii,:};
   a(ii,:) = allA(ii,:);
end
%%
for ii = 1:20
    load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\test\',num2str(ii),'.mat'])
    x{ii,:} = allX{ii,:}; y{ii,:} = allY{ii,:};
    a(ii,:) = allA(ii,:);
    lam(ii,:) = allLambda{1,1}.allErrAvg;
end
