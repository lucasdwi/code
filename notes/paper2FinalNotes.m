%% Grab preBinge data 
[data,samp,preFiles] = collateData(['C:\Users\Pythia\Documents\GreenLab'...
    '\data\paper2\preBingeCombined\'],{'base'},{'pow','coh'},'trl','rel');
% Only keep preBinge epochs corresponding to non-overlapping windows (i.e.
% [-5 0],[-10 -5],[-15 -10],[-20 -15],[-25 -20],[-30 -25],[-35 -30], and
% [-40 -35]) corresponds to cells 1,6,11,16,21,26,31, and 36 out of the
% original 61 overlapping windows.
preData = data{1,1}(:,1:5:36);
preSamp = samp{1,1}(:,1:5:36);
% Concatenate all 'preData'
preCat = cat(1,preData{:,1});
% Generate 20 80/20 splits of pre data
preTestX = cell(1,20);
thisTrainX = cell(1,20);
newPre = cell(1,20);
preTrainX = cell(1,20);
preTrainY = cell(1,20);
preTestY = cell(1,20);
% load notFeeding data to use as majority case in ADASYN
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'],'all')
for ii = 1:20
    % Generate random indices corresponding to 80/20 split
    rI = randperm(size(preCat,1),round(size(preCat,1)*.20));
    % Split data
    preTestX{ii} = preCat(rI,:);
    thisTrainX{ii} = preCat(~ismember(1:size(preCat,1),rI),:);
    % Extract not feeding data
    notFeed = all.trainX{ii}(all.trainY{ii}==0,:);
    % Use 'notFeed' as majority case in ADASYN to impute preFeeding up to
    % 175 per animal, or 300 total per animal; uses 2200 to ensure enough
    % samples (need 2100)
    [newPre{ii},~] = ADASYN([thisTrainX{ii};notFeed(randperm(...
        size(notFeed,1),2200),:)],[ones(size(thisTrainX{ii},1),1);...
        zeros(2200,1)]);
    % Concatenate real and imputed data together
    preTrainX{ii} = cat(1,thisTrainX{ii},newPre{ii});
    % Create Ys for training and testing
    preTrainY{ii} = ones(size(preTrainX{ii},1),1);
    preTestY{ii} = ones(size(preTestX{ii},1),1);
end
% Create data structure 'pre5' for data from first 5 seconds before feeding
pre5.trainX = preTrainX;
pre5.trainY = preTrainY;
pre5.testX = preTestX;
pre5.testY = preTestY;
% Prep test sets for data 40<x<5 seconds before feeding
c = 1;
pre40Test = cell(1,7);
for ii = 2:8
    pre40Test{c} = cat(1,preData{:,ii});
    c = c+1;
end
% Calculate distribution for when determining how many zeros need to be
% added to maintain natural ratio
d = sum(cellfun(@(x) size(x,1),preData(:,2:8)),1)./(sum(cellfun(...
    @(x) size(x,1),preData(:,2:8)),1)+sum(cellfun(@(x) size(x,1),...
    data{1,1}(:,end)))); 
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\',...
    'finalNew\preBingeTrainTest.mat'],'pre5','pre40Test','thisTrainX',...
    'newPre','d')
%% Grab bingeNot data that does not overlap preBinge data 
% Grabs all binge and rest trial data from baseline, dep24, dep48, and chow
% conditions; uses relPow and relCoh.
% [data,samp,~] = collateData(['C:\Users\Pythia\Documents\GreenLab'...
%     '\data\paper2\binge_notbinge\'],{'base';'dep24';'dep48';'chow'},...
%     {'pow','coh'},'trl','rel');
[data,samp,~] = collateData(['C:\Users\Pythia\Documents\GreenLab'...
    '\data\paper2\processed2\'],{'base';'dep24';'dep48';'chow'},...
    {'pow','coh'},'trl','rel');
% Preallocate
nonOverlapInd = cell(4,12);
bingeNotData = cell(size(data));
% Use binge start times to find potential overlap between preBinge
% (bingeStart-40seconds:bingeStart) and notBinge
for h = 1:size(samp,2)
    for ii = 1:size(samp{1,h},1)
        thisCat = [];
        for jj = 1:size(samp{1,h}{ii,1},1)
            thisCat = cat(2,thisCat,samp{1,h}{ii,1}(jj,1)-24000:...
                samp{1,h}{ii,1}(jj,1));
        end
        % Find all binge starts and stops which overlap with preBinge
        % samples (1 = overlap)
        startOverlap = sum(samp{1,h}{ii,2}(:,1) == unique(thisCat),2);
        stopOverlap = sum(samp{1,h}{ii,2}(:,2) == unique(thisCat),2);
        % Get indices of epochs which do not overlap either starts or stops
        nonOverlapInd{h,ii} = logicFind(0,startOverlap+stopOverlap,'==');
        % Subset notBinge data by 'nonOverlapInd'
        bingeNotData{1,h}{ii,1} = data{1,h}{ii,1};
        bingeNotData{1,h}{ii,2} = data{1,h}{ii,2}(nonOverlapInd{h,ii},:);
    end
end
% Rename orignal data, then replace with 'bingeNotData'
rawData = data;
% data = bingeNotData;
% Create 'trlDat' and 'avgDat' from 'bingeNotData'
trlDat = cell(2,size(bingeNotData,2));
avgDat = cell(size(trlDat));
for ii = 1:size(bingeNotData,2)
    for k = 1:size(bingeNotData{1,ii},2)
        % Concatenate all trials together
        trlDat{k,ii} = cat(1,bingeNotData{1,ii}{:,k});
        % First average all trials together per animal, then concatenate
        thisAvg = cellfun(@mean,bingeNotData{1,ii},'UniformOutput',0);
        avgDat{k,ii} = cat(1,thisAvg{:,k});
    end
end
% Save data
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\',...
%     'finalNew\bingeNotData.mat'],'avgDat','trlDat','data','rawData');
%% Grab bingeNot data, CONTAINS preBinge data
[data,~,~] = collateData(['C:\Users\Pythia\Documents\GreenLab'...
    '\data\paper2\binge_notbinge\'],{'base';'dep24';'dep48';'chow'},...
    {'pow','coh'},'trl','rel');
trlDat = cell(2,size(data,2));
for ii = 1:size(data,2)
    for k = 1:size(data{1,ii},2)
        trlDat{k,ii} = cat(1,data{1,ii}{:,k});
    end
end     
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
%     'finalNew\bingeNotPreData.mat'],'trlDat','data');
%% Plot distribution of changes in voracity
% Load voracity data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    '4conditionBingeSize.mat'])
% Calculate change in voracity
vorDiff = [(voracity(:,3)-voracity(:,1))./voracity(:,1);...
    (voracity(:,2)-voracity(:,1))./voracity(:,1)];
figure
scatter(ones(1,24),vorDiff.*100,200,'.k')
set(gca,'XTick',1,'XTickLabel','')
ylabel('Percent Change in Voracity (%)')
%% Find potential noise features
% Load data created above
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\',...
    'bingeNotData.mat'])
% load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\',...
%     'bingeRestData.mat'])
% Load voracity data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    '4conditionBingeSize.mat'])
% Correlate voracity with all features
vorDiff = [(voracity(:,3)-voracity(:,1))./voracity(:,1);...
    (voracity(:,2)-voracity(:,1))./voracity(:,1)];
% Get non-NaN indices
inds24 = logicFind(1,~isnan(voracity(:,2)),'==');
inds48 = logicFind(1,~isnan(voracity(:,3)),'==');
% Do same subtractions for features
dep48Feats = (avgDat{1,3}-avgDat{2,3});
dep24Feats = (avgDat{1,2}-avgDat{2,2});
baseFeats = (avgDat{1,1}-avgDat{2,1});
% Divide by baseline and concatenate features
feats = [(dep48Feats-baseFeats(inds48,:))./baseFeats(inds48,:);...
    (dep24Feats-baseFeats(inds24,:))./baseFeats(inds24,:)];
% Obtain r- and p-values using corrcoef
r = zeros(1,size(feats,2));
p = zeros(1,size(feats,2));
for fi = 1:size(feats,2)
    [thisR,thisP] = corrcoef(vorDiff(~isnan(vorDiff)),feats(:,fi));
    r(fi) = thisR(1,2)^2;
    p(fi) = thisP(1,2);
end
% Get indices of significant p-values to exclude those features
pInds = logicFind(0.05,p,'<=');
% Save data, now with pInds
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\',...
%     'finalNew\bingeNotData.mat'],'avgDat','trlDat','data','rawData',...
%     'pInds');
%% Plot features with significant correlations with changes in voracity
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
    '\finalNew\bingeNotData.mat'])
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
thiscorr = cell(length(pInds),1);
for ii = 1:length(pInds)
    figure
    % Convert features to %
    scatter(100.*vorDiff(~isnan(vorDiff)),100.*feats(:,pInds(ii)),'ok',...
        'Filled')
    lsline
    thiscorr{ii} = fitlm(100.*vorDiff(~isnan(vorDiff)),...
        feats(:,pInds(ii)));
    title(nameVect{pInds(ii)})
    text(0,0,['R^2 = ',num2str(round(thiscorr{ii}.Rsquared.Ordinary,2)),...
        newline,'p = ',num2str(round(thiscorr{ii}.Coefficients.pValue(2)...
        ,3))])
    xlabel('Change in Voracity (gm/ms)')
    ylabel('Percent Change in Feature')
end
%% Plot calories consumed per group
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    '4conditionBingeSize.mat'])
cmb = nchoosek(1:4,2);
stats = cell(1,6);
for ii = 1:6
    % Keep stats for dof and t
    [~,p(ii),~,stats{ii}] = ttest2(bingeCal(:,cmb(ii,1)),...
        bingeCal(:,cmb(ii,2)));
end
p = p.*6;
figure
hold on
for ii = 1:4
    scatter(ones(12,1).*ii,bingeCal(:,ii),200,'.k')
    plot([ii-0.25 ii+0.25],repmat(mean(bingeCal(:,ii),'omitnan'),1,2),'r-')
end
group = mat2cell(cmb,ones(1,6),2);
sigstar(group(p<0.05),p(p<0.05))
ylabel('Calories (kCal)')
set(gca,'XTick',1:4,'XTickLabel',{'Base','dep24','Dep48','Chow'})
%% Set up bingeNotRest data
% Get trialized data
[data,~,~] = collateData(['C:\Users\Pythia\Documents\GreenLab'...
    '\data\paper2\bingeRest\'],{'base';'dep24';'dep48';'chow'},...
    {'pow','coh'},'trl','rel');
trialData = cell(2,size(data,2));
for ii = 1:4
    for jj = 1:2
        trialData{jj,ii} = cat(1,data{ii}{:,jj});
    end
end
% Get average data
[data,~,~] = collateData(['C:\Users\Pythia\Documents\GreenLab'...
    '\data\paper2\bingeRest\'],{'base';'dep24';'dep48';'chow'},...
    {'pow','coh'},'avg','rel');
avgData = cell(2,size(data,2));
for ii = 1:4
    for jj = 1:2
        avgData{jj,ii} = cat(1,data{ii}{:,jj});
    end
end
% Grab pInds from 'bingeNotData.mat'
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\bingeNotData.mat'],'pInds')
% Save
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\bingeRest.mat'],'avgData','trialData','pInds')
%% Baseline binge size
% Run 'runBaseBingeSize.mat' and analyze results.
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseBingeSize\baseBingeSize.mat'])
main = {'Predicting Baseline Binge Size: Binge - Rest',...
    'Predicting Baseline Binge Size: Binge',...
    'Predicting Baseline Binge Size: Rest'};
for ii = 1:3
   doubleHist(real{ii,1}.err,perm{ii,1}.err(1:1000),'Main',main{ii},...
       'xlab','Mean Absolute Error (gm)'); 
end
%% Binge size change
% Run 'runBingeSizeChange.mat' and analyze results.
main = {'Predicting Binge Size Change: Binge - Rest',...
    'Predicting Binge Size Change: Binge',...
    'Predicting Binge Size Change: Rest'};
for ii = 1:3
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\bingeSizeChange\',num2str(ii),'.mat'])
    doubleHist(real{1,ii}.err,perm{1,ii}.err(1:1000),'Main',main{ii},...
       'xlab','Mean Absolute Error (%)'); 
end
%% Food palatability
% Run 'runPalat.mat' and analyze results.
main = {'Predicting Palatability: Binge - Rest',...
    'Predicting Palatability: Binge',...
    'Predicting Palatability: Rest'};
for ii = 1:3
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\palat\palat',num2str(ii),'.mat'])
    doubleHist(1-real.err,1-perm.err(1:1000),'Main',main{ii},...
       'xlab','Accuracy (%)'); 
end
%% Create baseline data set with 500 per animal
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\bingeNotData2.mat'])
% load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
%     '\bingeNotPreData.mat'])
[all,each,rnd] = evenDataSplit(data{1,1},500,100,'ADA',20); %#ok<ASGLU>
save('baseline500Each6000All50-50v2.mat','all','each','rnd','pInds')
% save('baselinePre500Each6000All50-50.mat','all','each','rnd')
%% Create baseline data set with 80/20 split per animal
% load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
%     '\bingeNotData.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\bingeNotPreData.mat'])
[all,each,rnd,weights] = evenDataSplit(data{1,1},0.8,0.2,'',20);%#ok<ASGLU>
% save('baseline500Each6000All.mat','all','each','rnd','weights','pInds')
save('baselinePre500Each6000All.mat','all','each','rnd','weights')
%% Create dep24 data set with 500 per animal
% load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
%     '\bingeNotData.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\bingeNotPreData.mat'])
[all,each,rnd,~] = evenDataSplit(data{1,2},500,100,'ADA',20);%#ok<ASGLU>
% save('dep24_500Each6000All50-50.mat','all','each','rnd','pInds')
save('dep24Pre_500Each6000All50-50.mat','all','each','rnd')
%% Create dep48 data set with 500 per animal
% load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
%     '\bingeNotData.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\bingeNotPreData.mat'])
[all,each,rnd,~] = evenDataSplit(data{1,3},500,100,'ADA',20);%#ok<ASGLU>
save('dep48Pre_500Each6000All50-50.mat','all','each','rnd')
%% Create chow data set with 500 per animal
% load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
%     '\bingeNotData.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\bingeNotPreData.mat'])
[all,each,rnd,~] = evenDataSplit(data{1,4},500,100,'ADA',20);
% save('chow500Each6000All50-50.mat','all','each','rnd','pInds')
save('chowPre500Each6000All50-50.mat','all','each','rnd')
%% Load learningCurve50-50 data
% Concatenate AUCs into matrix  100 (5 repetitions x 20 datasets) x 11
% subsample sets
auc = zeros(100,18);
for ii = 1:20
    load(['C:/Users/Pythia/Documents/GreenLab/data/paper2/analyzed/'...
        'finalNew/learningCurve50-50/',num2str(ii),'.mat'])
    auc(ii*5-4:ii*5,:) = concatData.auc;
end
% Get average and standard deviations of aucs
mA = mean(auc,1);
sA = std(auc,[],1);
% Plot learning curve
samp = [4:2:10,14:8:42,50:50:500];
scatterErr(samp,mA,sA,1)
plot(samp,mA,'k')
ylim([0.6 .9])
set(gca,'YTick',0.6:0.1:.9)
title('Learning Curve: Imputed')
xlabel('Training Set Sample Size (per Animal)')
ylabel('AUC')
% Find when the sample size is not significantly different from the full
% 500 per animal
for ii = 1:17
   [~,p(ii)] = ttest2(auc(:,ii),auc(:,18)); 
end
% Apply Bonferroni correction
p = p.*17;
% Plot square
plot(samp(logicFind(0.05,p,'>','first')),mA(logicFind(0.05,p,'>',...
    'first')),'sq','MarkerSize',9,'Color',[0.5 0.5 0.5])
%% Load learningCurve data weighted
% Concatenate AUCs into matrix  100 (5 repetitions x 20 datasets) x 17
% subsample sets
auc = zeros(100,17);
for ii = 1:20
    load(['C:/Users/Pythia/Documents/GreenLab/data/paper2/analyzed/'...
        'finalNew/learningCurve/',num2str(ii),'.mat'])
    auc(ii*5-4:ii*5,:) = concatData.auc;
end
% Get average and standard deviations of aucs
mA = mean(auc,1);
sA = std(auc,[],1);
% Plot learning curve
samp = [6:2:10,14:8:42,50:50:500];
scatterErr(samp,mA,sA,1)
plot(samp,mA,'k')
ylim([0.6 .9])
set(gca,'YTick',0.6:0.1:.9)
title('Learning Curve: Weighted')
xlabel('Training Set Sample Size (per Animal)')
ylabel('AUC')
% Find when the sample size is not significantly different from the full
% 500 per animal
for ii = 1:16
   [~,p(ii)] = ttest2(auc(:,ii),auc(:,17)); 
end
% Apply Bonferroni correction
p = p.*17;
% Plot square
plot(samp(logicFind(0.05,p,'>','first')),mA(logicFind(0.05,p,'>',...
    'first')),'sq','MarkerSize',9,'Color',[0.5 0.5 0.5])
%% Load learningCurveEach data
rocA = zeros(100,12,12);
for ii = 1:240
    load(['C:/Users/Pythia/Documents/GreenLab/data/paper2/analyzed/'...
       'finalNew/learningCurveEach2/',num2str(ii),'.mat'])
   animal = ceil(ii/20);
   iter = rem(ii,20);
   if iter == 0
       iter = 20;
   end
    rocA(iter*5-4:iter*5,1:12,animal) = eachData.auc;
end
% Get mean and standard deviation of AUCs
mA = squeeze(mean(rocA,1));
% sA = squeeze(std(rocA,[],1));
% Find break point for each animal
for ii = 1:12
   for jj = 1:11
      [~,p(ii,jj)] = ttest2(rocA(:,jj,ii),rocA(:,12,ii)); 
   end
end
% Apply Bonferroni correction
p = p.*132;
% Plot learning curve
samp = [30:20:100,150:50:500];
figure
hold on
for ii = 1:12
     plot(samp,mA(:,ii),'LineWidth',1,'Color',[0.5 0.5 0.5])
     plot(samp(logicFind(0.05,p(ii,:),'>','first')),mA(logicFind(0.05,...
         p(ii,:),'>','first'),ii),'sq','MarkerSize',7,'Color',...
         [0.5 0.5 0.5])
end
scatterErr(samp,mean(mA,2),std(mA,[],2),0)
plot(samp,mean(mA,2),'-k')
set(gca,'YTick',0.6:0.1:1)
title('Learning Curve: Each')
xlabel('Training Set Sample Size')
ylabel('AUC')
%% Load learningCurvePre data
% Concatenate AUCs into matrix  100
auc = zeros(100,16);
for ii = 1:20
    load(['C:/Users/Pythia/Documents/GreenLab/data/paper2/analyzed/'...
        'finalNew/learningCurvePre/',num2str(ii),'.mat'])
    auc(ii*5-4:ii*5,:) = preData.auc;
end
% Get average and standard deviations of aucs
mA = mean(auc,1);
sA = std(auc,[],1);
% Plot learning curve
samp = [50:100:950,1080,1800,2400,3000,3600,4200];
scatterErr(samp./12,mA,sA,1)
plot(samp./12,mA,'k')
ylim([0.6 .9])
set(gca,'YTick',0.6:0.1:.9)
title('Learning Curve: Pre')
xlabel('Training Set Sample Size (per Animal)')
ylabel('AUC')
% Find when the sample size is not significantly different from the full
% 150 per animal
for ii = 1:12
   [~,p(ii)] = ttest2(auc(:,ii),auc(:,13)); 
end
% Apply Bonferroni correction
p = p.*12;
% Plot square
plot(samp(logicFind(0.05,p,'>','first'))./12,mA(logicFind(0.05,p,'>',...
    'first')),'sq','MarkerSize',9,'Color',[0.5 0.5 0.5])
%% Load preBinge files
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\'...
    'preBinge\'])
% Preallocate
[preA,preARand] = deal(zeros(20,8));
% Hardcoded first dimension from know size of roc curves in concatData
[preX,preY] = deal(zeroes(21207,20));
[preRandX,preRandY] = deal(zeros(21222,20));
beta = zeros(20,58);
for ii = 1:20
   load([num2str(ii),'.mat'])
   preA(ii,:) = concatData{1,8}.auc;
   beta(ii,:) = concatData{1,1}.allBeta{1,1}.survBeta;
   preX(:,ii) = concatData{1,1}.acc{1,1}.x;
   preY(:,ii) = concatData{1,1}.acc{1,1}.y;
end
mBeta = mean(beta,1)';
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\'...
    'preBingeRand\'])
for ii = 1:20
   load([num2str(ii),'.mat'])
   preARand(ii,:) = concatData{1,8}.auc;
   preRandX(:,ii) = concatData{1,1}.acc{1,1}.x;
   preRandY(:,ii) = concatData{1,1}.acc{1,1}.y;
end
%% Plot pre->pre
figure
h(1) = plot(mean(preX,2),mean(preY,2),'k');
hold on
h(2) = plot(mean(preRandX,2),mean(preRandY,2),'--k');
h(3) = plot(NaN,NaN,'color','none');
legend(h,{['Actual: ',num2str(round(mean(preA(:,1)),2)),'\pm',...
    num2str(round(conf(preA(:,1)',0.95,'tail',2),2))],...
    ['Permuted: ',num2str(round(mean(preARand(:,1)),2)),'\pm',...
    num2str(round(conf(preARand(:,1)',0.95,'tail',2),2))],...
    ['d = ',num2str(round(distES(preA(:,1),preARand(:,1)),2))]},...
    'Location','SE')
%% Plot preBinge through time
% Get CI and run t-tests
preACI = conf(preA',0.95,'tail',2);
preARandCI = conf(preARand',0.95,'tail',2);
% T-test
for ii = 1:8  
    [~,p(ii)] = ttest2(preA(:,ii),preARand(:,ii));
end
% Apply Bonferroni correction
p = p*8;
% Plot with standard deviation
% scatterErr(1:8,mean(preA,1),std(preA,1),1)
% scatterErr(1:8,mean(preARand,1),std(preARand,1),0)
% Plot with 95% CI
scatterErr(1:8,mean(preA,1),preACI,1)
scatterErr(1:8,mean(preARand,1),preARandCI,0,'Line','--')
% Set up axes to show correct pre time on x and flip y location
set(gca,'XTick',1:8,'XTickLabel',-2.5:-5:-42.5,'XDir',...
    'Reverse','YAxisLocation','right')
% Extend x-axis to include 0
xlim([0.75 8])
xlabel('Time Before (s)')
ylabel('AUC')
title('Pre Binge Prediction')
%% Test f-nF lasso model on preBinge test data
% Load preBinge data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\preBingeTrainTest.mat'])
% Load f-nF data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'])
inds = 1:60;
inds = inds(~ismember(inds,pInds));
% Preallocate
[cPreX,cPreY] = deal(cell(20,8));
cPreA = zeros(20,8);
for ii = 1:20
    % Load model
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\concat\',num2str(ii),'.mat'],'concatData') 
    notFeedTest = all.testX{ii}(all.testY{ii}==0,inds);
    % Figure out number of preBinges to use
    numPre = round((size(notFeedTest,1)-(1-0.006).*size(notFeedTest,1))...
        ./(1-0.006));
    % Set up testing data for first window
    testX = [pre5.testX{1,ii}(randperm(size(pre5.testX{ii},1),numPre),...
        inds);notFeedTest];
    testY = [ones(numPre,1);zeros(size(notFeedTest,1),1)];
    % Test model
    prob = cvglmnetPredict(concatData.model,zscore(testX),['lambda_',...
        concatData.hist.cfg.minTerm],'response');
    [cPreX{ii,1},cPreY{ii,1},~,cPreA(ii,1)] = perfcurve(testY,prob,1);
    % Set up counter starting at 2 to account for first window done above
    c = 2;
    for jj = 1:7
        % Figure out number of preBinges to use
        numPre = round((size(notFeedTest,1)-(1-d(jj))*...
            size(notFeedTest,1))/(1-d(jj)));
        % Get random indices
        rInd = randperm(size(pre40Test{jj},1),numPre);
        % Set up testing data - concatenate preTest data with not feeding
        % data
        testX = [pre40Test{jj}(rInd,inds);notFeedTest];
        testY = [ones(numPre,1);zeros(size(notFeedTest,1),1)];
        % Test on all other windows
        [predY] = cvglmnetPredict(concatData.model,zscore(testX),...
            ['lambda_',concatData.hist.cfg.minTerm],'response');
        [cPreX{ii,c},cPreY{ii,c},~,cPreA(ii,c)] = perfcurve(testY,predY,1);
        c = c+1;
    end
end
% Get average and CI of 'cPreA'
cPreAM = mean(cPreA,1);
cPreACI = conf(cPreA',0.95,'tail',2);
% Plot concat -> pre with pre -> pre (with offset) and rand
scatterErr(1.1:8.1,cPreAM,cPreACI,1,'col',[0.5 0.5 0.5])
scatterErr(0.9:7.9,mean(preA,1),preACI,0)
scatterErr(1:8,mean(preARand,1),preARandCI,0,'Line','--')
title('concat -> pre')
set(gca,'XTick',1:8,'XTickLabel',-2.5:-5:-42.5,'XDir',...
    'Reverse','YAxisLocation','right')
xlim([0.75 8.1])
xlabel('Time Before (s)')
ylabel('AUC')
%% Test preBinge lasso model on f-nF test data
% Load f-nF data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'])
inds = 1:60;
inds = inds(~ismember(inds,pInds));
% Preallocate
[preCX,preCY] = deal(zeros(20,1201));
preCA = zeros(1,20);
for ii = 1:20
    % Load preBinge model
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\preBinge\',num2str(ii),'.mat'])
    testX = all.testX{ii}(:,inds);
    testY = all.testY{ii};
    prob = cvglmnetPredict(concatData{1}.model,zscore(testX),['lambda_',...
        concatData{1}.hist.cfg.minTerm],'response');
    [preCX(ii,:),preCY(ii,:),~,preCA(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/1200:1,'UseNearest',0);
end
% Average ROC
preCXM = mean(preCX,1);
preCYM = mean(preCY,1);
% Average and CI AUC
preCAM = mean(preCA,2);
preCACI = conf(preCA,0.95,'tail',2);
% Get effect size
d = distES(preCA,ccA);
% Compare to concatLasso ROC
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLasso.mat'],'ccXM','ccYM','ccA','ccAM','ccACI')
figure
hold on
h(1) = plot(preCXM,preCYM,'Color',[0.5 0.5 0.5]);
h(2) = plot(ccXM,ccYM,'-k');
h(3) = plot(NaN,NaN,'Color','None');
legend(h,{['Concat: ',num2str(round(ccAM,2)*100),'\pm',...
    num2str(round(ccACI,2)*100),'%'],['Pre: ',...
    num2str(round(preCAM,2)*100),'\pm',num2str(round(preCACI,2)*100),...
    '%'],['d = ',num2str(round(d,2))]},'Location','southeast')
title('preBinge vs. concatLasso')
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
%% Power or Coh: PreBinge-Binge-PostBinge
% Channels: SL,SR,CL,CR; SLSR,SLCL,SLCR,SRCL,SRCR,CLCR
% Freq: delta,theta,alpha,beta,lgamma,hgamma
clear chan pair freq
chan = [];
pair = 2;
freq = 1;
feat = 'd';
loc = 'slcl';
% Get preBinge and notBinge data
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    'preBingeCombined'],'base','in');
% Preallocate
[prePow,preCoh] = deal(zeros(size(files,2),61));
[notPow,notCoh] = deal(zeros(1,size(files,2)));
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'psdTrls')
    else
        load(files{ii},'coh')
    end
   for t = 1:61
       if isempty(pair)
           prePow(ii,t) = mean(psdTrls{t}.relPow(freq,chan,:),'omitnan');
       else
           preCoh(ii,t) = mean(coh{t}.rel(pair,freq,:),'omitnan');
       end
   end
   if isempty(pair)
       notPow(ii) = mean(psdTrls{1,62}.relPow(freq,chan,:),'omitnan');
   else
       notCoh(ii) = mean(coh{1,62}.rel(pair,freq,:));
   end
end

if isempty(pair)
    mPre = mean(prePow,1,'omitnan');
    sPre = std(prePow,[],1,'omitnan');
    mNot = mean(mean(notPow,1,'omitnan'));
else
    mPre = mean(preCoh,1,'omitnan');
    sPre = std(preCoh,[],1,'omitnan');
    mNot = mean(mean(notCoh,1,'omitnan'));
end
% Get Binge data
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    'binge'],'base','in');
% Preallocate
[bingePow,bingeCoh] = deal(zeros(size(files,2),31));
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'psdTrls')
    else
        load(files{ii},'coh')
    end
    for t = 1:31
        if isempty(pair)
            bingePow(ii,t) = mean(psdTrls{t}.relPow(freq,chan,:),...
                'omitnan');
        else
            bingeCoh(ii,t) = mean(coh{t}.rel(pair,freq,:),'omitnan');
        end
    end
end
if isempty(pair)
    mBinge = mean(bingePow,1,'omitnan');
    sBinge = std(bingePow,[],1,'omitnan');
else
    mBinge = mean(bingeCoh,1,'omitnan');
    sBinge = std(bingeCoh,[],1,'omitnan');
end

% Get PostBinge data
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    'postBinge'],'(e','ex','test','ex');
% Preallocate
[postPow,postCoh] = deal(zeros(size(files,2),61));
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'psdTrls')
    else
        load(files{ii},'coh')
    end
    for t = 1:61
        if isempty(pair)
            postPow(ii,t) = mean(psdTrls{t}.relPow(freq,chan,:),'omitnan');
        else
            postCoh(ii,t) = mean(coh{t}.rel(pair,freq,:),'omitnan');
        end
    end
end
if isempty(pair)
    mPost = mean(postPow,1,'omitnan');
    sPost = std(postPow,[],1,'omitnan');
else
    mPost = mean(postCoh,1,'omitnan');
    sPost = std(postCoh,[],1,'omitnan');
end

% Plot
figure
hold on
shadedErrorBar(1:61,fliplr(mPre.*100),fliplr(sPre.*100))
shadedErrorBar(62:92,fliplr(mBinge.*100),fliplr(sBinge.*100),...
    {'color',[0 0.45 0.74]})
shadedErrorBar(102:110,mPost(1:9).*100,sPost(1:9).*100,...
    {'color',[0 0.45 0.74]})
shadedErrorBar(111:162,mPost(10:61).*100,sPost(10:61).*100)
plot(1:92,ones(1,92).*mNot*100,'k')
plot(1:92,ones(1,92).*mean(mBinge)*100,'--','color',[0 0.45 0.74])
plot(102:162,ones(1,61).*mNot*100,'k')
plot(102:162,ones(1,61).*mean(mBinge)*100,'--','color',[0 0.45 0.74])
xlim([1 162])
set(gca,'XTick',[1:10:51,61.5,71:10:91,101,110.5,121:10:162],...
    'XTickLabel',[-62.5:10:-12.5,0,12.5:10:32.5,-12.5,0,12.5:10:52.5])
title([loc,' ',feat]);
ylabel(['% ',feat])
text(162,mean(mBinge)*100,'Binge','color',[0 0.45 0.74])
text(162,mNot*100,'Other')
xlabel('Time')
box off
%% Find approach times associated with binge sessions of all pre-feeding 
% windows
files{1,1} = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\'...
    'paper2\preBingeCombined'],'base');
files{1,2} = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\'...
    'paper2\preBingeCombined'],'dep24');
files{1,3} = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\'...
    'paper2\preBingeCombined'],'dep48');
files{1,4} = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\'...
    'paper2\preBingeCombined'],'chow');
app = cell(1,4);
for ii = 1:size(files,2)
    for jj = 1:size(files{1,ii},2)
        load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
            'preBingeCombined\',files{ii}{jj}],'trls')
        for m = 1:61
            bStart = trls{m}.sampleinfo(:,2)./400;
            parts = strsplit(files{ii}{jj},'_');
            load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\mat\',...
                parts{1},'_',parts{2}],'eventTs')
            sInd = logicFind(1,strcmp(eventTs.label,'Approach (Start)'),...
                '==');
            eInd = logicFind(1,strcmp(eventTs.label,'Approach (End)'),...
                '==');
            for k = 1:size(bStart,1)
                ind = logicFind(1,eventTs.t{1,eInd}+40>bStart(k) & ...
                    eventTs.t{1,eInd}<bStart(k),'==');
                if ~isempty(ind)
                    app{ii}{m,jj}(k) = eventTs.t{1,eInd}(ind(end))-...
                        eventTs.t{1,sInd}(ind(end));
                else
                    app{ii}{m,jj}(k) = NaN;
                end
            end
        end
    end
end
%% Find the percent of trials coming from binge sessions with approaches 
% that overlap with pre-feeding windows
eoi = 5:65;
perc = zeros(1,61);
for ii = 1:61
    this = [];
    for jj = 1:4
        this = [this,(cat(2,app{jj}{ii,:})<eoi(ii) & ...
            cat(2,app{jj}{ii,:})>eoi(ii)-5)]; %#ok<AGROW>
    end
    perc(ii) = sum(this)/numel(this);
end
figure
plot(2.5:62.5,perc.*100,'k')
set(gca,'xtick',2.5:10:62.5,'YTick',0:10:30)
xlim([0 45])
xlabel('Time before Feeding (sec)')
ylabel('Percent of Data from Approach Behavior (%)')
box off
%% Plot concatLog
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLog.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLogRand.mat'])
d = distES(ccLogA,ccLogRandA);
figure
hold on
% fill(ccLogXfill,ccLogYfill,[.8 .8 .8])
h(1) = plot(ccLogXM,ccLogYM,'-k');
% fill(ccLogRandXfill,ccLogRandYfill,[.8 .8 .8])
h(2) = plot(ccLogRandXM,ccLogRandYM,'--k');
h(3) = plot(NaN,NaN,'Color','none');
legend(h,{['Real: ',num2str(round(ccLogAM,2)*100),'\pm',...
    num2str(round(ccLogACI,2))],['Permuted: ',...
    num2str(round(ccLogRandAM,2)),'\pm',...
    num2str(round(ccLogRandACI,2))],...
    ['d = ',num2str(round(d,2))]},'Location','southeast')
ylim([0 1])
xlim([0 1])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Concat->Concat: Logistic')
%% Plot concatLasso
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLasso.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLassoRand.mat'])
d = distES(ccA,ccRandA);
figure
hold on
% fill(ccXfill,ccYfill,[.8 .8 .8])
h(1) = plot(ccXM,ccYM,'-k');
% fill(ccRandXfill,ccRandYfill,[.8 .8 .8])
h(2) = plot(ccRandXM,ccRandYM,'--k');
h(3) = plot(NaN,NaN,'Color','none');
legend(h,{['Real: ',num2str(round(ccAM,2)*100),'\pm',...
    num2str(round(ccACI,2)*100),'%'],['Permuted: ',...
    num2str(round(ccRandAM,2)*100),'\pm',...
    num2str(round(ccRandACI,2)*100),'%'],...
    ['d = ',num2str(round(d,2))]},'Location','southeast')
ylim([0 1])
xlim([0 1])
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Concat->Concat: Lasso')
%% Plot self predictions (diag) - logistic
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\eachLog.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\eachLogRand.mat'])
% Preallocate
[selfLogXM,selfLogYM] = deal(zeros(12,size(eeLogX{1},1)));
% Grab 12 distinguishable colors
cols = distinguishable_colors(12);
% Calculate average ROC curve for each and plot
% Preallocate
[thisX,thisY] = deal(zeros(101,20));
[selfLogRandXM,selfLogRandYM] = deal(zeros(12,101));
for ii = 1:12
   selfLogXM(ii,:) = mean(eeLogX{1,ii}(:,ii,:),3);
   selfLogYM(ii,:) = mean(eeLogY{1,ii}(:,ii,:),3);
   for jj = 1:20
       thisX(:,jj) = eeLogRandX{ii,jj}{ii};
       thisY(:,jj) = eeLogRandY{ii,jj}{ii};
   end
   selfLogRandXM(ii,:) = mean(thisX,2);
   selfLogRandYM(ii,:) = mean(thisY,2);
end
% Get effect size
d = distES(eeSelfLogA,eeSelfLogRandA);
figure
hold on
% Calculate and plot 1 SD fill for average self ROC
% [selfLogXfill,selfLogYfill] = avgFill(selfLogXM',selfLogYM',2,1);
% fill(selfLogXfill,selfLogYfill,[0.8 0.8 0.8])
% Plot self ROC
h(1) = plot(mean(selfLogXM,1),mean(selfLogYM,1),'-k','LineWidth',2);
% Plot rand self ROC
h(2) = plot(mean(selfLogRandXM,1),mean(selfLogRandYM,1),'--k',...
    'LineWidth',2);
% Plot nothing, used for legend space for effect size
h(3) = plot(NaN,NaN,'Color','none');
title('Each to Self: Logistic')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend(h,{['Real: ',num2str(round(eeSelfLogAM,2)*100),'\pm',...
    num2str(round(eeSelfLogACI,2)*100),'%'],['Permuted: ',...
    num2str(round(eeSelfLogRandAM,2)*100),'\pm',...
    num2str(round(eeSelfLogRandACI,2)*100),'%'],...
    ['d = ',num2str(round(d,2))]},'Location','southeast')
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
%% Plot self predictions (diag) - Lasso
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\eachLasso.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\eachLassoRand.mat'])
% Preallocate
[selfXM,selfYM] = deal(zeros(12,size(eeX{1},1)));
[selfRandXM,selfRandYM] = deal(zeros(12,101));
% Calculate average ROC curve for each and plot
for ii = 1:12
   selfXM(ii,:) = mean(eeX{1,ii}(:,ii,:),3);
   selfYM(ii,:) = mean(eeY{1,ii}(:,ii,:),3);
   [thisX,thisY] = deal(zeros(101,20));
   for jj = 1:20
       if isequal(eeRandX{ii,jj}{ii},[0;1])
           thisX(:,jj) = 0:1/100:1;
       else
           thisX(:,jj) = eeRandX{ii,jj}{ii};
       end
       if isequal(eeRandY{ii,jj}{ii},[0;1])
           thisY(:,jj) = 0:1/100:1;
       else
           thisY(:,jj) = eeRandY{ii,jj}{ii};
       end
   end
   selfRandXM(ii,:) = mean(thisX,2);
   selfRandYM(ii,:) = mean(thisY,2);
end
% Get effect size
d = distES(eeSelfA,eeSelfRandA);
figure
hold on
% Calculate and plot 1 SD fill for average self ROC
% [selfXfill,selfYfill] = avgFill(selfXM',selfYM',2,1);
% fill(selfXfill,selfYfill,[0.8 0.8 0.8])
% Plot self ROC
h(1) = plot(mean(selfXM,1),mean(selfYM,1),'-k','LineWidth',2);
% Plot rand self ROC
h(2) = plot(mean(selfRandXM,1),mean(selfRandYM,1),'--k','LineWidth',2);
% Plot nothing, used for legend space for effect size
h(3) = plot(NaN,NaN,'Color','none');
title('Each to Self: Lasso')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend(h,{['Real: ',num2str(round(eeSelfAM,2)*100),'\pm',...
    num2str(round(eeSelfACI,2)*100),'%'],['Permuted: ',...
    num2str(round(eeSelfRandAM,2)*100),'\pm',...
    num2str(round(eeSelfRandACI,2)*100),'%'],...
    ['d = ',num2str(round(d,2))]},'Location','southeast')
text(0.75,0.2,['AUC: ',num2str(round(eeSelfAM,2)*100),'\pm',...
    num2str(round(eeSelfACI,2)*100),'%'])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
%% Plot example of each-self vs. concat-each
% Preallocate
[eeLogXM,eeLogYM] = deal(zeros(size(eeLogX{1},1),12,12));
% Get average each-each
for ii = 1:12
   eeLogXM(:,:,ii) = mean(eeLogX{ii},3);
   eeLogYM(:,:,ii) = mean(eeLogY{ii},3);
end
% Get AUC differences
diffAUC(1,:) = diag(mean(eeLogA,3))-mean(ceLogA,1)';
% Plot both ROCs and area
x = 4;
figure
hold on
fill([eeLogXM(:,x,x);flipud(ceLogXM(x,:)')],[eeLogYM(:,x,x);...
    flipud(ceLogYM(x,:)')],[0.8 0.8 0.8]);
plot(eeLogXM(:,x,x),eeLogYM(:,x,x),'-k','LineWidth',2)
plot(ceLogXM(x,:),ceLogYM(x,:),'--k','LineWidth',2)
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
title('Concat vs. Self: Example')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend({['AUC Diff = ',num2str(round(diffAUC(x),2))],'Self','Concat'})
%% concatLogLOO
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLogLOO\'])
% Preallocate
a = zeros(20,12);
[x,y] = deal(zeros(101,12,20));
for ii = 1:20
   load([num2str(ii),'.mat'])
   a(ii,:) = concatData.auc;
   x(:,:,ii) = cat(2,concatData.rocX{:});
   y(:,:,ii) = cat(2,concatData.rocY{:});
end
% Get AUC differences
diffAUC(2,:) = diag(mean(eeLogA,3))-mean(a,1)';
aM = mean(a,1);
aCI = conf(aM,0.95,'tail',2);
xM = mean(x,3);
yM = mean(y,3);
uberXM = mean(xM,2);
uberYM = mean(yM,2);
%% concatLogLOORand
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLogLOORand\'])
% Preallocate
aRand = zeros(20,12);
[xRand,yRand] = deal(zeros(101,12,20));
for ii = 1:20
   load([num2str(ii),'.mat'])
   aRand(ii,:) = concatData.auc;
   xRand(:,:,ii) = cat(2,concatData.rocX{:});
   yRand(:,:,ii) = cat(2,concatData.rocY{:});
end
aRandM = mean(aRand,1);
aRandCI = conf(aRandM,0.95,'tail',2);
xRandM = mean(xRand,3);
yRandM = mean(yRand,3);
uberXRandM = mean(xRandM,2);
uberYRandM = mean(yRandM,2);
%%
% Get effect size
d = distES(aM,aRandM);
figure
hold on
h(1) = plot(uberXM,uberYM,'-k');
h(2) = plot(uberXRandM,uberYRandM,'--k');
h(3) = plot(NaN,NaN,'Color','none');
legend(h,{['Real: ',num2str(round(mean(aM),2)*100),'\pm',...
    num2str(round(aCI,2)*100),'%'],['Permuted: ',...
    num2str(round(mean(aRandM),2)*100),'\pm',...
    num2str(round(aRandCI,2)*100),'%'],['d = ',num2str(round(d,2))]},...
    'Location','SE')
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Logisitic LOO')
%% Plot distribution of auc differences
figure
scatter(ones(1,12),diffAUC(1,:),400,'.k')
hold on
scatter(ones(1,12).*2,diffAUC(2,:),400,[0.5 0.5 0.5],'.');
plot([0.5 2.5],[0 0],'--k')
for ii = 1:12
    plot([1 2],diffAUC(:,ii))
end
xlim([0.5 2.5])
set(gca,'XTick',1:2,'XTickLabel',{'All','LOO'})
title('AUC Differences: Self-Concat')
ylabel('AUC Difference')
%% All states model - permuted
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\indvGen.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'], 'pInds')
inds = 1:60;
inds = inds(~ismember(inds,pInds));

for ii = 1:20
    trainX = genTrainX{ii}(:,inds);
    trainY = genTrainY{ii}(randperm(length(genTrainY{ii})));
    testX = genTestX{ii}(:,inds);
    testY = genTestY{ii}(randperm(length(genTestY{ii})));
    mdl = fitglm(trainX,trainY,'distribution','binomial') ;
    prob = predict(mdl,testX);
    [~,~,~,rndA(ii)] = perfcurve(testY,prob,1);
end
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
%     '\allPermuted.mat'],'rndA')
%% logCmbs - all
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\logCmbs\all\'])
load('1.mat')
allA{1} = A;
load('2.mat')
allA{2} = A;
allA{3} = [];
for ii = 1:20
   load(['trip',num2str(ii),'.mat'])
   allA{3} = [allA{3};A(ii,:)];
end
monadM = mean(allA{1,1},1)';
dyadM = mean(allA{1,2},1)';
triadM = mean(allA{1,3},1)';
%% Monadic
% Use performance of permuted full logistic for cutoff
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\allPermuted.mat'])
[monadMS,monadInd] = sort(monadM,'descend');
% Sort all monad data by above indices
monad = allA{1,1}(:,monadInd);
for ii = 1:length(monadMS)
   [~,p(ii)] = ttest2(monad(:,ii),rndA); 
end
% Apply Bonferroni correction
p = p.*length(p);
% Find first monad that is not significantly different from permuted
worst = logicFind(0.05,p,'>=','first');
% Find teirs
monadTier = tier(monad);
% Remove teirs below worst and add in worst
monadTier(monadTier>=worst) = [];
monadTier = [monadTier,worst]';
% Get feature names
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
monadFeats = nameVect(monadInd)';
%% Dyadic
% Use performance of permuted full logistic for cutoff
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\allPermuted.mat'])
[dyadMS,dyadInd] = sort(dyadM,'descend');
% Sort all monad data by above indices
dyad = allA{1,2}(:,dyadInd);
for ii = 1:length(dyadMS)
   [~,p(ii)] = ttest2(dyad(:,ii),rndA); 
end
% Apply Bonferroni correction
p = p.*length(p);
% Find first monad that is not significantly different from permuted
worst = logicFind(0.05,p,'>=','first');
% Find teirs
dyadTier = tier(dyad);
% Remove teirs below worst and add in worst
dyadTier(dyadTier>=worst) = [];
dyadTier = [dyadTier,worst]';
% Get feature names
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
cmbs = nchoosek(1:58,2);
for ii = 1:size(cmbs,1)
   dyadFeats(ii,1) = nameVect(cmbs(dyadInd(ii),1));
   dyadFeats(ii,2) = nameVect(cmbs(dyadInd(ii),2));
end
%% Triadic
% Use performance of permuted full logistic for cutoff
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\allPermuted.mat'])
[triadMS,triadInd] = sort(triadM,'descend');
% Sort all monad data by above indices
triad = allA{1,3}(:,triadInd);
for ii = 1:length(triadMS)
   [~,p(ii)] = ttest2(triad(:,ii),rndA); 
end
% Apply Bonferroni correction
p = p.*length(p);
% Find first monad that is not significantly different from permuted
worst = logicFind(0.05,p,'>=','first');
% Find teirs
triadTier = tier(triad);
% Remove teirs below worst and add in worst
triadTier(triadTier>=worst) = [];
triadTier = [triadTier,worst]';
% Get feature names
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
cmbs = nchoosek(1:58,3);
for ii = 1:size(cmbs,1)
   triadFeats(ii,1) = nameVect(cmbs(triadInd(ii),1));
   triadFeats(ii,2) = nameVect(cmbs(triadInd(ii),2));
   triadFeats(ii,3) = nameVect(cmbs(triadInd(ii),3));
end
%% Get frequency of different frequency ranges in monads, dyads, and triads
freqs = {'d','t','a','b','lg','hg'};
for ii = 1:size(freqs,2)
   perc(ii,1) = sum(~cellfun(@isempty,regexp(monadFeats(1:(monadTier(end)-1),:),freqs{ii})))/(monadTier(end)-1);
   perc(ii,2) = sum(any(~cellfun(@isempty,regexp(dyadFeats(1:(dyadTier(end)-1),:),freqs{ii})),2))/(dyadTier(end)-1);
   perc(ii,3) = sum(any(~cellfun(@isempty,regexp(triadFeats(1:(triadTier(end)-1),:),freqs{ii})),2))/(triadTier(end)-1);
   
   topPerc(ii,1) = sum(~cellfun(@isempty,regexp(monadFeats(1:(monadTier(2)-1),:),freqs{ii})))/(monadTier(2)-1);
   topPerc(ii,2) = sum(any(~cellfun(@isempty,regexp(dyadFeats(1:(dyadTier(2)-1),:),freqs{ii})),2))/(dyadTier(2)-1);
   topPerc(ii,3) = sum(any(~cellfun(@isempty,regexp(triadFeats(1:(triadTier(2)-1),:),freqs{ii})),2))/(triadTier(2)-1);
end
% figure
% bar(topPerc'.*100,'stacked')
% ylabel('Cumulative Percent of Models (%)');
% set(gca,'XTickLabel',{'Monad','Dyad','Triad'})
% box off

%% logCmbs - base
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\logCmbs\base\'])
% Preallocate
[baseAuc,dep24Auc,dep48Auc,chowAuc] = deal(cell(1,2));
for ii = 1:2
    load([num2str(ii),'.mat'])
    baseAuc{ii} = A;
    dep24Auc{ii} = dep24A;
    dep48Auc{ii} = dep48A;
    chowAuc{ii} = chowA;
end
baseAuc{3} = [];
dep24Auc{3} = [];
dep48Auc{3} = [];
chowAuc{3} = [];
% Only take ii'th row to correct for stupid mistake in storing aucs
for ii = 1:20
    load(['trip',num2str(ii),'.mat'])
    baseAuc{3} = [baseAuc{3};A(ii,:)];
    dep24Auc{3} = [dep24Auc{3};dep24A(ii,:)];
    dep48Auc{3} = [dep48Auc{3};dep48A(ii,:)];
    chowAuc{3} = [chowAuc{3};chowA(ii,:)];
end
% Get means
baseM = cellfun(@(x) mean(x,1),baseAuc,'UniformOutput',0);
dep24M = cellfun(@(x) mean(x,1),dep24Auc,'UniformOutput',0);
dep48M = cellfun(@(x) mean(x,1),dep48Auc,'UniformOutput',0);
chowM = cellfun(@(x) mean(x,1),chowAuc,'UniformOutput',0);
% Get CIs
baseCI = cellfun(@(x) conf(x',0.95,'tail',2),baseAuc,'UniformOutput',0);
dep24CI = cellfun(@(x) conf(x',0.95,'tail',2),dep24Auc,'UniformOutput',0);
dep48CI = cellfun(@(x) conf(x',0.95,'tail',2),dep48Auc,'UniformOutput',0);
chowCI = cellfun(@(x) conf(x',0.95,'tail',2),chowAuc,'UniformOutput',0);
%% Get top 58 models for base cmbs
% Preallocate
cmbs = cell(1,3);
[topA,topADep24,topADep48,topAChow,topCI,topCIDep24,topCIDep48,...
    topCIChow] = deal(zeros(58,3));
topInd = cell(58,3);
sortInd = cell(1,3);
for ii  = 1:3
    cmbs{ii} = nchoosek(1:58,ii);
    % Sort each condition by baseline performance
    [sortX,sortInd{ii}] = sort(baseM{1,ii},'descend');
    sortXDep24 = dep24M{1,ii}(sortInd{ii});
    sortXDep48 = dep48M{1,ii}(sortInd{ii});
    sortXChow = chowM{1,ii}(sortInd{ii});
%     [sortXDep24,sortIndDep24] = sort(dep24M{1,ii},'descend');
%     [sortXDep48,sortIndDep48] = sort(dep48M{1,ii},'descend');
%     [sortXChow,sortIndChow] = sort(chowM{1,ii},'descend');
    % Get top performers
    topA(:,ii) = sortX(1:58);
    topADep24(:,ii) = sortXDep24(1:58);
    topADep48(:,ii) = sortXDep48(1:58);
    topAChow(:,ii) = sortXChow(1:58);
    
    topCI(:,ii) = baseCI{1,ii}(sortInd{ii}(1:58));
    topCIDep24(:,ii) = dep24CI{1,ii}(sortInd{ii}(1:58));
    topCIDep48(:,ii) = dep48CI{1,ii}(sortInd{ii}(1:58));
    topCIChow(:,ii) = chowCI{1,ii}(sortInd{ii}(1:58));
%     topCIDep24(:,ii) = dep24CI{1,ii}(sortIndDep24(1:58));
%     topCIDep48(:,ii) = dep48CI{1,ii}(sortIndDep48(1:58));
%     topCIChow(:,ii) = chowCI{1,ii}(sortIndChow(1:58));
    for jj = 1:58
        topInd{jj,ii} = cmbs{ii}(sortInd{ii}(jj),:);
    end
end
% Get variable names
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'],'pInds')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
nameVect = nameVect(inds);
% Preallocate
monad = cell(58,1);
dyad = cell(58,2);
triad = cell(58,3);
for ii = 1:58
    monad{ii,1} = nameVect{topInd{ii,1}};
    for jj = 1:2
        dyad{ii,jj} = nameVect{topInd{ii,2}(jj)};
    end
    for jj = 1:3
        triad{ii,jj} = nameVect{topInd{ii,3}(jj)};
    end
end
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
%     '\logCmbs.mat'],'baseM','baseCI','dep24M','dep24CI','dep48M',...
%     'dep48CI','chowM','chowCI','topA','topCI','monad','dyad','triad')
%% Plot baseline model 'topA' with 'topCI'
% Load performance of baseline concat logistic model
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLog.mat'],'ccLogAM')
% Load betas from conact lasso
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLasso.mat'],'betaM')
% Get indcies of betas that are always in model
lassoInd = logicFind(1,betaM(sortInd{1}),'==');
% Load logCmb performance
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\logCmbs.mat'],'topA','topCI')
% Uni
scatterErr(1:58,topA(:,1),topCI(:,1),1)
% Highlight lasso indices
scatterErr(lassoInd,topA(lassoInd,1),topCI(lassoInd,1),0,'col','r')
hold on
% Plot concat performance
plot([0 60],[ccLogAM ccLogAM],'-k')
% Plot .5 line
plot([0 60],[0.5 0.5],'--k')
title('Monadic Logistic')
ylabel('AUC')
xlabel('Feature')
% Dyad
scatterErr(1:58,topA(:,2),topCI(:,2),1)
hold on
% Plot concat performance
plot([0 60],[ccLogAM ccLogAM],'-k')
% Plot .5 line
% plot([0 60],[0.5 0.5],'--k')
title('Dyadic Logistic')
ylabel('AUC')
xlabel('Feature')
% Triad
scatterErr(1:58,topA(:,3),topCI(:,3),1)
hold on
% Plot concat performance
plot([0 60],[ccLogAM ccLogAM],'-k')
% Plot .5 line
% plot([0 60],[0.5 0.5],'--k')
title('Triadic Logistic')
ylabel('AUC')
xlabel('Feature')
%% Plot performance across conditions
% Concatenate data
means = [topA(1,1),topADep24(1,1),topADep48(1,1),topAChow(1,1),topA(1,2)...
    ,topADep24(1,2),topADep48(1,2),topAChow(1,2),topA(1,3),...
    topADep24(1,3),topADep48(1,3),topAChow(1,3)];
cis = [topCI(1,1),topCIDep24(1,1),topCIDep48(1,1),topCIChow(1,1),...
    topCI(1,2),topCIDep24(1,2),topCIDep48(1,2),topCIChow(1,2),topCI(1,3)...
    ,topCIDep24(1,3),topCIDep48(1,3),topCIChow(1,3)];
scatterErr(1:12,means,cis,1)
set(gca,'XTick',1:12,'XTickLabel',repmat({'Base','Dep24','Dep48','Chow'}...
    ,1,3))
xtickangle(90)
title('Logistics Across Condition')
xlabel('Test Set')
ylabel('AUC')
%% Plot log, lasso, and logCmbs performance together
% Load logistic performance
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLog.mat'],'ccLogAM','ccLogACI')
% Load lasso performance
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLasso.mat'],'ccAM','ccACI','featM')
% Load logCmbs performance
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\logCmbs.mat'],'topA','topCI')
% Combine data for one-stop shopping
mData = [topA(1,:),ccAM,ccLogAM];
ciData = [topCI(1,:),ccACI,ccLogACI];
% Plot
scatterErr(1:5,mData,ciData,1)
% Alter x axis
set(gca,'XTickLabel',{'Monad','Dyad','Triad',['Lasso (',...
    num2str(round(featM)),')'],'Logistic'},'XTick',1:5)
% Add a little space before first point
xlim([0.75 5])
title('Performance vs. Feature #')
ylabel('AUC')
%% Apply concatLog baseline model to each other state
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'],'pInds')
% Get usable indices
inds = 1:60;
inds = inds(~ismember(inds,pInds));
% Preallocate
[dep24X,dep24Y,dep48X,dep48Y,chowX,chowY] = deal(cell(1,12));
[dep24A,dep48A,chowA] = deal(zeros(1,12));
% Load dep24 data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\dep24_500Each6000All50-50.mat'])
for ii = 1:20
    % Load logistic model
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\concatLog\',num2str(ii),'.mat'])
    % Create 'testX' and 'testY'
    testX = all.testX{ii}(:,inds);
    testY = all.testY{ii};
    % Calculate probabilities
    prob = predict(concatData.model,testX);
    % Get ROC curve data
    [dep24X{ii},dep24Y{ii},~,dep24A(ii)] = perfcurve(testY,prob,1);
    % Clear out 'testX' and 'testY'
    clear testX testY
end
% Load dep48 data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\dep48_500Each6000All50-50.mat'])
for ii = 1:20
    % Load logistic model
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\concatLog\',num2str(ii),'.mat'])
    % Create 'testX' and 'testY'
    testX = all.testX{ii}(:,inds);
    testY = all.testY{ii};
    % Calculate probabilities
    prob = predict(concatData.model,testX);
    % Get ROC curve data
    [dep48X{ii},dep48Y{ii},~,dep48A(ii)] = perfcurve(testY,prob,1);
    % Clear out 'testX' and 'testY'
    clear testX testY
end
% Load chow data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\chow500Each6000All50-50.mat'])
for ii = 1:20
    % Load logistic model
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\concatLog\',num2str(ii),'.mat'])
    % Create 'testX' and 'testY'
    testX = all.testX{ii}(:,inds);
    testY = all.testY{ii};
    % Calculate probabilities
    prob = predict(concatData.model,testX);
    % Get ROC curve data
    [chowX{ii},chowY{ii},~,chowA(ii)] = perfcurve(testY,prob,1);
    % Clear out 'testX' and 'testY'
    clear testX testY
end
%% Prepare for plotting
% Get average ROC curve
dep24XM = mean(cat(2,dep24X{:}),2);
dep24YM = mean(cat(2,dep24Y{:}),2);
% Get ROC 1 SD fill
[dep24Xfill,dep24Yfill] = avgFill(cat(2,dep24X{:}),cat(2,dep24Y{:}),2);
% Get 95% CI of AUC
dep24ACI = conf(dep24A,0.95,'tail',2);

% Get average ROC curve
dep48XM = mean(cat(2,dep48X{:}),2);
dep48YM = mean(cat(2,dep48Y{:}),2);
% Get ROC 1 SD fill
[dep48Xfill,dep48Yfill] = avgFill(cat(2,dep48X{:}),cat(2,dep48Y{:}),2);
% Get 95% CI of AUC
dep48ACI = conf(dep48A,0.95,'tail',2);

% Get average ROC curve
chowXM = mean(cat(2,chowX{:}),2);
chowYM = mean(cat(2,chowY{:}),2);
% Get ROC 1 SD fill
[chowXfill,chowYfill] = avgFill(cat(2,chowX{:}),cat(2,chowY{:}),2);
% Get 95% CI of AUC
chowACI = conf(chowA,0.95,'tail',2);
%% Plot
figure
hold on
% Add concat-concat line
plot(ccLogXM,ccLogYM,'-k','LineWidth',2)
% fill(dep24Xfill,dep24Yfill,'b','FaceAlpha',0.5)
plot(dep24XM,dep24YM,'-k','LineWidth',1)
% fill(dep48Xfill,dep48Yfill,'r','FaceAlpha',0.5)
plot(dep48XM,dep48YM,'--k','LineWidth',1.5)
% fill(chowXfill,chowYfill,'y','FaceAlpha',0.5)
plot(chowXM,chowYM,':k','LineWidth',1.5)
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend({['Base: ',num2str(round(mean(ccLogA),2)),'\pm',...
    num2str(round(dep24ACI,3)*100),'%'],...
    ['Dep24: ',num2str(round(mean(dep24A),2)),'\pm',...
    num2str(round(dep24ACI,3)*100),'%'],...
    ['Dep48: ',num2str(round(mean(dep48A),2)),'\pm',...
    num2str(round(dep48ACI,3)*100),'%'],...
    ['Chow: ',num2str(round(mean(chowA),2)),'\pm',...
    num2str(round(chowACI,3)*100),'%']},'Location','southeast')
title('Baseline Model across Conditions')
ylim([0 1])
xlim([0 1])
%% Apply eachLog baseline model to each other state
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'],'pInds')
% Get usable indices
inds = 1:60;
inds = inds(~ismember(inds,pInds));
% Preallocate
[dep24EachX,dep24EachY] = deal(cell(20,12));
[dep48EachX,dep48EachY,chowEachX,chowEachY] = deal(cell(20,9));
dep24EachA = zeros(20,12);
[dep48EachA,chowEachA] = deal(zeros(20,9));
% Load dep24 data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\dep24_500Each6000All50-50.mat'])
animal = 1;
iter = 1;
for ii = 1:240
    % Load logistic model
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\eachLog\',num2str(ii),'.mat'])
    % Create 'testX' and 'testY'
    testX = each.testX{iter,animal}(:,inds);
    testY = each.testY{iter,animal};
    % Calculate probabilities
    prob = predict(selfData.model,testX);
    % Get ROC curve data
    [dep24EachX{iter,animal},dep24EachY{iter,animal},~,...
        dep24EachA(iter,animal)] = perfcurve(testY,prob,1);
    % Clear out 'testX' and 'testY'
    clear testX testY
    % Add one to 'iter'
    iter = iter + 1;
    % If 'iter' will be 21, then move 'animal' forward one and reset 'iter'
    if iter == 21
        animal = animal + 1;
        iter = 1;
    end
end
% Load dep48 data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\dep48_500Each6000All50-50.mat'])
% Make sure to only load those models for which data exist in other states;
% corresponding to [1:3,5:7,9:11]
animal = 1;
iter = 1;
for ii = [1:60,81:140,161:220]
    % Load logistic model
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\eachLog\',num2str(ii),'.mat'])
    % Create 'testX' and 'testY'
    testX = each.testX{iter,animal}(:,inds);
    testY = each.testY{iter,animal};
    % Calculate probabilities
    prob = predict(selfData.model,testX);
    % Get ROC curve data
    [dep48EachX{iter,animal},dep48EachY{iter,animal},~,...
        dep48EachA(iter,animal)] = perfcurve(testY,prob,1);
    % Clear out 'testX' and 'testY'
    clear testX testY
    % Add one to 'iter'
    iter = iter + 1;
    % If 'iter' will be 21, then move 'animal' forward one and reset 'iter'
    if iter == 21
        animal = animal + 1;
        iter = 1;
    end
end
% Load chow data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\chow500Each6000All50-50.mat'])
% Make sure to only load those models for which data exist in other states;
% corresponding to [1:3,5:8,10:11]
animal = 1;
iter = 1;
for ii = [1:60,81:160,181:220]
    % Load logistic model
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\eachLog\',num2str(ii),'.mat'])
    % Create 'testX' and 'testY'
    testX = each.testX{iter,animal}(:,inds);
    testY = each.testY{iter,animal};
    % Calculate probabilities
    prob = predict(selfData.model,testX);
    % Get ROC curve data
    [chowEachX{iter,animal},chowEachY{iter,animal},~,...
        chowEachA(iter,animal)] = perfcurve(testY,prob,1);
    % Clear out 'testX' and 'testY'
    clear testX testY
    % Add one to 'iter'
    iter = iter + 1;
    % If 'iter' will be 21, then move 'animal' forward one and reset 'iter'
    if iter == 21
        animal = animal + 1;
        iter = 1;
    end
end
%% Prepare for plotting
% Get average AUC
dep24EachAM = mean(dep24EachA,1);
dep48EachAM = mean(dep48EachA,1);
chowEachAM = mean(chowEachA,1);
% Preallocate
[dep24EachACI,dep48EachACI,chowEachACI] = deal(zeros(1,size(dep24EachA,...
    2)));
% Get 95% CI of AUC
for ii = 1:size(dep24EachA,2)
    dep24EachACI(ii) = conf(dep24EachA(:,ii),0.95,'tail',2);
end
for ii = 1:size(dep48EachA,2)
    dep48EachACI(ii) = conf(dep48EachA(:,ii),0.95,'tail',2);
end
for ii = 1:size(chowEachA,2)
    chowEachACI(ii) = conf(chowEachA(:,ii),0.95,'tail',2);
end
%% Prep all state data
% Preallocate
[genTrainX,genTrainY,genTestX,genTestY] = deal(cell(1,20));
% Prep animal indices that have data across all states (exact indices will
% change depending on state due to drop out)
inds = [1:3,5:7,10,11];
% Load baseline data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'])
% Pull out individual data
indTrainX = each.trainX(:,inds);
indTrainY = each.trainY(:,inds);
indTestX = each.testX(:,inds);
indTestY = each.testY(:,inds);
% Concatenate into generalized data
for ii = 1:20
    genTrainX{ii} = cat(1,each.trainX{ii,inds});
    genTrainY{ii} = cat(1,each.trainY{ii,inds});
    genTestX{ii} = cat(1,each.testX{ii,inds});
    genTestY{ii} = cat(1,each.testY{ii,inds});
end
% % Pull out 125 samples from each animal
% permInd = randi(size(each.trainX{1,1},1),125,20);
% genTrainX = cell(1,20);
% genTrainY = cell(1,20);
% for ii = 1:20
%     for jj = 1:8
%        genTrainX{ii} = cat(1,genTrainX{ii},...
%            each.trainX{ii,inds(jj)}(permInd(:,ii),:));
%        genTrainY{ii} = cat(1,genTrainY{ii},...
%            each.trainY{ii,inds(jj)}(permInd(:,ii),:));
%     end
% end

% Load dep24 data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\dep24_500Each6000All50-50.mat'])
% Prep indices
inds = [1:3,5:7,10,11];
% Generate random indices
% permInd = randi(size(each.trainX{1,1},1),125,20);
for ii = 1:20
    for jj = 1:8
        indTrainX{ii,jj} = cat(1,indTrainX{ii,jj},each.trainX{ii,jj});
        indTrainY{ii,jj} = cat(1,indTrainY{ii,jj},each.trainY{ii,jj});
        indTestX{ii,jj} = cat(1,indTestX{ii,jj},each.testX{ii,jj});
        indTestY{ii,jj} = cat(1,indTestY{ii,jj},each.testY{ii,jj});
        % Random sample
%         genTrainX{ii} = cat(1,genTrainX{ii},...
%             each.trainX{ii,inds(jj)}(permInd(:,ii),:));
%         genTrainY{ii} = cat(1,genTrainY{ii},...
%             each.trainY{ii,inds(jj)}(permInd(:,ii),:));
    end
    genTrainX{ii} = cat(1,genTrainX{ii},cat(1,each.trainX{ii,inds}));
    genTrainY{ii} = cat(1,genTrainY{ii},cat(1,each.trainY{ii,inds}));
    genTestX{ii} = cat(1,genTestX{ii},cat(1,each.testX{ii,inds}));
    genTestY{ii} = cat(1,genTestY{ii},cat(1,each.testY{ii,inds}));
end
% Load dep48 data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\dep48_500Each6000All50-50.mat'])
% Prep indices
inds = [1:6,8,9];
% Generate random indices
% permInd = randi(size(each.trainX{1,1},1),125,20);
for ii = 1:20
    for jj = 1:8
        indTrainX{ii,jj} = cat(1,indTrainX{ii,jj},each.trainX{ii,jj});
        indTrainY{ii,jj} = cat(1,indTrainY{ii,jj},each.trainY{ii,jj});
        indTestX{ii,jj} = cat(1,indTestX{ii,jj},each.testX{ii,jj});
        indTestY{ii,jj} = cat(1,indTestY{ii,jj},each.testY{ii,jj});
        % Random sample
%         genTrainX{ii} = cat(1,genTrainX{ii},...
%             each.trainX{ii,inds(jj)}(permInd(:,ii),:));
%         genTrainY{ii} = cat(1,genTrainY{ii},...
%             each.trainY{ii,inds(jj)}(permInd(:,ii),:));
    end
    genTrainX{ii} = cat(1,genTrainX{ii},cat(1,each.trainX{ii,inds}));
    genTrainY{ii} = cat(1,genTrainY{ii},cat(1,each.trainY{ii,inds}));
    genTestX{ii} = cat(1,genTestX{ii},cat(1,each.testX{ii,inds}));
    genTestY{ii} = cat(1,genTestY{ii},cat(1,each.testY{ii,inds}));
end
% Load chow data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\chow500Each6000All50-50.mat'])
% Prep indices
inds = [1:6,8,9];
% Generate random indices
% permInd = randi(size(each.trainX{1,1},1),125,20);
for ii = 1:20
    for jj = 1:8
        indTrainX{ii,jj} = cat(1,indTrainX{ii,jj},each.trainX{ii,jj});
        indTrainY{ii,jj} = cat(1,indTrainY{ii,jj},each.trainY{ii,jj});
        indTestX{ii,jj} = cat(1,indTestX{ii,jj},each.testX{ii,jj});
        indTestY{ii,jj} = cat(1,indTestY{ii,jj},each.testY{ii,jj});
        % Random sample
%         genTrainX{ii} = cat(1,genTrainX{ii},...
%             each.trainX{ii,inds(jj)}(permInd(:,ii),:));
%         genTrainY{ii} = cat(1,genTrainY{ii},...
%             each.trainY{ii,inds(jj)}(permInd(:,ii),:));
    end
    genTrainX{ii} = cat(1,genTrainX{ii},cat(1,each.trainX{ii,inds}));
    genTrainY{ii} = cat(1,genTrainY{ii},cat(1,each.trainY{ii,inds}));
    genTestX{ii} = cat(1,genTestX{ii},cat(1,each.testX{ii,inds}));
    genTestY{ii} = cat(1,genTestY{ii},cat(1,each.testY{ii,inds}));
end
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
%     'finalNew\indvGen.mat'],'genTestX','genTestY','genTrainX',...
%     'genTrainY','indTestX','indTestY','indTrainX','indTrainY')
%% Build baseline pop model, test on other states (indvGen); pop->pop
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\indvGen.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'],'pInds')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
% Preallocate
[popBaseBaseA,popBaseDep24A,popBaseDep48A,popBaseChowA] = deal(zeros(1,...
    20));
for ii = 1:20
    % Build baseline model
    mdl = fitglm(genTrainX{ii}(1:4000,inds),genTrainY{ii}(1:4000),...
        'distribution','binomial');
    % Base
    prob = predict(mdl,genTestX{ii}(1:800,inds));
    [~,~,~,popBaseBaseA(ii)] = perfcurve(genTestY{ii}(1:800),prob,1);
    % Dep24
    prob = predict(mdl,genTestX{ii}(801:1600,inds));
    [~,~,~,popBaseDep24A(ii)] = perfcurve(genTestY{ii}(801:1600),prob,1);
    % Dep48
    prob = predict(mdl,genTestX{ii}(1601:2400,inds));
    [~,~,~,popBaseDep48A(ii)] = perfcurve(genTestY{ii}(1601:2400),prob,1);
    % Chow
    prob = predict(mdl,genTestX{ii}(2401:3200,inds));
    [~,~,~,popBaseChowA(ii)] = perfcurve(genTestY{ii}(2401:3200),prob,1);
end
%% Build baseline ind model, test on other states (indvGen); ind->ind
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\indvGen.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'],'pInds')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
% Preallocate
[indBaseBaseA,indBaseDep24A,indBaseDep48A,indBaseChowA] = ...
    deal(zeros(8,8,20));
for ii = 1:20
    for jj = 1:8
        % Build baseline model
        mdl = fitglm(indTrainX{ii,jj}(1:500,inds),...
            indTrainY{ii,jj}(1:500),'distribution','binomial');
        % Test on each
        for k = 1:8
            % Base
            prob = predict(mdl,indTestX{ii,jj}(1:100,inds));
            [~,~,~,indBaseBaseA(jj,k,ii)] = perfcurve(...
                indTestY{ii,k}(1:100),prob,1);
            % Dep24
            prob = predict(mdl,indTestX{ii,jj}(101:200,inds));
            [~,~,~,indBaseDep24A(jj,k,ii)] = perfcurve(...
                indTestY{ii,k}(101:200),prob,1);
            % Dep48
            prob = predict(mdl,indTestX{ii,jj}(201:300,inds));
            [~,~,~,indBaseDep48A(jj,k,ii)] = perfcurve(...
                indTestY{ii,k}(201:300),prob,1);
            % Chow
            prob = predict(mdl,indTestX{ii,jj}(301:400,inds));
            [~,~,~,indBaseChowA(jj,k,ii)] = perfcurve(...
                indTestY{ii,k}(301:400),prob,1);
        end
    end
end
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
%     '\baseAllStates.mat'],'popBaseBaseA','popBaseDep24A','popBaseDep48A'...
%     ,'popBaseChowA','indBaseBaseA','indBaseDep24A','indBaseDep48A',...
%     'indBaseChowA')
%% Load models built from all states both pop and ind
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\'...
    'allEach'])
% Preallocate
[popAllA,popBaseA,popDep24A,popDep48A,popChowA] = deal(zeros(1,20));
[indAllA,indBaseA,indDep24A,indDep48A,indChowA] = deal(zeros(8,8,20));
for ii = 1:20
   load([num2str(ii),'.mat'])
   % Population - pop->pop
   popAllA(ii) = all.PopA;
%    popAllX(ii,:) = all.PopX;
%    popAllY(ii,:) = all.PopY;
   popBaseA(ii) = base.A;
%    popBaseX(ii,:) = base.X;
%    popBaseY(ii,:) = base.Y;
   popDep24A(ii) = dep24.A;
%    popDep24X(ii,:) = dep24.X;
%    popDep24Y(ii,:) = dep24.Y;
   popDep48A(ii) = dep48.A;
%    popDep48X(ii,:) = dep48.X;
%    popDep48Y(ii,:) = dep48.Y;
   popChowA(ii) = chow.A;
%    popChowX(ii,:) = chow.X;
%    popChowY(ii,:) = chow.Y;

   % Individual - ind->ind
   indAllA(:,:,ii) = indAll.A;
%    indAllX(ii,:) = indAll.X;
%    indAllY(ii,:) = indAll.Y;
   indBaseA(:,:,ii) = indBase.A;
%    indBaseX(ii,:) = indBase.X;
%    indBaseY(ii,:) = indBase.Y;
   indDep24A(:,:,ii) = indDep24.A;
%    indDep24X(ii,:) = indDep24.X;
%    indDep24Y(ii,:) = indDep24.Y;
   indDep48A(:,:,ii) = indDep48.A;
%    indDep48X(ii,:) = indDep48.X;
%    indDep48Y(ii,:) = indDep48.Y;
   indChowA(:,:,ii) = indChow.A;
%    indChowX(ii,:) = indChow.X;
%    indChowY(ii,:) = indChow.Y;
end
%% Plot scatterErr of AUCs
% Models built from baseline
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseAllStates.mat'])
% Population - baseline means
baseBasePopAM = mean(popBaseBaseA);
dep24BasePopAM = mean(popBaseDep24A);
dep48BasePopAM = mean(popBaseDep48A);
chowBasePopAM = mean(popBaseChowA);
% Individual - baseline means
% Preallocate
[baseBaseDiag,baseDep24Diag,baseDep48Diag,baseChowDiag] = ...
    deal(zeros(20,8));
for ii =1:20
   baseBaseDiag(ii,:) = diag(indBaseBaseA(:,:,ii));
   baseDep24Diag(ii,:) = diag(indBaseDep24A(:,:,ii));
   baseDep48Diag(ii,:) = diag(indBaseDep48A(:,:,ii));
   baseChowDiag(ii,:) = diag(indBaseChowA(:,:,ii));
end
baseBaseIndAM = mean(mean(baseBaseDiag));
dep24BaseIndAM = mean(mean(baseDep24Diag));
dep48BaseIndAM = mean(mean(baseDep48Diag));
chowBaseIndAM = mean(mean(baseChowDiag));
% baseBaseIndAM = mean(mean(indBaseBaseA));
% dep24BaseIndAM = mean(mean(indBaseDep24A));
% dep48BaseIndAM = mean(mean(indBaseDep48A));
% chowBaseIndAM = mean(mean(indBaseChowA));
% Concatenate means
baseM = [baseBasePopAM,baseBaseIndAM,dep24BasePopAM,dep24BaseIndAM,...
    dep48BasePopAM,dep48BaseIndAM,chowBasePopAM,chowBaseIndAM];
% Population - baseline confidence intervals
baseBasePopACI = conf(popBaseBaseA,.95,'tail',2);
dep24BasePopACI = conf(popBaseDep24A,.95,'tail',2);
dep48BasePopACI = conf(popBaseDep48A,.95,'tail',2);
chowBasePopACI = conf(popBaseChowA,.95,'tail',2);
% Individual - baseline confidence intervals
baseBaseIndACI = conf(mean(baseBaseDiag,1),.95,'tail',2);
dep24BaseIndACI = conf(mean(baseDep24Diag,1),.95,'tail',2);
dep48BaseIndACI = conf(mean(baseDep48Diag,1),.95,'tail',2);
chowBaseIndACI = conf(mean(baseChowDiag,1),.95,'tail',2);
% baseBaseIndACI = conf(mean(indBaseBaseA,1),.95,'tail',2);
% dep24BaseIndACI = conf(mean(indBaseDep24A,1),.95,'tail',2);
% dep48BaseIndACI = conf(mean(indBaseDep48A,1),.95,'tail',2);
% chowBaseIndACI = conf(mean(indBaseChowA,1),.95,'tail',2);
% Concatenate confidence intervals
baseCI = [baseBasePopACI,baseBaseIndACI,dep24BasePopACI,dep24BaseIndACI,...
    dep48BasePopACI,dep48BaseIndACI,chowBasePopACI,chowBaseIndACI];

% Models built from all states
% Population - all means
baseAllPopAM = mean(popBaseA);
dep24AllPopAM = mean(popDep24A);
dep48AllPopAM = mean(popDep48A);
chowAllPopAM = mean(popChowA);
% Individual - all means (just grab self, i.e. diagonal)
% Preallocate
[baseDiag,dep24Diag,dep48Diag,chowDiag] = deal(zeros(20,8));
for ii = 1:20
    baseDiag(ii,:) = diag(indBaseA(:,:,ii));
    dep24Diag(ii,:) = diag(indDep24A(:,:,ii));
    dep48Diag(ii,:) = diag(indDep48A(:,:,ii));
    chowDiag(ii,:) = diag(indChowA(:,:,ii));
end
baseAllIndAM = mean(mean(baseDiag,1));
dep24AllIndAM = mean(mean(dep24Diag,1));
dep48AllIndAM = mean(mean(dep48Diag,1));
chowAllIndAM = mean(mean(chowDiag,1));
% Concatenate means
allM = [baseAllPopAM,baseAllIndAM,dep24AllPopAM,dep24AllIndAM,...
    dep48AllPopAM,dep48AllIndAM,chowAllPopAM,chowAllIndAM];
% Population - baseline confidence intervals
baseAllPopACI = conf(popBaseA,.95,'tail',2);
dep24AllPopACI = conf(popDep24A,.95,'tail',2);
dep48AllPopACI = conf(popDep48A,.95,'tail',2);
chowAllPopACI = conf(popChowA,.95,'tail',2);
% Individual - baseline confidence intervals
baseAllIndACI = conf(mean(baseDiag,1),.95,'tail',2);
dep24AllIndACI = conf(mean(dep24Diag,1),.95,'tail',2);
dep48AllIndACI = conf(mean(dep48Diag,1),.95,'tail',2);
chowAllIndACI = conf(mean(chowDiag,1),.95,'tail',2);
% Concatenate confidence intervals
allCI = [baseAllPopACI,baseAllIndACI,dep24AllPopACI,dep24AllIndACI,...
    dep48AllPopACI,dep48AllIndACI,chowAllPopACI,chowAllIndACI];
%%
% Run t-tests
[~,p(1),~,stats{1}] = ttest2(mean(baseBaseDiag,1),popBaseBaseA);
[~,p(2),~,stats{2}] = ttest2(mean(baseDep24Diag,1),popBaseDep24A);
[~,p(3),~,stats{3}] = ttest2(mean(baseDep48Diag,1),popBaseDep48A);
[~,p(4),~,stats{4}] = ttest2(mean(baseChowDiag,1),popBaseChowA);
[~,p(5),~,stats{5}] = ttest2(mean(baseDiag,1),popBaseA);
[~,p(6),~,stats{6}] = ttest2(mean(dep24Diag,1),popDep24A);
[~,p(7),~,stats{7}] = ttest2(mean(dep48Diag,1),popDep48A);
[~,p(8),~,stats{8}] = ttest2(mean(chowDiag,1),popChowA);

% [~,p(9),~,stats{9}] = ttest2(popBaseBaseA,popBaseA);
% [~,p(10),~,stats{10}] = ttest2(mean(baseBaseDiag,1),mean(baseDiag,1));
% [~,p(11),~,stats{11}] = ttest2(popBaseDep24A,popDep24A);
% [~,p(12),~,stats{12}] = ttest2(mean(baseDep24Diag,1),mean(dep24Diag,1));
% [~,p(13),~,stats{13}] = ttest2(popBaseDep48A,popDep48A);
% [~,p(14),~,stats{14}] = ttest2(mean(baseDep48Diag,1),mean(dep48Diag,1));
% [~,p(15),~,stats{15}] = ttest2(popBaseChowA,popChowA);
% [~,p(16),~,stats{16}] = ttest2(mean(baseChowDiag,1),mean(chowDiag,1));
% Apply correction
p1 = p(1:4).*8;
p2 = p(5:8).*8;
% Set up comparison groups
comp = {[1,2],[3,4],[5,6],[7,8]};
% Remove p-values and 'comp' groups if p>0.05
comp1 = comp(p1<0.05);
p1 = p1(p1<0.05);
comp2 = comp(p2<0.05);
p2 = p2(p2<0.05);
data1 = {popBaseBaseA',mean(baseBaseDiag)',popBaseDep24A',...
    mean(baseDep24Diag)',popBaseDep48A',mean(baseDep48Diag)',...
    popBaseChowA',mean(baseChowDiag)'};
data2 = {popBaseA',mean(baseDiag)',popDep24A',mean(dep24Diag)',...
    popDep48A',mean(dep48Diag)',popChowA',mean(chowDiag)'};
%% Plot
figure
subplot(1,2,1)
hold on
% Plot raw data with jitter
% plotSpread(data1,'distributionmarker','.','distributioncolors','k')
% Plot mean and 95% CI
scatterErr(1:2:8,baseM(1:2:8),baseCI(1:2:8),0,'mark','o','col','k')
scatterErr(2:2:8,baseM(2:2:8),baseCI(2:2:8),0,'mark','o','col',...
    [0.5 0.5 0.5])
% Plot significance bars
sigstar(comp1,p1,1);
xlim([0.5 8.5])
ylabel('AUC')
set(gca,'XTick',1:8,'XTickLabel',{'popBase','indBase','popDep24',...
    'indDep24','popDep48','indDep48','popChow','indChow'})
xtickangle(90)
title('Baseline')

subplot(1,2,2)
% Plot raw data with jitter
% plotSpread(data2,'distributionmarker','.','distributioncolors','k')
% Plot mean and 95% CI
scatterErr(1:2:8,allM(1:2:8),allCI(1:2:8),0,'mark','o','col','k')
scatterErr(2:2:8,allM(2:2:8),allCI(2:2:8),0,'mark','o','col',[0.5 0.5 0.5])

% Plot significance bars
sigstar(comp2,p2,1);
xlim([0.5 8.5])
ylabel('AUC')
set(gca,'XTick',1:8,'XTickLabel',{'popBase','indBase','popDep24',...
    'indDep24','popDep48','indDep48','popChow','indChow'})
xtickangle(90)
title('All')

%% Leave N Out
for ii = 1:11
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
        'finalNew\LNO\',num2str(ii),'.mat'])
    a{ii} = lno.a;
end
aM = cell2mat(cellfun(@(x) mean(mean(x,2)),a,'UniformOutput',0));
aCI = cell2mat(cellfun(@(x) conf(mean(x,2)',0.95,'tail',2),a,...
    'UniformOutput',0));
% aS = cell2mat(cellfun(@(x) std(x),a,'UniformOutput',0));
% Get data for concatLog models
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLog.mat'],'ccLogA','ccLogAM','ccLogACI')
figure
scatterErr(1:12,[ccLogAM,aM],[ccLogACI,aCI],0)
% scatterErr(1:12,[ccLogAM,aM],[std(ccLogA),aS],1)
set(gca,'XTick',1:12,'XTickLabel',0:11)
xlabel('Number of Animals Left Out')
ylabel('AUC')
%% Load uni, dy, and tri
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\'...
    'figure7\'])
[uniBase,uniDep24,uniDep48,uniChow] = deal(zeros(8,58,20));
[dyBase,dyDep24,dyDep48,dyChow] = deal(zeros(8,1653,20));
[triBase,triDep24,triDep48,triChow] = deal(zeros(8,30856,20));
for ii = 1:60
    load([num2str(ii),'.mat'])
    c = ceil(ii/20);
    num = rem(ii,20);
    if num == 0
        num = 20;
    end
    if c == 1
        uniBase(:,:,num) = baseA;
        uniDep24(:,:,num) = dep24A;
        uniDep48(:,:,num) = dep48A;
        uniChow(:,:,num) = chowA;
    elseif c == 2
       dyBase(:,:,num) = baseA;
       dyDep24(:,:,num) = dep24A;
       dyDep48(:,:,num) = dep48A;
       dyChow(:,:,num) = chowA;
    elseif c == 3
        triBase(:,:,num) = baseA;
        triDep24(:,:,num) = dep24A;
        triDep48(:,:,num) = dep48A;
        triChow(:,:,num) = chowA;
    end
end
%% Load log (58) and lasso
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\'...
    'allEachLasso\'])
[lassoBase,lassoDep24,lassoDep48,lassoChow] = deal(zeros(8,20));
for ii = 1:20
   load([num2str(ii),'.mat'],'indBase','indDep24','indDep48','indChow',...
       'indAll') 
   lassoBase(:,ii) = diag(indBase.A);
   lassoDep24(:,ii) = diag(indDep24.A);
   lassoDep48(:,ii) = diag(indDep48.A);
   lassoChow(:,ii) = diag(indChow.A);
end
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\'...
    'allEach\'])
[logBase,logDep24,logDep48,logChow] = deal(zeros(8,20));
for ii = 1:20
   load([num2str(ii),'.mat'],'indBase','indDep24','indDep48','indChow',...
       'indAll') 
   logBase(:,ii) = diag(indBase.A);
   logDep24(:,ii) = diag(indDep24.A);
   logDep48(:,ii) = diag(indDep48.A);
   logChow(:,ii) = diag(indChow.A);
end
%% Get means and confidence intervals
uniMeans = cellfun(@(x) mean(mean(x,3),1),{uniBase,uniDep24,uniDep48,...
    uniChow},'UniformOutput',0);
uniCI = cellfun(@(x) conf(mean(x,3)',0.95,'tail',2),{uniBase,uniDep24,...
    uniDep48,uniChow},'UniformOutput',0);
dyMeans = cellfun(@(x) mean(mean(x,3),1),{dyBase,dyDep24,dyDep48,dyChow}...
    ,'UniformOutput',0);
dyCI = cellfun(@(x) conf(mean(x,3)',0.95,'tail',2),{dyBase,dyDep24,...
    dyDep48,dyChow},'UniformOutput',0);
triMeans = cellfun(@(x) mean(mean(x,3),1),{triBase,triDep24,triDep48,...
    triChow},'UniformOutput',0);
triCI = cellfun(@(x) conf(mean(x,3)',0.95,'tail',2),{triBase,triDep24,...
    triDep48,triChow},'UniformOutput',0);
lassoMeans = cellfun(@(x) mean(mean(x,2)),{lassoBase,lassoDep24,...
    lassoDep48,lassoChow},'UniformOutput',0);
lassoCI = cellfun(@(x) conf(mean(x,2)',0.95,'tail',2),{lassoBase,...
    lassoDep24,lassoDep48,lassoChow},'UniformOutput',0);
logMeans = cellfun(@(x) mean(mean(x,2)),{logBase,logDep24,logDep48,...
    logChow},'UniformOutput',0);
logCI = cellfun(@(x) conf(mean(x,2)',0.95,'tail',2),{logBase,logDep24,...
    logDep48,logChow},'UniformOutput',0);
% Get top performing features according to baseline performance
[~,uniInd] = max(uniMeans{1,1});
[~,dyInd] = max(dyMeans{1,1});
[~,triInd] = max(triMeans{1,1});
topA = [uniMeans{1,1}(uniInd),uniMeans{1,2}(uniInd),...
    uniMeans{1,3}(uniInd),uniMeans{1,4}(uniInd);dyMeans{1,1}(dyInd),...
    dyMeans{1,2}(dyInd),dyMeans{1,3}(dyInd),dyMeans{1,4}(dyInd);...
    triMeans{1,1}(triInd),triMeans{1,2}(triInd),triMeans{1,3}(triInd),...
    triMeans{1,4}(triInd);cell2mat(lassoMeans);cell2mat(logMeans)];
topACI = [uniCI{1,1}(uniInd),uniCI{1,2}(uniInd),uniCI{1,3}(uniInd),...
    uniCI{1,4}(uniInd);dyCI{1,1}(dyInd),dyCI{1,2}(dyInd),...
    dyCI{1,3}(dyInd),dyCI{1,4}(dyInd);triCI{1,1}(triInd),...
    triCI{1,2}(triInd),triCI{1,3}(triInd),triCI{1,4}(triInd);...
    cell2mat(lassoCI);cell2mat(logCI)];
% Plot each state independently
main = {'Base','Dep24','Dep48','Chow'};
figure
for ii = 1:4
    subplot(2,2,ii)
    scatterErr(1:5,topA(:,ii),topACI(:,ii),0)
    title(main{ii})
    ylim([0.4 1])
    set(gca,'XTick',1:5,'XTickLabel',{'Monad','Dyad','Triad','Lasso',...
        'Log'})
    xtickangle(45)
end
% Plot all
scatterErr(1:20,reshape(topA',1,20),reshape(topACI',1,20),1)
set(gca,'XTick',1:20,'XTickLabel',repmat({'Base','Dep24','Dep48','Chow'}...
    ,1,4))
xtickangle(90)
ylabel('AUC')
% Run t-tests
[~,p(1),~,stats{1}] = ttest2(mean(logBase,2),mean(uniBase(:,uniInd,:),3));
[~,p(2),~,stats{2}] = ttest2(mean(logDep24,2),mean(uniDep24(:,uniInd,:),...
    3));
[~,p(3),~,stats{3}] = ttest2(mean(logDep48,2),mean(uniDep48(:,uniInd,:),...
    3));
[~,p(4),~,stats{4}] = ttest2(mean(logChow,2),mean(uniChow(:,uniInd,:),3));

[~,p(5),~,stats{5}] = ttest2(mean(logBase,2),mean(dyBase(:,dyInd,:),3));
[~,p(6),~,stats{6}] = ttest2(mean(logDep24,2),mean(dyDep24(:,dyInd,:),3));
[~,p(7),~,stats{7}] = ttest2(mean(logDep48,2),mean(dyDep48(:,dyInd,:),3));
[~,p(8),~,stats{8}] = ttest2(mean(logChow,2),mean(dyChow(:,dyInd,:),3));

[~,p(9),~,stats{9}] = ttest2(mean(logBase,2),mean(triBase(:,triInd,:),3));
[~,p(10),~,stats{10}] = ttest2(mean(logDep24,2),...
    mean(triDep24(:,triInd,:),3));
[~,p(11),~,stats{11}] = ttest2(mean(logDep48,2),...
    mean(triDep48(:,triInd,:),3));
[~,p(12),~,stats{12}] = ttest2(mean(logChow,2),...
    mean(triChow(:,triInd,:),3));

[~,p(13),~,stats{13}] = ttest2(mean(logBase,2),mean(lassoBase,2));
[~,p(14),~,stats{14}] = ttest2(mean(logDep24,2),mean(lassoDep24,2));
[~,p(15),~,stats{15}] = ttest2(mean(logDep48,2),mean(lassoDep48,2));
[~,p(16),~,stats{16}] = ttest2(mean(logChow,2),mean(lassoChow,2));
% Correct p-values
p = p.*16;
% Groups
group = {[1,17],[2,18],[3,19],[4,20],[5,17],[6,18],[7,19],[8,20],[9,17],...
    [10,18],[11,19],[12,20],[13,17],[14,18],[15,19],[16,20]};
group = group(p<0.05);
p = p(p<0.05);
sigstar(group,p)
%%
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
    'finalNew\indvGen.mat'])
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
    'finalNew\baseline500Each6000All50-50.mat'],'pInds')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
for ii = 1:20
    for jj = 1:8
        trainX = indTrainX{ii,jj}(:,inds);
        trainX = trainX(:,[6,21]);
        trainY = indTrainY{ii,jj};
        dyad{ii,jj} = fitglm(trainX,trainY,'distribution','binomial');
        trainX = indTrainX{ii,jj}(:,inds);
        trainX = trainX(:,[3,12,54]);
        trainY = indTrainY{ii,jj};
        triad{ii,jj} = fitglm(trainX,trainY,'distribution','binomial');
    end
end
%% Extract betas from models
dyadBeta = cellfun(@(x) table2array(x.Coefficients(2:3,1)),dyad,'UniformOutput',0);
clear dyadSign
dyadSign = zeros(2,8,20);
for ii = 1:20
   dyadSign(:,:,ii) = cat(2,dyadBeta{ii,:}); 
end
dyadSign = dyadSign>0;
test = mean(dyadSign,3);