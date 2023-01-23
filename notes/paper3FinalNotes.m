%% Grab all water and alcohol data
% Cross compare feature lists
waterFeat = names({'ILL','CA1L','PL','SL','PR','CA1R','ILR','SR'},...
    {'d','t','a','b','lg','hg'});
alcFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
% Find overlap
inds = zeros(1,60);
for ii = 1:60
    ind = logicFind(alcFeat{ii},waterFeat,'==');
    if ~isempty(ind)
        inds(ii) = ind;
    end
end
% Doesn't find PR-SL coherence since in waterFeat it is SL-PR
inds(43:48) = 157:162;
% Determine overlap between alcFeat and waterAlcFeat (same features,
% different order)
waterAlcFeat = names({'PL','SL','PR','SR'},{'d','t','a','b','lg','hg'});
inds2 = zeros(1,60);
for ii = 1:60
    ind = logicFind(alcFeat{ii},waterAlcFeat,'==');
    if ~isempty(ind)
        inds2(ii) = ind;
    end
end
% Doesn't find PR-SL coherence since in waterFeat it is SL-PR
inds2(43:48) = 43:48;
%__________________________________________________________________________
% Grab water only data
[wData,~,wFiles] = collateData(['D:\paper3\'...
    'waterProcessed\'],{'.mat'},{'pow','coh'},'trl','rel');
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),wData{1,1}(:,1),'uniformoutput',0));
nFiles = size(wData{1},1);
[animal,day] = deal(cell(nFiles,1));
for ii = 1:nFiles
    parts = strsplit(wFiles{1}{ii},'_');
    animal{ii,1} = repmat({['AH',parts{2}]},n(ii,1),1);
    day{ii,1} = repmat({parts{3}},n(ii,1),1); %#ok
end
% Concatenate and z-score water data
% allWater = zscore(cat(1,wData{1}{:,1}));
% Concatenate
allWater = cat(1,wData{1}{:,1});
% Limit and tabulate water data
water = [array2table(allWater(:,inds)),cat(1,animal{:}),cat(1,day{:})];
% Add groups (0 = water)
water(:,63) = array2table(zeros(size(water,1),1));
water.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];

% Get notWater data
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),wData{1,1}(:,2),'uniformoutput',0));
nFiles = size(wData{1},1);
[animal,day] = deal(cell(nFiles,1));
for ii = 1:nFiles
    parts = strsplit(wFiles{1}{ii},'_');
    animal{ii,1} = repmat({['AH',parts{2}]},n(ii,1),1);
    day{ii,1} = repmat({parts{3}},n(ii,1),1); %#ok
end
% Concatenate and z-score water data
% allWater = zscore(cat(1,wData{1}{:,1}));
% Concatenate
allNotWater = cat(1,wData{1}{:,2});
% Limit and tabulate water data
notWater = [array2table(allNotWater(:,inds)),cat(1,animal{:}),...
    cat(1,day{:})];
% Add groups (0 = water)
notWater(:,63) = array2table(zeros(size(notWater,1),1));
notWater.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];
%__________________________________________________________________________
% Grab alcohol only data
[aData,~,aFiles] = collateData('D:\paper3\drinkNot\',{'.mat'},...
    {'pow','coh'},'trl','rel');
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),aData{1,1}(:,1),'uniformoutput',0));
nFiles = size(aData{1},1);
[animal,day] = deal(cell(nFiles,1));
for ii = 1:nFiles
    parts = strsplit(aFiles{1}{ii},'_');
    animal{ii,1} = repmat({parts{1}},n(ii,1),1); %#ok
    day{ii,1} = repmat({parts{2}},n(ii,1),1); %#ok
end
% Concatenate
allAlc = cat(1,aData{1}{:,1});
% Limit and tabulate alcohol data
alc = [array2table(allAlc(:,1:60)),cat(1,animal{:}),cat(1,day{:})];
% Add groups (1 = alcohol)
alc(:,63) = array2table(ones(size(alc,1),1));
alc.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];
% Grab not alcohol data
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),aData{1,1}(:,2),'uniformoutput',0));
[animal,day] = deal(cell(nFiles,1));
for ii = 1:size(aData{1},1)
    parts = strsplit(aFiles{1}{ii},'_');
    animal{ii,1} = repmat({parts{1}},n(ii,1),1); %#ok
    day{ii,1} = repmat({parts{2}},n(ii,1),1); %#ok
end
% Concatenate
allNotAlc = cat(1,aData{1}{:,2});
% Limit and tabulate alcohol data
notAlc = [array2table(allNotAlc(:,1:60)),cat(1,animal{:}),cat(1,day{:})];
% Add groups (1 = alcohol)
notAlc(:,63) = array2table(zeros(size(notAlc,1),1));
notAlc.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];
%__________________________________________________________________________
% Grab drinkNot data (taken from waterAlcohol cohort but set to look at
% data corresponding to neither behavior)
[dataNot,~,files] = collateData('D:\paper3\waterAlcohol\notDrink\',...
    {'.mat'},{'pow','coh'},'trl','rel');
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),dataNot{1,1}(:,1),'uniformoutput',0));
nFiles = size(dataNot{1},1);
[animal,day] = deal(cell(nFiles,1));
for ii = 1:nFiles
    parts = strsplit(files{1}{ii},'_');
    animal{ii,1} = repmat({parts{1}},n(ii,1),1);
    day{ii,1} = repmat({parts{3}},n(ii,1),1); %#ok
end
% Concatenate not alcohol data
allDataNot = cat(1,dataNot{1}{:,1});
% Limit and tabulate alcohol data
notDrink = [array2table(allDataNot(:,1:60)),cat(1,animal{:}),cat(1,day{:})];
% Add groups (0 = not)
notDrink(:,63) = array2table(zeros(size(notDrink,1),1));
notDrink.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];
% Combine into notD
notD = [notAlc;notDrink];
%__________________________________________________________________________
% Grab water/alcohol data
[awData,~,awFiles] = collateData('D:\paper3\waterAlcohol\processed\',...
    {'.mat'},{'pow','coh'},'trl'...
    ,'rel');
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),awData{1,1},'uniformoutput',0));
nFiles = size(awData{1},1);
[animal,day] = deal(cell(nFiles,1));
for ii = 1:nFiles
    parts = strsplit(awFiles{1}{ii},'_');
    for jj = 1:2
        animal{ii,jj} = repmat({parts{1}},n(ii,jj),1);
        day{ii,jj} = repmat({parts{3}},n(ii,jj),1); %#ok
    end
end
% Concatenate
allWatAlc = [cat(1,awData{1}{:,1});cat(1,awData{1}{:,2})];
% Limit and tabulate water/alcohol data
watAlc = [array2table(allWatAlc(:,inds2)),[cat(1,animal{:,1});...
    cat(1,animal{:,2})],[cat(1,day{:,1});cat(1,day{:,2})]];
% Add groups (1 = alcohol; 0 = water)
watAlc(:,63) = array2table([zeros(sum(n(:,1)),1);ones(sum(n(:,2)),1)]);
watAlc.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];
%__________________________________________________________________________
% Group data together (animal x condition [alcohol,water,neither])
ids = unique([unique(notD.ID);unique(water.ID);unique(alc.ID);...
    unique(watAlc.ID)]);
nIDs = numel(ids);
% Preallocate
samps = zeros(nIDs,3);
for ii = 1:nIDs
    % Alcohol drinking
    samps(ii,1) = sum(strcmp(watAlc.ID,ids{ii}) & ...
        watAlc.group == 1) + sum(strcmp(alc.ID,ids{ii}) & ...
        alc.group == 1);
    % Water drinking
    samps(ii,2) = sum(strcmp(watAlc.ID,ids{ii}) & ...
        watAlc.group == 0) + sum(strcmp(water.ID,ids{ii}) ...
        & water.group == 1);
    % Neither
    samps(ii,3) = sum(strcmp(notD.ID,ids{ii}) & notD.group == 0);
end
theseIDs = ids(sum(samps>0,2)==3);
theseSamps = samps(sum(samps>0,2)==3,:);
% Grab and combine data
nTheseIDs = numel(theseIDs);
% Preallocate
[thisAlc,thisWat,thisNeither] = deal(cell(nTheseIDs,1));
for ii = 1:nTheseIDs
    thisAlc{ii,1} = [watAlc(strcmp(watAlc.ID,theseIDs{ii}) & ...
        watAlc.group == 1,:); alc(strcmp(alc.ID,theseIDs{ii}) & ...
        alc.group == 1,:)];
    thisWat{ii,1} = [watAlc(strcmp(watAlc.ID,theseIDs{ii}) & ...
        watAlc.group == 0,:); water(strcmp(water.ID,theseIDs{ii}) ...
        & water.group == 1,:)];
    thisNeither{ii,1} = notD(strcmp(notD.ID,theseIDs{ii}) & ...
        notD.group == 0,:);
end
% Combine
allData = [thisAlc,thisWat,thisNeither]; %#ok
save('D:\paper3\analyzed\final\waterAlcoholNotData.mat','alc','allData',...
    'notD','watAlc','water','notWater')
%% Plot drinking amounts
load('D:\paper3\analyzed\final\drinkAmounts.mat')
figure
plotSpread(other)
hold on
plotSpread(these,'distributionColors','r')
ylabel('EtOH consumed (g/kg)')
%% and feeding amounts
load('D:\paper2\4conditionBingeSize.mat')
figure
plotSpread(bingeCal(:,2))
ylabel('sweet-fat food consumed (kCal)')
%% Feeding data
[data,~,~] = collateData('F:\paper2\preBinge_240sec_redo\',...
    {'dep24';'chow'},{'pow','coh'},'trl','rel');
both = [1:3,5:8,10,11];
allFeedData(:,1) = data{1,1}(both,242);
allFeedData(:,2) = data{1,2}(:,242);
bingeSamps = cellfun(@(x) size(x,1),allFeedData);
for ii = 1:9
    allFeedData{ii,3} = [data{1}{both(ii),241};data{2}{ii,241}];
end
% Set up all y vectors
allY = cell(size(allFeedData));
for ii = 1:size(allFeedData,1)
    for jj = 1:size(allFeedData,2)
        allY{ii,jj} = ones(bingeSamps(ii,jj),1);
        % Set neither column and chow to zero by subtracting one
        if jj ~= 1
           allY{ii,jj} = allY{ii,jj} - 1; 
        end
    end
end
notFeed = [data{1}(both,241),data{2}(:,241)];
save('F:\paper3\analyzed\final\24sweetChowNotDataNew.mat',...
    'allFeedData','allY','bingeSamps','notFeed')
%% Generalized model
% N.B.: run this with the data from 1feedWeight or 1feedWeightLOO
% Data from discovery (runBingeAlcoholOther)
steps = 1;
nSteps = numel(steps);
nGen = 232;
nF = 150;
nD = 100;
% Prealloate
[genAcc,fAcc,dAcc,genAccR] = deal(zeros(81,nSteps));
[genX,genY,genXR,genYR] = deal(zeros(nSteps,81,nGen));
[fX,fY] = deal(zeros(nSteps,81,nF));
[dX,dY] = deal(zeros(nSteps,81,nD));
c = 1;
for jj = steps
    cd(['D:\paper3\analyzed\final\genModel\',num2str(jj),'feedWeightLOO'])
    for ii = 1:81
        % Load data and calculate ROCs
        load(['feedDrinkOtherLOO',num2str(ii),'.mat'],'acc','accR',...
            'drinkAcc','feedAcc','hist','histR','feedRocX','feedRocY',...
            'drinkRocX','drinkRocY')
        [x,y,~,genAcc(ii,c)] = perfcurve(hist.cfg.naive.testY,...
            acc{1}.pred,1);
        genX(c,ii,:) = interp1(linspace(0,1,numel(x)),x,...
            linspace(0,1,nGen));
        genY(c,ii,:) = interp1(linspace(0,1,numel(y)),y,...
            linspace(0,1,nGen));
        [xR,yR,~,genAccR(ii,c)] = perfcurve(histR.cfg.naive.testY,...
            accR{1}.pred,1);
        genXR(c,ii,:) = interp1(linspace(0,1,numel(xR)),xR,...
            linspace(0,1,nGen));
        genYR(c,ii,:) = interp1(linspace(0,1,numel(yR)),yR,...
            linspace(0,1,nGen));
        fAcc(ii,c) = feedAcc;
        fX(c,ii,:) = interp1(linspace(0,1,numel(feedRocX)),feedRocX,...
            linspace(0,1,nF));
        fY(c,ii,:) = interp1(linspace(0,1,numel(feedRocY)),feedRocY,...
            linspace(0,1,nF));
        dAcc(ii,c) = drinkAcc;
        dX(c,ii,:) = interp1(linspace(0,1,numel(drinkRocX)),drinkRocX,...
            linspace(0,1,nD));
        dY(c,ii,:) = interp1(linspace(0,1,numel(drinkRocY)),drinkRocY,...
            linspace(0,1,nD));
    end
    c = c+1;
end
% Plot results of gen model ()
this = 1;
figure
hold on
plot(mean(squeeze(genX(this,:,:)),1),mean(squeeze(genY(this,:,:)),1))
plot(mean(squeeze(genXR(this,:,:)),1),mean(squeeze(genYR(this,:,:)),1))
plot(mean(squeeze(dX(this,:,:)),1),mean(squeeze(dY(this,:,:)),1))
plot(mean(squeeze(fX(this,:,:)),1),mean(squeeze(fY(this,:,:)),1))
legend({['Gen: ',num2str(round(mean(genAcc(:,this)),2)),'\pm',...
    num2str(round(conf(genAcc(:,this)',0.95),2))],...
    ['GenP: ',num2str(round(mean(genAccR(:,this)),2)),'\pm',...
    num2str(round(conf(genAccR(:,this)',0.95),3))],...
    ['Drink: ',num2str(round(mean(dAcc(:,this)),2)),'\pm',...
    num2str(round(conf(dAcc(:,this)',0.95),2))],...
    ['Feed: ',num2str(round(mean(fAcc(:,this)),2)),'\pm',...
    num2str(round(conf(fAcc(:,this)',0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title(['generalized model: ',num2str(steps(this)),'feedWeight'])
% % Plot AUC across steps
% figure
% hold on
% plot(steps,[mean(genAcc,1);mean(genAccR,1);mean(fAcc,1);mean(dAcc,1)]',...
%     '-o')
% legend({'Gen','GenP','Feed','Drink'})
% xlabel('feed data weight relative to drink weight = 1')
% ylabel('AUC')
% title('Feed weight lasso')
% 
% figure
% hold on
% shadedErrorBar(steps,mean(fAcc),std(fAcc),[],1)
% shadedErrorBar(steps,mean(dAcc),std(dAcc),[],1)
%% 80:20 feed vs. other
cd D:\paper3\analyzed\final\feedOtherNew\
for ii = 1:100
   load(['feedOther',num2str(ii),'.mat'])
   [thisX,thisY,~,feedA(ii)] = perfcurve(...
       hist.cfg.naive.testY,acc{1}.pred{1},1);
   [thisXR,thisYR,~,feedAR(ii)] = perfcurve(...
       histR.cfg.naive.testY,accR{1}.pred{1},1);
   feedX(ii,:) = interp1(linspace(0,1,numel(thisX)),...
       thisX,linspace(0,1,81));
   feedY(ii,:) = interp1(linspace(0,1,numel(thisY)),...
       thisY,linspace(0,1,81));
   feedXR(ii,:) = interp1(linspace(0,1,numel(thisXR)),...
       thisXR,linspace(0,1,81));
   feedYR(ii,:) = interp1(linspace(0,1,numel(thisYR)),...
       thisYR,linspace(0,1,81));
   % Test model on sweet-fat vs. chow and neither independently
   %        [thisX,thisY,~,chowA(ii,jj)] = perfcurve(...
   %            hist{jj}.cfg.naive.testY(1:indSplit(jj,2)),...
   %            acc{jj}{1}.pred(1:indSplit(jj,2)),1);
   [thisX,thisY,~,chowA(ii)] = perfcurve(...
      hist.cfg.naive.testY(testDataSet==1|testDataSet==2),...
       acc{1}.pred{1}(testDataSet==1|testDataSet==2),1);
   chowX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
       linspace(0,1,81));
   chowY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
       linspace(0,1,81));
   %        [thisX,thisY,~,neitherFA(ii,jj)] = perfcurve(...
   %            hist{jj}.cfg.naive.testY([1:indSplit(jj,1),...
   %            indSplit(jj,2)+1:end]),acc{jj}{1}.pred([1:indSplit(jj,1),...
   %            indSplit(jj,2)+1:end]),1);
   [thisX,thisY,~,neitherFA(ii)] = perfcurve(...
       hist.cfg.naive.testY(testDataSet==1|testDataSet==3)...
       ,acc{1}.pred{1}(testDataSet==1|testDataSet==3),1);
   neitherFX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
       linspace(0,1,81));
   neitherFY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
       linspace(0,1,81));
end
%%
chowP = (sum(feedAR>mean(chowA))+1)/(numel(feedAR)+1);
neitherP = (sum(feedAR>mean(neitherFA))+1)/(numel(feedAR)+1);
figure
hold on
plot(mean(feedX,1),mean(feedY,1))
plot(mean(chowX,1),mean(chowY,1))
plot(mean(neitherFX,1),mean(neitherFY,1))
plot(mean(feedXR,1),mean(feedXR,1))
legend({['feed: ',num2str(round(mean(feedA),2)),'\pm',...
    num2str(round(conf(feedA,0.95),2))],...
    ['sf vs. chow: ',num2str(round(mean(chowA),2)),'\pm',...
    num2str(round(conf(chowA,0.95),2)),' ',num2str(chowP)],...
    ['sf vs. neither: ',num2str(round(mean(neitherFA),2)),'\pm',...
    num2str(round(conf(neitherFA,0.95),2)),' ',num2str(neitherP)],...
    ['permuted: ',num2str(round(mean(feedAR),2)),'\pm',...
    num2str(round(conf(feedAR,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('feed: 80:20')
%% 80:20 drink vs. other
cd D:\paper3\analyzed\final\drinkOtherNew\
for ii = 1:100
   load(['drinkOther',num2str(ii),'.mat'])
   [thisX,thisY,~,drinkA(ii)] = perfcurve(...
       hist.cfg.naive.testY,acc{1}.pred{1},1);
   [thisXR,thisYR,~,drinkAR(ii)] = perfcurve(...
       histR.cfg.naive.testY,accR{1}.pred{1},1);
   drinkX(ii,:) = interp1(linspace(0,1,numel(thisX)),...
       thisX,linspace(0,1,81));
   drinkY(ii,:) = interp1(linspace(0,1,numel(thisY)),...
       thisY,linspace(0,1,81));
   drinkXR(ii,:) = interp1(linspace(0,1,numel(thisXR)),...
       thisXR,linspace(0,1,81));
   drinkYR(ii,:) = interp1(linspace(0,1,numel(thisYR)),...
       thisYR,linspace(0,1,81));
   % Test on water
   [thisX,thisY,~,waterA(ii)] = perfcurve(...
      hist.cfg.naive.testY(testDataSet==1|testDataSet==2),...
       acc{1}.pred{1}(testDataSet==1|testDataSet==2),1);
   waterX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
       linspace(0,1,81));
   waterY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
       linspace(0,1,81));
   % Test on neither
   [thisX,thisY,~,neitherDA(ii)] = perfcurve(...
       hist.cfg.naive.testY(testDataSet==1|testDataSet==3)...
       ,acc{1}.pred{1}(testDataSet==1|testDataSet==3),1);
   neitherDX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
       linspace(0,1,81));
   neitherDY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
       linspace(0,1,81));
end
%%
waterP = (sum(drinkAR>mean(waterA))+1)/(numel(drinkAR)+1);
neitherP = (sum(drinkAR>mean(neitherDA))+1)/(numel(drinkAR)+1);
figure
hold on
plot(mean(drinkX,1),mean(drinkY,1))
plot(mean(waterX,1),mean(waterY,1))
plot(mean(neitherDX,1),mean(neitherDY,1))
plot(mean(drinkXR,1),mean(drinkXR,1))
legend({['drink: ',num2str(round(mean(drinkA),2)),'\pm',...
    num2str(round(conf(drinkA,0.95),2))],...
    ['EtOH vs. water: ',num2str(round(mean(waterA),2)),'\pm',...
    num2str(round(conf(waterA,0.95),2)),' ',num2str(waterP)],...
    ['EtOH vs. neither: ',num2str(round(mean(neitherDA),2)),'\pm',...
    num2str(round(conf(neitherDA,0.95),2)),' ',num2str(neitherP)],...
    ['permuted: ',num2str(round(mean(drinkAR),2)),'\pm',...
    num2str(round(conf(drinkAR,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('drink: 80:20')
%%
cd D:\paper3\analyzed\final\drinkOtherLOONew\
% cd D:\paper3\analyzed\final\drinkOtherAll\
[x,y,xR,yR] = deal(cell(100,9));
% Set up indices of where EtOH samples and water samples end 
indSplit = [4,15;6,26;40,43;40,48;40,60;3,5;5,10;40,51;5,25];
for ii = 1:100
   load(['drinkOtherLOO',num2str(ii),'.mat'],'acc','accR','hist','histR')
%    load(['drinkOtherAll',num2str(ii),'.mat'],'acc','accR','hist','histR')
   for jj = 1:9
       [thisX,thisY,~,drinkA(ii,jj)] = perfcurve(...
           hist{jj}.cfg.naive.testY,acc{jj}{1}.pred{1},1);
       [thisXR,thisYR,~,drinkAR(ii,jj)] = perfcurve(...
           histR{jj}.cfg.naive.testY,accR{jj}{1}.pred{1},1);
       drinkX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),...
           thisX,linspace(0,1,81));
       drinkY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),...
           thisY,linspace(0,1,81));
       drinkXR(ii,jj,:) = interp1(linspace(0,1,numel(thisXR)),...
           thisXR,linspace(0,1,81));
       drinkYR(ii,jj,:) = interp1(linspace(0,1,numel(thisYR)),...
           thisYR,linspace(0,1,81));
       % Test model on EtOH vs. water and neither independently
       [thisX,thisY,~,waterA(ii,jj)] = perfcurve(...
           hist{jj}.cfg.naive.testY(1:indSplit(jj,2)),...
           acc{jj}{1}.pred{1}(1:indSplit(jj,2)),1);
       waterX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
           linspace(0,1,81));
       waterY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
           linspace(0,1,81));
       [thisX,thisY,~,neitherDA(ii,jj)] = perfcurve(...
           hist{jj}.cfg.naive.testY([1:indSplit(jj,1),...
           indSplit(jj,2)+1:end]),acc{jj}{1}.pred{1}([1:indSplit(jj,1),...
           indSplit(jj,2)+1:end]),1);
       neitherDX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
           linspace(0,1,81));
       neitherDY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
           linspace(0,1,81));
   end
end
%% LOO feed vs. other - full model
cd D:\paper3\analyzed\final\feedOtherLOO\
for ii = 1:100
   load(['feedOtherLOO',num2str(ii),'.mat'],'acc','accR','hist',...
       'histR','testDataSet')
   for jj = 1:9
       [thisX,thisY,~,feedA(ii,jj)] = perfcurve(...
           hist{jj}.cfg.naive.testY,acc{jj}{1}.pred{1},1);
       [thisXR,thisYR,~,feedAR(ii,jj)] = perfcurve(...
           histR{jj}.cfg.naive.testY,accR{jj}{1}.pred{1},1);
       feedX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),...
           thisX,linspace(0,1,81));
       feedY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),...
           thisY,linspace(0,1,81));
       feedXR(ii,jj,:) = interp1(linspace(0,1,numel(thisXR)),...
           thisXR,linspace(0,1,81));
       feedYR(ii,jj,:) = interp1(linspace(0,1,numel(thisYR)),...
           thisYR,linspace(0,1,81));
       % Test model on sweet-fat vs. chow and neither independently
%        [thisX,thisY,~,chowA(ii,jj)] = perfcurve(...
%            hist{jj}.cfg.naive.testY(1:indSplit(jj,2)),...
%            acc{jj}{1}.pred(1:indSplit(jj,2)),1);
       [thisX,thisY,~,chowA(ii,jj)] = perfcurve(...
           hist{jj}.cfg.naive.testY(testDataSet{jj}==1|testDataSet{jj}==2),...
           acc{jj}{1}.pred{1}(testDataSet{jj}==1|testDataSet{jj}==2),1);
       chowX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
           linspace(0,1,81));
       chowY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
           linspace(0,1,81));
%        [thisX,thisY,~,neitherFA(ii,jj)] = perfcurve(...
%            hist{jj}.cfg.naive.testY([1:indSplit(jj,1),...
%            indSplit(jj,2)+1:end]),acc{jj}{1}.pred([1:indSplit(jj,1),...
%            indSplit(jj,2)+1:end]),1);
       [thisX,thisY,~,neitherFA(ii,jj)] = perfcurve(...
           hist{jj}.cfg.naive.testY(testDataSet{jj}==1|testDataSet{jj}==3)...
           ,acc{jj}{1}.pred{1}(testDataSet{jj}==1|testDataSet{jj}==3),1);
       neitherFX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
           linspace(0,1,81));
       neitherFY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
           linspace(0,1,81));
   end
end
%% prepare LOO data for plotting
% Reshape and mean ROC axes 
mFeedX = mean(reshape(feedX,900,81));
mFeedY = mean(reshape(feedY,900,81));
mFeedXR = mean(reshape(feedXR,900,81));
mFeedYR = mean(reshape(feedYR,900,81));
mChowX = mean(reshape(chowX,900,81));
mChowY = mean(reshape(chowY,900,81));
mNeitherFX = mean(reshape(neitherFX,900,81));
mNeitherFY = mean(reshape(neitherFY,900,81));

mDrinkX = mean(reshape(drinkX,900,81));
mDrinkY = mean(reshape(drinkY,900,81));
mDrinkXR = mean(reshape(drinkXR,900,81));
mDrinkYR = mean(reshape(drinkYR,900,81));
mWaterX = mean(reshape(waterX,900,81));
mWaterY = mean(reshape(waterY,900,81));
mNeitherDX = mean(reshape(neitherDX,900,81));
mNeitherDY = mean(reshape(neitherDY,900,81));
% Split drink test by sex
mInd = [4:5,8:9]; fInd = [1:3,6:7];
mDrinkMX = mean(reshape(drinkX(:,mInd,:),400,81));
mDrinkMY = mean(reshape(drinkY(:,mInd,:),400,81));
mDrinkFX = mean(reshape(drinkX(:,fInd,:),500,81));
mDrinkFY = mean(reshape(drinkY(:,fInd,:),500,81));
% Reshape accuracy matrices to vectors
fA = reshape(feedA,1,900);
fAR = reshape(feedAR,1,900);
cA = reshape(chowA,1,900);
fnA = reshape(neitherFA,1,900);
dA = reshape(drinkA,1,900);
dAR = reshape(drinkAR,1,900);
wA = reshape(waterA,1,900);
dnA = reshape(neitherDA,1,900);
% Split drink test by sex
dFA = reshape(drinkA(:,fInd),1,500);
dMA = reshape(drinkA(:,mInd),1,400);
% save('F:\paper3\analyzed\final\feedAlcoholOtherFullModelsLOO2.mat',...
%     'mFeedX','mFeedY','mFeedXR','mFeedYR','fA','fAR','mDrinkX',...
%     'mDrinkY','mDrinkXR','mDrinkYR','dA','dAR')
% save('F:\paper3\analyzed\final\feedAlcoholOtherFullModels.mat',...
%     'mFeedX','mFeedY','mFeedXR','mFeedYR','fA','fAR','mDrinkX',...
%     'mDrinkY','mDrinkXR','mDrinkYR','dA','dAR')
%% Plot ROCs
figure
hold on
plot(mFeedX,mFeedY)
plot(mFeedXR,mFeedYR)
plot(mDrinkX,mDrinkY)
plot(mDrinkXR,mDrinkYR)
legend({['feed: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],...
    ['feed perm: ',num2str(round(mean(fAR),2)),'\pm',...
    num2str(round(conf(fAR,0.95),2))],...
    ['drink: ',num2str(round(mean(dA),2)),'\pm',...
    num2str(round(conf(dA,0.95),2))],...
    ['drink perm: ',num2str(round(mean(dAR),2)),'\pm',...
    num2str(round(conf(dAR,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('behavior vs. not+other: all features')
% Plot auc histograms
doubleHist(fA,fAR)
legend({['Real: \mu = ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],...
    ['Permuted: \mu = ',num2str(round(mean(fAR),2)),'\pm',...
    num2str(round(conf(fAR,0.95),2))]},...
    'Location','nw');
title('Feed vs. not')

doubleHist(dA,dAR)
legend({['Real: \mu = ',num2str(round(mean(dA),2)),'\pm',...
    num2str(round(conf(dA,0.95),2))],...
    ['Permuted: \mu = ',num2str(round(mean(dAR),2)),'\pm',...
    num2str(round(conf(dAR,0.95),2))]},...
    'Location','nw');
title('Drink vs. not')
%% Split and compare test sets from above
chowP = (sum(fAR>mean(cA))+1)/(numel(fAR)+1);
neitherP = (sum(fAR>mean(fnA))+1)/(numel(fAR)+1);
figure
hold on
plot(mFeedX,mFeedY)
plot(mChowX,mChowY)
plot(mNeitherFX,mNeitherFY)
plot(mFeedXR,mFeedYR)
legend({['feed: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],...
    ['sf vs. chow: ',num2str(round(mean(cA),2)),'\pm',...
    num2str(round(conf(cA,0.95),2))],...
    ['sf vs. neither: ',num2str(round(mean(fnA),2)),'\pm',...
    num2str(round(conf(fnA,0.95),2))],...
    ['Permuted: ',num2str(round(mean(fAR),2)),'\pm',...
    num2str(round(conf(fAR,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('feed: LOO')

waterP = (sum(dAR>mean(wA))+1)/(numel(dAR)+1);
neitherP = (sum(dAR>mean(dnA))+1)/(numel(dAR)+1);
figure
hold on
plot(mDrinkX,mDrinkY)
plot(mWaterX,mWaterY)
plot(mNeitherDX,mNeitherDY)
plot(mDrinkXR,mDrinkYR)
legend({['drink: ',num2str(round(mean(dA),2)),'\pm',...
    num2str(round(conf(dA,0.95),2))],...
    ['EtOH vs. water: ',num2str(round(mean(wA),2)),'\pm',...
    num2str(round(conf(wA,0.95),2)),' ',num2str(waterP)],...
    ['EtOH vs. neither: ',num2str(round(mean(dnA),2)),'\pm',...
    num2str(round(conf(dnA,0.95),2)),' ',num2str(neitherP)],...
    ['drink perm: ',num2str(round(mean(dAR),2)),'\pm',...
    num2str(round(conf(dAR,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('drink: LOO')
%% LOO feed vs. other - subset features
% N.B.: run with or without LOO
% Feeding
cd D:\paper3\analyzed\final\feedOtherLOO\
% cd D:\paper3\analyzed\final\feedOther\
for ii = 1:100
   load(['feedOtherLOO',num2str(ii),'.mat'],'acc','accR','hist','histR')
%    load(['feedOther',num2str(ii),'.mat'],'acc','accR','hist','histR')
   for jj = 1:9
       [thisX,thisY,~,feedA(ii,jj)] = perfcurve(...
           hist{jj}.cfg.naive.testY,acc{jj}{1}.pred,1);
       [thisXR,thisYR,~,feedAR(ii,jj)] = perfcurve(...
           histR{jj}.cfg.naive.testY,accR{jj}{1}.pred,1);
       feedX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),...
           thisX,linspace(0,1,81));
       feedY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),...
           thisY,linspace(0,1,81));
       feedXR(ii,jj,:) = interp1(linspace(0,1,numel(thisXR)),...
           thisXR,linspace(0,1,81));
       feedYR(ii,jj,:) = interp1(linspace(0,1,numel(thisYR)),...
           thisYR,linspace(0,1,81));
   end
   % Apply to drink data
   load(['D:\paper3\analyzed\final\drinkOtherLOO\drinkOther',num2str(ii),...
       '.mat'],'hist')
%    load(['D:\paper3\analyzed\final\drinkOtherLOO\drinkOtherLOO',num2str(ii),...
%        '.mat'],'hist')
   for jj = 1:9
       pred = cvglmnetPredict(acc{jj}{1}.mdl{1},...
           hist{jj}.cfg.naive.testX,'lambda_1se','response');
       [thisX,thisY,~,feedDrinkA(ii,jj)] = perfcurve(hist{jj}.cfg.naive.testY,...
           pred,1);
       feedDrinkX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
           linspace(0,1,81));
       feedDrinkY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
           linspace(0,1,81));
   end
end
% Drinking
cd D:\paper3\analyzed\final\drinkOtherLOO\
% cd D:\paper3\analyzed\final\drinkOther\
[x,y,xR,yR] = deal(cell(100,9));
for ii = 1:100
   load(['drinkOtherLOO',num2str(ii),'.mat'],'acc','accR','hist','histR')
%    load(['drinkOther',num2str(ii),'.mat'],'acc','accR','hist','histR')
   for jj = 1:9
       [thisX,thisY,~,drinkA(ii,jj)] = perfcurve(...
           hist{jj}.cfg.naive.testY,acc{jj}{1}.pred,1);
       [thisXR,thisYR,~,drinkAR(ii,jj)] = perfcurve(...
           histR{jj}.cfg.naive.testY,accR{jj}{1}.pred,1);
       drinkX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),...
           thisX,linspace(0,1,81));
       drinkY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),...
           thisY,linspace(0,1,81));
       drinkXR(ii,jj,:) = interp1(linspace(0,1,numel(thisXR)),...
           thisXR,linspace(0,1,81));
       drinkYR(ii,jj,:) = interp1(linspace(0,1,numel(thisYR)),...
           thisYR,linspace(0,1,81));
   end
   % Apply to feed data
    load(['D:\paper3\analyzed\final\feedOtherLOO\feedOtherLOO',...
        num2str(ii),'.mat'],'hist')
%    load(['D:\paper3\analyzed\final\feedOther\feedOther',num2str(ii),...
%        '.mat'],'hist')
   for jj = 1:9
       pred = cvglmnetPredict(acc{jj}{1}.mdl{1},...
           hist{jj}.cfg.naive.testX,'lambda_1se','response');
       [thisX,thisY,~,drinkFeedA(ii,jj)] = perfcurve(...
           hist{jj}.cfg.naive.testY,pred,1);
       drinkFeedX(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
           linspace(0,1,81));
       drinkFeedY(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
           linspace(0,1,81));
   end
end
% Reshape and mean ROC axes 
mFeedX = mean(reshape(feedX,900,81));
mFeedY = mean(reshape(feedY,900,81));
mFeedXR = mean(reshape(feedXR,900,81));
mFeedYR = mean(reshape(feedYR,900,81));

mDrinkX = mean(reshape(drinkX,900,81));
mDrinkY = mean(reshape(drinkY,900,81));
mDrinkXR = mean(reshape(drinkXR,900,81));
mDrinkYR = mean(reshape(drinkYR,900,81));

mfdX = mean(reshape(feedDrinkX,900,81));
mfdY = mean(reshape(feedDrinkY,900,81));
mdfX = mean(reshape(drinkFeedX,900,81));
mdfY = mean(reshape(drinkFeedY,900,81));
% Reshape accuracy matrices to vectors
fA = reshape(feedA,1,900);
fAR = reshape(feedAR,1,900);
dA = reshape(drinkA,1,900);
dAR = reshape(drinkAR,1,900);
dfA = reshape(drinkFeedA,1,900);
fdA = reshape(feedDrinkA,1,900);
% save('F:\paper3\analyzed\final\feedAlcoholOtherSubsetModelsLOO.mat',...
%     'mFeedX','mFeedY','mFeedXR','mFeedYR','fA','fAR','mDrinkX',...
%     'mDrinkY','mDrinkXR','mDrinkYR','dA','dAR','mfdX','mfdY','mdfX',...
%      'mdfY','dfA','fdA')
% save('F:\paper3\analyzed\final\feedAlcoholOtherSubsetModels.mat',...
%     'mFeedX','mFeedY','mFeedXR','mFeedYR','fA','fAR','mDrinkX',...
%     'mDrinkY','mDrinkXR','mDrinkYR','dA','dAR','mfdX','mfdY','mdfX',...
%      'mdfY','dfA','fdA')
%% Plot ROCs
% N.B.: run this with the data from feedAlcoholOtherSubsetModels or
% feedAlcoholOtherSubsetModelsLOO
load('D:\paper3\analyzed\final\feedAlcoholOtherSubsetModels.mat')
% load('D:\paper3\analyzed\final\feedAlcoholOtherSubsetModelsLOO.mat')
% Within models
% ROCs
figure
hold on
plot(mFeedX,mFeedY)
plot(mFeedXR,mFeedYR)
plot(mDrinkX,mDrinkY)
plot(mDrinkXR,mDrinkYR)
legend({['feed: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],...
    ['feed perm: ',num2str(round(mean(fAR),2)),'\pm',...
    num2str(round(conf(fAR,0.95),2))],...
    ['drink: ',num2str(round(mean(dA),2)),'\pm',...
    num2str(round(conf(dA,0.95),2))],...
    ['drink perm: ',num2str(round(mean(dAR),2)),'\pm',...
    num2str(round(conf(dAR,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('behavior vs. not+other: subset features')
% Across models
% ROCs
figure
hold on
plot(mfdX,mfdY)
plot(mdfX,mdfY)
legend({['feed->drink: ',num2str(round(mean(fdA),2)),'\pm',...
    num2str(round(conf(fdA,0.95),2))],...
    ['drink->feed: ',num2str(round(mean(dfA),2)),'\pm',...
    num2str(round(conf(dfA,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('behavior vs. not+other; subset features')
%% Load single dataset models (food and drink) - data from runBingeOther 
% and runAlcoholOther - full models (all features); N.B. different features
% Preallocate
[fX,fY,fXR,fYR] = deal(zeros(100,142));
[dX,dY,dXR,dYR] = deal(zeros(100,94));
[fA,fAR,dA,dAR] = deal(zeros(1,100));
% Feed
for ii = 1:100
    load(['D:\paper3\analyzed\final\feedOtherAll\feedOtherAll',...
        num2str(ii),'.mat'],'hist','histR','acc','accR')
    [fX(ii,:),fY(ii,:),~,fA(ii)] = perfcurve(hist.cfg.naive.testY,...
        acc{1}.pred,1);
    [thisX,thisY,~,fAR(ii)] = perfcurve(histR.cfg.naive.testY,...
        accR{1}.pred,1);
    fXR(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,142));
    fYR(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,142));
end
% Drink
for ii = 1:100
    load(['D:\paper3\analyzed\final\drinkOtherAll\drinkOtherAll',...
        num2str(ii),'.mat'],'hist','histR','acc','accR')
    [thisX,thisY,~,dA(ii)] = perfcurve(hist.cfg.naive.testY,acc{1}.pred,...
        1);
    dX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,94));
    dY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,94));
    [thisX,thisY,~,dAR(ii)] = perfcurve(histR.cfg.naive.testY,...
        accR{1}.pred,1);
    dXR(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,94));
    dYR(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,94));
end
% save('F:\paper3\analyzed\final\feedAlcoholOtherFullModels.mat','fX',...
%     'fY','fA','fXR','fYR','fAR','dX','dY','dA','dXR','dYR','dAR')
% Plot ROCs
figure
hold on
plot(mean(fX,1),mean(fY,1))
plot(mean(fXR,1),mean(fXR,1))
plot(mean(dX,1),mean(dY,1))
plot(mean(dXR,1),mean(dYR,1))
legend({['feed: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],...
    ['feed perm: ',num2str(round(mean(fAR),2)),'\pm',...
    num2str(round(conf(fAR,0.95),2))],...
    ['drink: ',num2str(round(mean(dA),2)),'\pm',...
    num2str(round(conf(dA,0.95),2))],...
    ['drink perm: ',num2str(round(mean(dAR),2)),'\pm',...
    num2str(round(conf(dAR,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('behavior vs. not+other: all features')
% Plot AUC hists
figure
hold on
h(1) = histogram(fA,'normalization','probability','FaceAlpha',1);
h(2) = histogram(fAR,'normalization','probability','FaceAlpha',1);
h(3) = histogram(dA,'normalization','probability','FaceAlpha',1);
xax = get(gca,'xlim');
len = xax(2)-xax(1);
bin = len/40;
set(h(1),'BinWidth',bin);
set(h(2),'BinWidth',bin);
set(h(3),'BinWidth',bin);
legend({['feed: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],...
    ['feed perm: ',num2str(round(mean(fAR),2)),'\pm',...
    num2str(round(conf(fAR,0.95),2))],...
    ['drink: ',num2str(round(mean(dA),2)),'\pm',...
    num2str(round(conf(dA,0.95),2))]},'location','nw')
xlabel('AUC'); ylabel('probability')
title('behavior vs. not+other: all features')
%% Load single dataset models (food and drink) - data from runBingeOther 
% and runAlcoholOther - subset models (shared features)
% N.B.: run with LOO or without
% Preallocate
[fX,fY,fXR,fYR,dfX,dfY] = deal(zeros(100,142));
[dX,dY,dXR,dYR,fdX,fdY] = deal(zeros(100,94));
[fA,fAR,fdA,dA,dAR,dfA] = deal(zeros(1,100));
% Feed
for ii = 1:100
    load(['D:\paper3\analyzed\final\feedOther\feedOther',num2str(ii),...
        '.mat'])
%     load(['D:\paper3\analyzed\final\feedOther\feedOtherLOO',...
%         num2str(ii),'.mat'])
    [fX(ii,:),fY(ii,:),~,fA(ii)] = perfcurve(hist.cfg.naive.testY,...
        acc{1}.pred,1);
    [thisX,thisY,~,fAR(ii)] = perfcurve(histR.cfg.naive.testY,...
        accR{1}.pred,1);
    fXR(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,142));
    fYR(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,142));
    % Apply to drink data
    load(['D:\paper3\analyzed\final\drinkOther\drinkOther',num2str(ii),...
        '.mat'],'hist')
%     load(['D:\paper3\analyzed\final\drinkOtherlLOO\drinkOtherLOO',...
%             num2str(ii),'.mat'],'hist')
    pred = cvglmnetPredict(acc{1}.mdl{1},hist.cfg.naive.testX,...
        'lambda_1se','response');
    [fdX(ii,:),fdY(ii,:),~,fdA(ii)] = perfcurve(hist.cfg.naive.testY,...
        pred,1);
end
% Drink
for ii = 1:100
    load(['D:\paper3\analyzed\final\drinkOther\drinkOther',num2str(ii),...
        '.mat'])
%     load(['D:\paper3\analyzed\final\drinkOtherLOO\drinkOtherLOO',...
%         num2str(ii),'.mat'])
    [thisX,thisY,~,dA(ii)] = perfcurve(hist.cfg.naive.testY,acc{1}.pred,...
        1);
    dX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,94));
    dY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,94));
    [thisX,thisY,~,dAR(ii)] = perfcurve(histR.cfg.naive.testY,...
        accR{1}.pred,1);
    if isequal(thisX,[0;1])
        dXR(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
            linspace(0,1,94));
        dYR(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
            linspace(0,1,94));
    else
        dXR(ii,:) = thisX; dYR(ii,:) = thisY;
    end
    % Apply to feed data
    load(['D:\paper3\analyzed\final\feedOther\feedOther',num2str(ii),...
        '.mat'],'hist')
%     load(['D:\paper3\analyzed\final\feedOtherLOO\feedOtherLOO',...
%         num2str(ii),'.mat'],'hist')
    pred = cvglmnetPredict(acc{1}.mdl{1},hist.cfg.naive.testX,...
        'lambda_1se','response');
    [thisX,thisY,~,dfA(ii)] = perfcurve(hist.cfg.naive.testY,pred,1);
    dfX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,...
        142));
    dfY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,...
        142));
end
% save('F:\paper3\analyzed\final\feedAlcoholOtherSubsetModels.mat','fX',...
%     'fY','fA','fXR','fYR','fAR','fdX','fdY','fdA','dX','dY','dA','dXR',...
%     'dYR','dAR','dfX','dfY','dfA')
% save('F:\paper3\analyzed\final\feedAlcoholOtherSubsetModelsLOO.mat',...
%     'fX','fY','fA','fXR','fYR','fAR','fdX','fdY','fdA','dX','dY','dA',...
%     'dXR','dYR','dAR','dfX','dfY','dfA')
%% Plot alcoholOther and feedOther; subset models (shared features)
% Load either 80/20 or LOO
load('D:\paper3\analyzed\final\feedAlcoholOtherSubsetModels.mat')
% load('D:\paper3\analyzed\final\feedAlcoholOtherSubsetModelsLOO.mat')
% Within models
% ROCs
figure
hold on
plot(mean(fX,1),mean(fY,1))
plot(mean(fXR,1),mean(fXR,1))
plot(mean(dX,1),mean(dY,1))
plot(mean(dXR,1),mean(dYR,1))
legend({['feed: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],...
    ['feed perm: ',num2str(round(mean(fAR),2)),'\pm',...
    num2str(round(conf(fAR,0.95),2))],...
    ['drink: ',num2str(round(mean(dA),2)),'\pm',...
    num2str(round(conf(dA,0.95),2))],...
    ['drink perm: ',num2str(round(mean(dAR),2)),'\pm',...
    num2str(round(conf(dAR,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('behavior vs. not+other: subset features')
% AUC dists
figure
hold on
h(1) = histogram(fA,'normalization','probability','FaceAlpha',1);
h(2) = histogram(fAR,'normalization','probability','FaceAlpha',1);
h(3) = histogram(dA,'normalization','probability','FaceAlpha',1);
xax = get(gca,'xlim');
len = xax(2)-xax(1);
bin = len/40;
set(h(1),'BinWidth',bin);
set(h(2),'BinWidth',bin);
set(h(3),'BinWidth',bin);
legend({['feed: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],...
    ['feed perm: ',num2str(round(mean(fAR),2)),'\pm',...
    num2str(round(conf(fAR,0.95),2))],...
    ['drink: ',num2str(round(mean(dA),2)),'\pm',...
    num2str(round(conf(dA,0.95),2))]},'location','ne')
xlabel('AUC'); ylabel('probability')
title('behavior vs. not+other: subset features')
% Across models
% ROCs
figure
hold on
plot(mean(fdX),mean(fdY))
plot(mean(dfX),mean(dfY))
legend({['feed->drink: ',num2str(round(mean(fdA),2)),'\pm',...
    num2str(round(conf(fdA,0.95),2))],...
    ['drink->feed: ',num2str(round(mean(dfA),2)),'\pm',...
    num2str(round(conf(dfA,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('behavior vs. not+other; subset features')
% AUC dists
figure
hold on
h(1) = histogram(fdA,'normalization','probability','FaceAlpha',1);
h(2) = histogram(dfA,'normalization','probability','FaceAlpha',1);
xax = get(gca,'xlim');
len = xax(2)-xax(1);
bin = len/40;
set(h(1),'BinWidth',bin);
set(h(2),'BinWidth',bin);
legend({['feed->drink: ',num2str(round(mean(fdA),2)),'\pm',...
    num2str(round(conf(fdA,0.95),2))],...
    ['drink->feed: ',num2str(round(mean(dfA),2)),'\pm',...
    num2str(round(conf(dfA,0.95),2))]},'location','nw')
xlabel('AUC'); ylabel('probability')
title('behavior vs. not+other; subset features')
%% Single feature comparison - LOO
% Get single features of 1feedWeight (equal weight models)
% Preallocate
[s,a] = deal(zeros(81,18));
for ii = 1:81
    load(['D:/paper3/analyzed/final/genModel/1feedWeightLOO/'...
        'feedDrinkOtherLOO',num2str(ii),'.mat'],'sA','sSign')
    a(ii,:) = sA;
    s(ii,:) = sSign;
end
% [~,p] = ttest(a-.5);
% pAdj = p*18;
% sigInds = pAdj<0.05;
sM = mean(s);
aM = mean(a);
[aSort,sortInd] = sort(aM,'descend');
fFeat = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
feedShellInds = [1:12,25:30];
fFeat = fFeat(feedShellInds);
fFeat = fFeat(sortInd); %#ok
% Get single features of feedOther and drinkOther
[sF,fA,sD,dA] = deal([]);
for ii = 1:100
    load(['D:/paper3/analyzed/final/feedOtherLOO/feedOtherLOO',...
        num2str(ii),'.mat'],'s','sFA')
    % Average over the nine different models from leaving each animal out
    sF(ii,:) = mean(s,1);
    fA(ii,:) = mean(sFA,1);
    load(['D:/paper3/analyzed/final/drinkOtherLOO/drinkOtherLOO',...
        num2str(ii),'.mat'],'s','sDA')
    sD(ii,:) = mean(s,1);
    dA(ii,:) = mean(sDA,1);
end
[sortFA,sortFAind] = sort(mean(fA,1),'descend');
faFeat = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
faFeat = faFeat(feedShellInds);
faFeat = faFeat(sortFAind)';
[sortDA,sortDAind] = sort(mean(dA,1),'descend');
daFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
drinkShellInds = [13:24,55:60];
daFeat = daFeat(drinkShellInds);
daFeat = daFeat(sortDAind)';
sFM = sign(mean(sF,1));
sDM = sign(mean(sD,1));
% Highlight top two features from generalized model
sigInds = sortInd(1:2);
figure
hold on
h(1) = plot(mean(fA(:,1:12)).*sFM(1:12),mean(dA(:,1:12)).*sDM(1:12),'ok');
h(2) = plot(mean(fA(:,13:end)).*sFM(13:end),...
    mean(dA(:,13:end)).*sDM(13:end),'sk');
h(3) = scatter(mean(fA(:,sigInds)).*sFM(sigInds),...
    mean(dA(:,sigInds)).*sDM(sigInds),'o','filled');
plot([-1 1],[0.5 0.5],'k')
plot([-1 1],[-0.5 -0.5],'k')
plot([0.5 0.5],[-1 1],'k')
plot([-0.5 -0.5],[-1 1],'k')
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
xlabel('feed auc')
ylabel('drink auc')
title('single feature comparison: LOO')
legend(h,{'power','coherence','top gen features'},'location','north')
%% Single feature comparison - all (80/20)
% Get single features of 1feedWeight (equal weight models)
% Preallocate
[s,a] = deal(zeros(81,18));
for ii = 1:81
    load(['D:/paper3/analyzed/final/genModel/1feedWeight/'...
        'feedDrinkOther',num2str(ii),'.mat'],'sA','sSign')
    a(ii,:) = sA;
    s(ii,:) = sSign;
end
% [~,p] = ttest(a-.5);
% pAdj = p*18;
% sigInds = pAdj<0.05;
sM = mean(s);
aM = mean(a);
[aSort,sortInd] = sort(aM,'descend');
feat = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
feedShellInds = [1:12,25:30];
feat = feat(feedShellInds);
feat = feat(sortInd); %#ok
% Get single features of feedOther and drinkOther
[sF,fA,sD,dA] = deal([]);
for ii = 1:100
    load(['D:/paper3/analyzed/final/feedOtherAll/feedOtherAll',...
        num2str(ii),'.mat'],'s','sFA')
    % Average over the nine different models from leaving each animal out
    sF(ii,:) = mean(s,1);
    fA(ii,:) = mean(sFA,1);
    load(['D:/paper3/analyzed/final/drinkOtherAll/drinkOtherAll',...
        num2str(ii),'.mat'],'s','sDA')
    sD(ii,:) = mean(s,1);
    dA(ii,:) = mean(sDA,1);
end
fA = fA(:,feedShellInds);
[sortFA,sortFAind] = sort(mean(fA,1),'descend');
faFeat = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
faFeat = faFeat(feedShellInds);
faFeat = faFeat(sortFAind)';

drinkShellInds = [13:24,55:60];
dA = dA(:,drinkShellInds);
[sortDA,sortDAind] = sort(mean(dA,1),'descend');
daFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
daFeat = daFeat(drinkShellInds);
daFeat = daFeat(sortDAind)';

sFM = sign(mean(sF,1));
sFM = sFM(feedShellInds);
sDM = sign(mean(sD,1));
sDM = sDM(drinkShellInds);

% Highlight top two features from generalized model
sigInds = sortInd(1:2);
figure
hold on
h(1) = plot(mean(fA(:,1:12)).*sFM(1:12),mean(dA(:,1:12)).*sDM(1:12),'ok');
h(2) = plot(mean(fA(:,13:end)).*sFM(13:end),...
    mean(dA(:,13:end)).*sDM(13:end),'sk');
h(3) = scatter(mean(fA(:,sigInds)).*sFM(sigInds),...
    mean(dA(:,sigInds)).*sDM(sigInds),'o','filled');
plot([-1 1],[0.5 0.5],'k')
plot([-1 1],[-0.5 -0.5],'k')
plot([0.5 0.5],[-1 1],'k')
plot([-0.5 -0.5],[-1 1],'k')
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
xlabel('feed auc')
ylabel('drink auc')
title('single feature comparison')
legend(h,{'power','coherence','top gen features'},'location','north')
%% single features with permuted and groups
for ii = 1:100
    load(['H:\Shared drives\dwielDoucetteLab\data\paper3\analyzed\final\feedOther\','feedOther',num2str(ii),'.mat'])
    fA(ii) = acc{1}.acc;
    fAP(ii) = accR{1}.acc;
    fAs(ii,:) = sFA;
    fAsP(ii,:) = sFAR;
    gFA(ii,:) = groupFA;
    load(['H:\Shared drives\dwielDoucetteLab\data\paper3\analyzed\final\drinkOther\','drinkOther',num2str(ii),'.mat'])
    dA(ii) = acc{1}.acc;
    dAP(ii) = accR{1}.acc;
    dAs(ii,:) = sDA;
    dAsP(ii,:) = sDAR;
    gDA(ii,:) = groupDA;
end
%%
for ii = 1:18
    fP(ii) = (sum(fAsP(:,ii)>mean(fAs(:,ii)))+1)/(101);
    dP(ii) = (sum(dAsP(:,ii)>mean(dAs(:,ii)))+1)/(101);
end
%% feature groups
load('H:\Shared drives\dwielDoucetteLab\data\paper3\analyzed\final\24sweetChowNotDataNew.mat', 'R')
feedShellInds = [1:12,25:30];
Rfeed = R(feedShellInds,feedShellInds);
faFeat = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
faFeat = faFeat(feedShellInds);
for ii = 1:18
    groupFeed{ii} = faFeat(abs(Rfeed(:,ii))>0.5);
end
mgfa = mean(gFA)';
load('H:\Shared drives\dwielDoucetteLab\data\paper3\analyzed\final\waterAlcoholNotData.mat', 'R')
drinkShellInds = [13:24,55:60];
Rdrink = R(drinkShellInds,drinkShellInds);
daFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
daFeat = daFeat(drinkShellInds);
for ii = 1:18
    groupDrink{ii} = daFeat(abs(Rdrink(:,ii))>0.5);
end
mgda = mean(gDA)';
%% correlations within data sets - feeding vs. not; drinking vs. not
load('H:\Shared drives\dwielDoucetteLab\data\paper3\analyzed\final\24sweetChowNotDataNew.mat')
allRF = corrcoef(cat(1,allFeedData{:,1},allFeedData{:,2},...
    allFeedData{:,3}));
sfR = corrcoef(cat(1,allFeedData{:,1}));
nfR = corrcoef(cat(1,allFeedData{:,3}));
for k = 1:200
    for ii = 1:9
        thisN = min(cellfun(@(x) size(x,1),allFeedData(ii,:)));
        thisSub{ii,1} = allFeedData{ii,1}(randperm(size(...
            allFeedData{ii,1},1),thisN),:);
        thisSub{ii,2} = allFeedData{ii,2}(randperm(size(...
            allFeedData{ii,2},1),thisN),:);
        thisSub{ii,3} = allFeedData{ii,3}(randperm(size(...
            allFeedData{ii,3},1),thisN),:);
    end
    subAllRF(:,:,k) = corrcoef(cat(1,thisSub{:,1},thisSub{:,2},...
        thisSub{:,3}));
    subSFR(:,:,k) = corrcoef(cat(1,thisSub{:,1}));
    subNFR(:,:,k) = corrcoef(cat(1,thisSub{:,3}));
end
%%
feedShellInds = [1:12,25:30];
figure
subplot(1,3,1)
pcolor(padarray(allRF(feedShellInds,feedShellInds),[1 1],'post'))
title('all')
subplot(1,3,2)
pcolor(padarray(sfR(feedShellInds,feedShellInds),[1 1],'post'))
title('sf')
subplot(1,3,3)
pcolor(padarray(nfR(feedShellInds,feedShellInds),[1 1],'post'))
title('not')
colormap('viridis')
figure
subplot(1,3,1)
pcolor(padarray(subAllRF(feedShellInds,feedShellInds),[1 1],'post'))
title('sub all')
subplot(1,3,2)
pcolor(padarray(subSFR(feedShellInds,feedShellInds),[1 1],'post'))
title('sub sf')
subplot(1,3,3)
pcolor(padarray(subNFR(feedShellInds,feedShellInds),[1 1],'post'))
title('sub not')
colormap('viridis')
%%
load('H:\Shared drives\dwielDoucetteLab\data\paper3\analyzed\final\waterAlcoholNotData.mat')
%%
drinkShellInds = [13:24,55:60];
this = cat(1,allData{:,1},allData{:,2},allData{:,3});
allRD = corrcoef(table2array(this(:,drinkShellInds)));
this = cat(1,allData{:,1});
edR = corrcoef(table2array(this(:,drinkShellInds)));
this = cat(1,allData{:,2});
wdR = corrcoef(table2array(this(:,drinkShellInds)));
this = cat(1,allData{:,3});
ndR = corrcoef(table2array(this(:,drinkShellInds)));
c = 1;
this = cell(1,3);
for k = 1:200
    for ii = 1:9
        thisN = min(cellfun(@(x) size(x,1),allData(ii,:)));
        thisSub{ii,1} = table2array(allData{ii,1}(randperm(size(...
            allData{ii,1},1),thisN),drinkShellInds));
        thisSub{ii,2} = table2array(allData{ii,2}(randperm(size(...
            allData{ii,2},1),thisN),drinkShellInds));
        thisSub{ii,3} = table2array(allData{ii,3}(randperm(size(...
            allData{ii,3},1),thisN),drinkShellInds));
    end
    subAllRD(:,:,k) = corrcoef(cat(1,thisSub{:,1},thisSub{:,2},...
        thisSub{:,3}));
    subEDR(:,:,k) = corrcoef(cat(1,thisSub{:,1}));
    subWDR(:,:,k) = corrcoef(cat(1,thisSub{:,2}));
    subNDR(:,:,k) = corrcoef(cat(1,thisSub{:,3}));
    for r = 1:25
        for ii = 1:9
            for jj = 1:18
                for m = 1:3
                    this{ii,m}(:,jj) = thisSub{ii,m}(randperm(size(...
                        thisSub{ii,m},1)),jj);
                end
            end
        end
        subAllRDPerm(:,:,c) = corrcoef(cat(1,this{:,1},this{:,2},this{:,3}));
        subEDRPerm(:,:,c) = corrcoef(cat(1,this{:,1}));
        subWDRPerm(:,:,c) = corrcoef(cat(1,this{:,2}));
        subNDRPerm(:,:,c) = corrcoef(cat(1,this{:,3}));
        c=c+1;
    end
%     subEvN(:,:,k) = corr(subNDR(:,:,k),subEDR(:,:,k),'type','spearman');
%     subEvW(:,:,k) = corr(subEDR(:,:,k),subWDR(:,:,k),'type','spearman');
%     for jj = 1:25
%         thisR = randperm(18);
%         subEvNperm(:,:,c) = corr(subNDR(thisR,thisR,k),...
%             subEDR(:,:,k),'type','spearman');
%         subEvWPerm(:,:,c) = corr(subEDR(thisR,thisR,k),...
%             subWDR(:,:,k),'type','spearman');
%         c = c+1;
%     end
end
%% remove 1's and only look at upper half
subAllRDPerm(subAllRDPerm == 1) = NaN;
subEDRPerm(subEDRPerm == 1) = NaN;
subWDRPerm(subWDRPerm == 1) = NaN;
subNDRPerm(subNDRPerm == 1) = NaN;

subAllRDPermUpper = [];
subEDRPermUpper = [];
subWDRPermUpper = [];
subNDRPermUpper = [];
for ii = 1:5000
    subAllRDPermUpper(:,:,ii) = triu(subAllRDPerm(:,:,ii),1);
    subEDRPermUpper(:,:,ii) = triu(subEDRPerm(:,:,ii),1);
    subWDRPermUpper(:,:,ii) = triu(subAllRDPerm(:,:,ii),1);
    subNDRPermUpper(:,:,ii) = triu(subNDRPerm(:,:,ii),1);
end
subAllRDPermUpper(subAllRDPermUpper==0) = NaN;
subEDRPermUpper(subEDRPermUpper==0) = NaN;
subWDRPermUpper(subWDRPermUpper==0) = NaN;
subNDRPermUpper(subNDRPermUpper==0) = NaN;
%%
figure
subplot(1,4,1)
histogram(reshape(subAllRDPermUpper,1,numel(subAllRDPermUpper)),'NumBins',30);
xlim([-0.55 0.55])
subplot(1,4,2)
histogram(reshape(subEDRPermUpper,1,numel(subEDRPermUpper)),'NumBins',30);
xlim([-0.55 0.55])
subplot(1,4,3)
histogram(reshape(subWDRPermUpper,1,numel(subWDRPermUpper)),'NumBins',30);
xlim([-0.55 0.55])
subplot(1,4,4)
histogram(reshape(subNDRPermUpper,1,numel(subNDRPermUpper)),'NumBins',30);
xlim([-0.55 0.55])
%%
figure
subplot(1,4,1)
this = allRD;
this(this==1) = NaN;
this(abs(this)<0.5) = NaN;
pcolor(padarray(this,[1 1],'post'))
title('all')
subplot(1,4,2)
this = edR;
this(this==1) = NaN;
this(abs(this)<0.5) = NaN;
pcolor(padarray(this,[1 1],'post'))
title('etoh')
subplot(1,4,3)
this = wdR;
this(this==1) = NaN;
this(abs(this)<0.5) = NaN;
pcolor(padarray(this,[1 1],'post'))
title('water')
subplot(1,4,4)
this = ndR;
this(this==1) = NaN;
this(abs(this)<0.5) = NaN;
pcolor(padarray(this,[1 1],'post'))
title('not')
colormap('viridis')
%%
figure
subplot(1,4,1)
this = mean(subAllRD,3);
this(abs(this)<0.5) = NaN;
this(this == 1) = NaN;
pcolor(padarray(this,[1 1],'post'))
title('sub all')
subplot(1,4,2)
this = mean(subEDR,3);
this(abs(this)<0.5) = NaN;
this(this == 1) = NaN;
pcolor(padarray(this,[1 1],'post'))
title('sub etoh')
subplot(1,4,3)
this = mean(subWDR,3);
this(abs(this)<0.5) = NaN;
this(this == 1) = NaN;
pcolor(padarray(this,[1 1],'post'))
title('sub water')
subplot(1,4,4)
this = mean(subNDR,3);
this(abs(this)<0.5) = NaN;
this(this == 1) = NaN;
pcolor(padarray(this,[1 1],'post'))
title('sub not')
colormap('viridis')
%% differences
% ethanol from water
for ii = 1:5000
    evwPerm(:,:,ii) = subWDRPerm(:,:,ii)-subEDRPerm(:,:,ii);
    evwPermUpper(:,:,ii) = triu(evwPerm(:,:,ii),1);
end
evwPermUpper(evwPermUpper==0) = NaN;
figure
histogram(reshape(evwPermUpper,1,numel(evwPermUpper)),'NumBins',30)
x = 0.1:0.1:1;
for ii = 1:numel(x)
    ps(ii) = (sum(abs(evwPermUpper)>x(ii),[1,2,3])+1)/5001;
end
evw = mean(subWDR,3)-mean(subEDR,3);
figure
this = evw;
this(abs(this)<=0.7) = NaN;
pcolor(padarray(this,[1,1],'post'))
colormap viridis
%%
drinkShellInds = [13:24,55:60];
[sortDA,sortDAind] = sort(mean(dAs,1),'descend');
daFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
daFeat = daFeat(drinkShellInds);
daFeat = daFeat(sortDAind)';

feedShellInds = [1:12,25:30];
[sortFA,sortFAind] = sort(mean(fAs,1),'descend');
faFeat = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
faFeat = faFeat(feedShellInds);
faFeat = faFeat(sortFAind)';
%% alcohol permuted features
a2 = [];
for ii = 1:153
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm2/drinkOtherPerm2_',num2str(ii),'.mat'],...
        'acc')
    a2(ii) = acc{1}.acc;
end
a3 = [];
for ii = 1:816
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm3/drinkOtherPerm3_',num2str(ii),'.mat'],...
        'acc')
    a3(ii) = acc{1}.acc;
end
a4 = [];
for ii = 1:3060
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm4/drinkOtherPerm4_',num2str(ii),'.mat'],...
        'acc')
    a4(ii) = acc{1}.acc;
end
a5 = [];
for ii = 1:8568
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm5/drinkOtherPerm5_',num2str(ii),'.mat'],...
        'acc')
    a5(ii) = acc{1}.acc;
end
a6 = [];
c = 1;
for ii = 1:10
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm6/drinkOtherPerm6_',num2str(ii),'.mat'],...
        'acc')
    for jj = 1:numel(acc)
        a6(c) = acc{jj}{1}.acc;
        c = c+1;
    end
end
a7 = [];
c = 1;
for ii = 1:10
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm7/drinkOtherPerm7_',num2str(ii),'.mat'],...
        'acc')
    for jj = 1:numel(acc)
        a7(c) = acc{jj}{1}.acc;
        c = c+1;
    end
end
a8 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm8/drinkOtherPerm8_',num2str(ii),'.mat'],...
        'acc')
    a8(ii) = acc{1}.acc;
end
a9 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm9/drinkOtherPerm9_',num2str(ii),'.mat'],...
        'acc')
    a9(ii) = acc{1}.acc;
end
a10 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm10/drinkOtherPerm10_',num2str(ii),'.mat'],...
        'acc')
    a10(ii) = acc{1}.acc;
end
a11 = [];
for ii = 1:1000
    disp(ii)
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm11/drinkOtherPerm11_',num2str(ii),'.mat'],...
        'acc')
    a11(ii) = acc{1}.acc;
end
%%
a12 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm12/drinkOtherPerm12_',num2str(ii),'.mat'],...
        'acc')
    a12(ii) = acc{1}.acc;
end
%%
a13 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm13/drinkOtherPerm13_',num2str(ii),'.mat'],...
        'acc')
    a13(ii) = acc{1}.acc;
end
%%
a14 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/drinkOtherPerm14/drinkOtherPerm14_',num2str(ii),'.mat'],...
        'acc')
    a14(ii) = acc{1}.acc;
end
%%
figure
hold on
plot(2:14,[mean(a2) mean(a3) mean(a4) mean(a5) mean(a6) ...
    mean(a7) mean(a8) mean(a9) mean(a10) mean(a11) mean(a12) ...
    mean(a13) mean(a14)],'ok')
plot([2 2],[mean(a2)-2*std(a2) mean(a2)+2*std(a2)],'-k')
plot([3 3],[mean(a3)-2*std(a3) mean(a3)+2*std(a3)],'-k')
plot([4 4],[mean(a4)-2*std(a4) mean(a4)+2*std(a4)],'-k')
plot([5 5],[mean(a5)-2*std(a5) mean(a5)+2*std(a5)],'-k')
plot([6 6],[mean(a6)-2*std(a6) mean(a6)+2*std(a6)],'-k')
plot([7 7],[mean(a7)-2*std(a7) mean(a7)+2*std(a7)],'-k')
plot([8 8],[mean(a8)-2*std(a8) mean(a8)+2*std(a8)],'-k')
plot([9 9],[mean(a9)-2*std(a9) mean(a9)+2*std(a9)],'-k')
plot([10 10],[mean(a10)-2*std(a10) mean(a10)+2*std(a10)],'-k')
plot([11 11],[mean(a11)-2*std(a11) mean(a11)+2*std(a11)],'-k')
plot([12 12],[mean(a12)-2*std(a12) mean(a12)+2*std(a12)],'-k')
plot([13 13],[mean(a13)-2*std(a13) mean(a13)+2*std(a13)],'-k')
plot([14 14],[mean(a14)-2*std(a14) mean(a14)+2*std(a14)],'-k')
%% feeding permuted features
f2 = [];
for ii = 1:153
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm2/feedOtherPerm2_',num2str(ii),'.mat'],...
        'acc')
    f2(ii) = acc{1}.acc;
end
f3 = [];
for ii = 1:816
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm3/feedOtherPerm3_',num2str(ii),'.mat'],...
        'acc')
    f3(ii) = acc{1}.acc;
end
f4 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm4/feedOtherPerm4_',num2str(ii),'.mat'],...
        'acc')
    f4(ii) = acc{1}.acc;
end
f5 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm5/feedOtherPerm5_',num2str(ii),'.mat'],...
        'acc')
    f5(ii) = acc{1}.acc;
end
f6 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm6/feedOtherPerm6_',num2str(ii),'.mat'],...
        'acc')
    f6(ii) = acc{1}.acc;
end
f7 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm7/feedOtherPerm7_',num2str(ii),'.mat'],...
        'acc')
    f7(ii) = acc{1}.acc;
end
f8 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm8/feedOtherPerm8_',num2str(ii),'.mat'],...
        'acc')
    f8(ii) = acc{1}.acc;
end
f9 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm9/feedOtherPerm9_',num2str(ii),'.mat'],...
        'acc')
    f9(ii) = acc{1}.acc;
end
f10 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm10/feedOtherPerm10_',num2str(ii),'.mat'],...
        'acc')
    f10(ii) = acc{1}.acc;
end
f11 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm11/feedOtherPerm11_',num2str(ii),'.mat'],...
        'acc')
    f11(ii) = acc{1}.acc;
end
f12 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm12/feedOtherPerm12_',num2str(ii),'.mat'],...
        'acc')
    f12(ii) = acc{1}.acc;
end
f13 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm13/feedOtherPerm13_',num2str(ii),'.mat'],...
        'acc')
    f13(ii) = acc{1}.acc;
end
f14 = [];
for ii = 1:1000
    load(['H:/Shared drives/dwielDoucetteLab/data/paper3/analyzed/',...
        'final/feedOtherPerm14/feedOtherPerm14_',num2str(ii),'.mat'],...
        'acc')
    f14(ii) = acc{1}.acc;
end
%%
figure
hold on
plot(2:14,[mean(f2) mean(f3) mean(f4) mean(f5) mean(f6) ...
    mean(f7) mean(f8) mean(f9) mean(f10) mean(f11) mean(f12) ...
    mean(f13) mean(f14)],'ok')
plot([2 2],[mean(f2)-2*std(f2) mean(f2)+2*std(f2)],'-k')
plot([3 3],[mean(f3)-2*std(f3) mean(f3)+2*std(f3)],'-k')
plot([4 4],[mean(f4)-2*std(f4) mean(f4)+2*std(f4)],'-k')
plot([5 5],[mean(f5)-2*std(f5) mean(f5)+2*std(f5)],'-k')
plot([6 6],[mean(f6)-2*std(f6) mean(f6)+2*std(f6)],'-k')
plot([7 7],[mean(f7)-2*std(f7) mean(f7)+2*std(f7)],'-k')
plot([8 8],[mean(f8)-2*std(f8) mean(f8)+2*std(f8)],'-k')
plot([9 9],[mean(f9)-2*std(f9) mean(f9)+2*std(f9)],'-k')
plot([10 10],[mean(f10)-2*std(f10) mean(f10)+2*std(f10)],'-k')
plot([11 11],[mean(f11)-2*std(f11) mean(f11)+2*std(f11)],'-k')
plot([12 12],[mean(f12)-2*std(f12) mean(f12)+2*std(f12)],'-k')
plot([13 13],[mean(f13)-2*std(f13) mean(f13)+2*std(f13)],'-k')
plot([14 14],[mean(f14)-2*std(f14) mean(f14)+2*std(f14)],'-k')
%% Dep 24 vs. chow
% load('D:\paper3\24sweetChowNotData.mat')
% [thisX,thisW,thisY] = deal(cell(9,2));
% % Grab 40 samples of sweet and chow; weight as needed
% for ii = 1:9
%     for jj = 1:2
%         thisN = size(allFeedData{ii,jj},1);
%         if thisN<40
%             thisX{ii,jj} = allFeedData{ii,jj};
%             thisW{ii,jj} = repmat(40/thisN,thisN,1);
%             if jj == 1
%                 thisY{ii,jj} = ones(thisN,1);
%             else
%                 thisY{ii,jj} = zeros(thisN,1);
%             end
%         else
%             thisX{ii,jj} = allFeedData{ii,jj}(randperm(thisN,40),:);
%             thisW{ii,jj} = ones(40,1);
%             % First column is sweet/fat food
%             if jj == 1
%                 thisY{ii,jj} = ones(40,1);
%             else
%                 thisY{ii,jj} = zeros(40,1);
%             end
%         end
%     end
% end
% for ii = 1:9
%     % LOO
%     trainInd = ~ismember(1:9,ii);
%     trainX = cat(1,thisX{trainInd,1},thisX{trainInd,2});
%     trainY = cat(1,thisY{trainInd,1},thisY{trainInd,2});
%     trainW = cat(1,thisW{trainInd,1},thisW{trainInd,2});
%     testX = cat(1,thisX{ii,1},thisX{ii,2});
%     testY = cat(1,thisY{ii,1},thisY{ii,2});
%     cfg = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se',trainW);
%     [~,lam,beta,fits,acc{ii},hist{ii}] = lassoNet(trainX,trainY,...
%         'binomial','class',1,10,1,cfg);
%     
%     cfg = lassoNetCfg({testX,testY},[],'y','y','n',100,'1se',trainW);
%     [~,lamP,betaP,fitsP,accP{ii},histP{ii}] = lassoNet(trainX,trainY,...
%         'binomial','class',1,10,1,cfg);
% end
% save(['~/data/depChow',num2str(n),'.mat'],'acc','hist','accP','histP')
%% Plot dep vs. chow
% for ii = 1:100
%     load(['D:\paper3\analyzed\final\depChow\depChow',num2str(ii),'.mat'])
%     for jj = 1:9
%         a(ii,jj) = acc{jj}{1}.acc;
%         x(ii,jj,:) = interp1(linspace(0,1,numel(acc{jj}{1}.x)),...
%             acc{jj}{1}.x,linspace(0,1,81));
%         y(ii,jj,:) = interp1(linspace(0,1,numel(acc{jj}{1}.y)),...
%             acc{jj}{1}.y,linspace(0,1,81));
%         aP(ii,jj) = accP{jj}{1}.acc;
%         xP(ii,jj,:) = interp1(linspace(0,1,numel(accP{jj}{1}.x)),...
%             accP{jj}{1}.x,linspace(0,1,81));
%         yP(ii,jj,:) = interp1(linspace(0,1,numel(accP{jj}{1}.y)),...
%             accP{jj}{1}.y,linspace(0,1,81));
%     end
% end
% %
% figure
% hold on
% plot(squeeze(mean(mean(x,1),2)),squeeze(mean(mean(y,1),2)))
% plot(squeeze(mean(mean(xP,1),2)),squeeze(mean(mean(yP,1),2)))
% legend({['LOO: ',num2str(round(mean(mean(a)),2)),'\pm',...
%     num2str(round(conf(reshape(a,1,900),.95),2))],['permuted: ',...
%     num2str(round(mean(mean(aP)),2)),'\pm',...
%     num2str(round(conf(reshape(aP,1,900),.95),2))]})
% title('sweet vs. chow')
% set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
% box off
%% Water vs. alcohol
% load('waterAlcoholNotData2.mat','allData')
% [thisX,thisW,thisY] = deal(cell(9,2));
% % Grab 40 samples of sweet and chow; weight as needed
% for ii = 1:9
%     for jj = 1:2
%         thisN = size(allData{ii,jj},1);
%         if thisN<40
%             thisX{ii,jj} = allFeedData{ii,jj};
%             thisW{ii,jj} = repmat(40/thisN,thisN,1);
%             if jj == 1
%                 thisY{ii,jj} = ones(thisN,1);
%             else
%                 thisY{ii,jj} = zeros(thisN,1);
%             end
%         else
%             thisX{ii,jj} = allFeedData{ii,jj}(randperm(thisN,40),:);
%             thisW{ii,jj} = ones(40,1);
%             % First column is water
%             if jj == 1
%                 thisY{ii,jj} = zeros(40,1);
%             else
%                 thisY{ii,jj} = ones(40,1);
%             end
%         end
%     end
%     % LOO
%     trainInd = ~ismember(1:9,ii);
%     trainX = cat(1,thisX{trainInd,1},thisX{trainInd,2});
%     trainY = cat(1,thisY{trainInd,1},thisY{trainInd,2});
%     trainW = cat(1,thisW{trainInd,1},thisW{trainInd,2});
%     testX = cat(1,thisX{ii,1},thisX{ii,2});
%     testY = cat(1,thisY{ii,1},thisY{ii,2});
%     cfg = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se',trainW);
%     [~,lam,beta,fits,acc{ii},hist{ii}] = lassoNet(trainX,trainY,...
%         'binomial','class',1,10,1,cfg);
%     
%     cfg = lassoNetCfg({testX,testY},[],'y','y','n',100,'1se',trainW);
%     [~,lamP,betaP,fitsP,accP{ii},histP{ii}] = lassoNet(trainX,trainY,...
%         'binomial','class',1,10,1,cfg);
% end
%% Plot alc vs. water
% for ii = 1:100
%     load(['D:\paper3\analyzed\final\alcWater\waterAlcohol',num2str(ii),...
%         '.mat'])
%     for jj = 1:9
%         a(ii,jj) = acc{jj}{1}.acc;
%         x(ii,jj,:) = interp1(linspace(0,1,numel(acc{jj}{1}.x)),...
%             acc{jj}{1}.x,linspace(0,1,81));
%         y(ii,jj,:) = interp1(linspace(0,1,numel(acc{jj}{1}.y)),...
%             acc{jj}{1}.y,linspace(0,1,81));
%         aP(ii,jj) = accP{jj}{1}.acc;
%         xP(ii,jj,:) = interp1(linspace(0,1,numel(accP{jj}{1}.x)),...
%             accP{jj}{1}.x,linspace(0,1,81));
%         yP(ii,jj,:) = interp1(linspace(0,1,numel(accP{jj}{1}.y)),...
%             accP{jj}{1}.y,linspace(0,1,81));
%     end
% end
% %
% figure
% hold on
% plot(squeeze(mean(mean(x,1),2)),squeeze(mean(mean(y,1),2)))
% plot(squeeze(mean(mean(xP,1),2)),squeeze(mean(mean(yP,1),2)))
% legend({['LOO: ',num2str(round(mean(mean(a)),2)),'\pm',...
%     num2str(round(conf(reshape(a,1,900),.95),2))],['permuted: ',...
%     num2str(round(mean(mean(aP)),2)),'\pm',...
%     num2str(round(conf(reshape(aP,1,900),.95),2))]})
% title('alcohol vs. water')
% set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
% box off
%% Build drinkNot models split by sex and tested across sex with leave one 
% animal out testing
[data,~,~] = collateData('F:\paper3\drinkNot\',{'.mat'},...
    {'pow','coh'},'trl','rel');
% Split
male = data{1}([1:18,24:25],:);
female = data{1}([21:23,26:45],:);
% Set up animal number vector
animal{1} = [1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5];
animal{2} = [1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5];
% save('D:\paper3\analyzed\final\maleFemale.mat','male','female','animal')
%% Male female drink models
cd D:/paper3/analyzed/final/maleFemaleDrinkSubset/
drinkShellInds = [13:24,55:60];
for ii = 1:100
    load(['maleFemale',num2str(ii),'.mat'])
    % Fix typo
%     maleMod80.testX = maleMod80.textX;
    % Female
    female80AUC(ii) = femaleMod80.a;
    female80x(ii,:) = femaleMod80.x;
    female80y(ii,:) = femaleMod80.y;
    % Female -> male
    prob = cvglmnetPredict(femaleMod80.mdl,...
        zscore(table2array(maleMod80.testX(:,drinkShellInds))),...
        'lambda_1se','response');
    [x,y,~,fm80A(ii)] = perfcurve(maleMod80.testY,prob,1);
    fm80X(ii,:) = interp1(linspace(0,1,numel(x)),x,linspace(0,1,30));
    fm80Y(ii,:) = interp1(linspace(0,1,numel(y)),y,linspace(0,1,30));
    % LOO
    for jj = 1:10
        % Female
        femaleLOOAUC(jj,ii) = femaleMod{jj}.auc;
        femaleLOOX(jj,ii,:) = interp1(linspace(0,1,...
            numel(femaleMod{jj}.x)),femaleMod{jj}.x,linspace(0,1,30));
        femaleLOOY(jj,ii,:) = interp1(linspace(0,1,...
            numel(femaleMod{jj}.y)),femaleMod{jj}.y,linspace(0,1,30));
        % Female -> male
        fmA(ii,jj,:) = femaleMod{jj}.mAuc;
        for k = 1:10
            thisX(k,:) = interp1(linspace(0,1,...
                numel(femaleMod{jj}.mX{k})),...
                femaleMod{jj}.mX{k},linspace(0,1,30));
            thisY(k,:) = interp1(linspace(0,1,...
                numel(femaleMod{jj}.mY{k})),...
                femaleMod{jj}.mY{k},linspace(0,1,30));
        end
        fmX(ii,jj,:) = mean(thisX,1);
        fmY(ii,jj,:) = mean(thisY,1);
        % Male
        maleLOOAUC(jj,ii) = maleMod{jj}.auc;
        maleLOOX(jj,ii,:) = interp1(linspace(0,1,...
            numel(maleMod{jj}.x)),maleMod{jj}.x,linspace(0,1,30));
        maleLOOY(jj,ii,:) = interp1(linspace(0,1,...
            numel(maleMod{jj}.y)),maleMod{jj}.y,linspace(0,1,30));
        % Male -> female
        mfA(ii,jj,:) = maleMod{jj}.fAuc;
        for k = 1:10
            thisX(k,:) = interp1(linspace(0,1,...
                numel(maleMod{jj}.fX{k})),maleMod{jj}.fX{k},...
                linspace(0,1,30));
            thisY(k,:) = interp1(linspace(0,1,...
                numel(maleMod{jj}.fY{k})),maleMod{jj}.fY{k},...
                linspace(0,1,30));
        end
        mfX(ii,jj,:) = mean(thisX,1);
        mfY(ii,jj,:) = mean(thisY,1);
    end
    % Male
    male80AUC(ii) = maleMod80.a;
    male80x(ii,:) = maleMod80.x;
    male80y(ii,:) = maleMod80.y;   
    % Male -> female
    prob = cvglmnetPredict(maleMod80.mdl,zscore(table2array(...
        femaleMod80.testX(:,drinkShellInds))),'lambda_1se','response');
    [x,y,~,mf80A(ii)] = perfcurve(femaleMod80.testY,prob,1);
    mf80X(ii,:) = interp1(linspace(0,1,numel(x)),x,linspace(0,1,30));
    mf80Y(ii,:) = interp1(linspace(0,1,numel(y)),y,linspace(0,1,30));
end
%% Permuted
cd('D:/paper3/analyzed/final/maleFemaleDrinkPermutedSubset')
clear thisX thisY
for ii = 1:252
    load(['maleFemaleP',num2str(ii),'.mat'])
    fP(ii,:) = femaleP.a;
    for jj = 1:10
        thisX(jj,:) = interp1(linspace(0,1,numel(femaleP.x{jj})),...
            femaleP.x{jj},linspace(0,1,20));
        thisY(jj,:) = interp1(linspace(0,1,numel(femaleP.y{jj})),...
            femaleP.y{jj},linspace(0,1,20));
    end
    fPX(ii,:) = mean(thisX,1);
    fPY(ii,:) = mean(thisY,1);
    f80P(ii) = female80P.a;
    f80PX(ii,:) = interp1(linspace(0,1,numel(female80P.x)),female80P.x,linspace(0,1,30));
    f80PY(ii,:) = interp1(linspace(0,1,numel(female80P.y)),female80P.y,linspace(0,1,30));

    mP(ii,:) = maleP.a;
    for jj = 1:10
        thisX(jj,:) = interp1(linspace(0,1,numel(maleP.x{jj})),maleP.x{jj},linspace(0,1,20));
        thisY(jj,:) = interp1(linspace(0,1,numel(maleP.y{jj})),maleP.y{jj},linspace(0,1,20));
    end
    mPX(ii,:) = mean(thisX,1);
    mPY(ii,:) = mean(thisY,1);
    m80P(ii) = female80P.a;
    m80PX(ii,:) = interp1(linspace(0,1,numel(male80P.x)),male80P.x,linspace(0,1,30));
    m80PY(ii,:) = interp1(linspace(0,1,numel(male80P.y)),male80P.y,linspace(0,1,30));
end
%% apply model with both sexes to each sex
cd D:\paper3\analyzed\final\drinkOtherNew\
maleID = {'AH11','AH12','AH16','AH17','AH20','AH38','AH39','AH87',...
    'AH88','AH90'};
femaleID = {'AH1','AH3','AH4','AH5','AH8','AH26','AH29','AH30','AH80',...
    'AH81'};
for ii = 1:100
    % load general model
    load(['drinkOther',num2str(ii),'.mat'])
    fInds = [];
    for jj = 1:numel(femaleID)
        fInds(:,jj) = strcmp(animalID,femaleID{jj});
    end
    fInds = any(fInds,2);
    mInds = [];
    for jj = 1:numel(maleID)
        mInds(:,jj) = strcmp(animalID,maleID{jj});
    end
    mInds = any(mInds,2);
    % gen -> F
    prob = cvglmnetPredict(acc{1}.mdl{1},zscore(hist.cfg.naive.testX(fInds,:)),...
        'lambda_1se','response');
    [thisX,thisY,~,fA(ii)] = perfcurve(hist.cfg.naive.testY(fInds,:),...
        prob,1);
    fX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,50));
    fY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,50));
    % gen -> M
    prob = cvglmnetPredict(acc{1}.mdl{1},zscore(hist.cfg.naive.testX(mInds,:)),...
        'lambda_1se','response');
    [thisX,thisY,~,mA(ii)] = perfcurve(hist.cfg.naive.testY(mInds,:),...
        prob,1);
    mX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,50));
    mY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,50));
    a(ii) = acc{1}.acc;
    aR(ii) = accR{1}.acc;
end
%%
figure
hold on
plot(mean(fX),mean(fY))
plot(mean(mX),mean(mY))
legend({['female: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2)),' ',...
    num2str((sum(aR>mean(fA))+1)/(numel(aR)+1))],['male: ',...
    num2str(round(mean(mA),2)),'\pm',num2str(round(conf(mA,0.95),2)),...
    ' ',num2str((sum(aR>mean(mA))+1)/(numel(aR)+1))]})
%% male and female plots
% single plot - 80:20
figure
hold on
plot(mean(female80x),mean(female80y),'--k')
plot(mean(f80PX),mean(f80PY),'-k')
plot(mean(male80x),mean(male80y),'--','color',[0.7,0.7,0.7])
plot(mean(m80PX),mean(m80PY),'color',[0.7,0.7,0.7])
plot(mean(fX),mean(fY),'b')
plot(mean(mX),mean(mY),'r')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['female: ',num2str(round(mean(female80AUC),2)),'\pm',...
    num2str(round(conf(female80AUC,0.95),2)),' ',...
    num2str((sum(f80P>mean(female80AUC))+1)/(numel(f80P)+1))],...
    ['f perm: ',num2str(round(mean(f80P),2)),'\pm',...
    num2str(round(conf(f80P,0.95),2))],...
    ['male: ',num2str(round(mean(male80AUC),2)),'\pm',...
    num2str(round(conf(male80AUC,0.95),2)),' ',...
    num2str((sum(m80P>mean(male80AUC))+1)/(numel(m80P)+1))],...
    ['m perm: ',num2str(round(mean(m80P),2)),'\pm',...
    num2str(round(conf(m80P,0.95),2))],...
    ['drink->female: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2)),' ',...
    num2str((sum(aR>mean(fA))+1)/(numel(aR)+1))],['drink->male: ',...
    num2str(round(mean(mA),2)),'\pm',num2str(round(conf(mA,0.95),2)),...
    ' ',num2str((sum(aR>mean(mA))+1)/(numel(aR)+1))]})
box off
%%
figure
hold on
plot(mean(female80x),mean(female80y))
plot(mean(f80PX),mean(f80PY))
title(['female 80:20; p = ',...
    num2str((sum(f80P>mean(female80AUC))+1)/(numel(f80P)+1))])
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['female: ',num2str(round(mean(female80AUC),2)),'\pm',...
    num2str(round(conf(female80AUC,0.95),2))],['permuted: ',...
    num2str(round(mean(f80P),2)),'\pm',...
    num2str(round(conf(f80P,0.95),2))]})
box off

figure
hold on
plot(squeeze(mean(femaleLOOX,[1,2])),squeeze(mean(femaleLOOY,[1,2])))
plot(mean(fPX),mean(fPY))
title(['female LOO; p = ',...
    num2str((sum(reshape(fP,1,numel(fP))>mean(femaleLOOAUC,[1,2]))+1)/...
    (numel(fP)+1))])
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['female: ',num2str(round(mean(femaleLOOAUC,[1,2]),2)),'\pm',...
    num2str(round(conf(reshape(femaleLOOAUC,1,numel(femaleLOOAUC)),0.95)...
    ,2))],['permuted: ',num2str(round(mean(fP,[1,2]),2)),'\pm',...
    num2str(round(conf(reshape(fP,1,numel(fP)),0.95),2))]})
box off

figure
hold on
plot(mean(male80x),mean(male80y))
plot(mean(m80PX),mean(m80PY))
title(['male 80:20; p = ',...
    num2str((sum(m80P>mean(male80AUC))+1)/(numel(m80P)+1))])
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['male: ',num2str(round(mean(male80AUC),2)),'\pm',...
    num2str(round(conf(male80AUC,0.95),2))],['permuted: ',...
    num2str(round(mean(m80P),2)),'\pm',...
    num2str(round(conf(m80P,0.95),2))]})
box off

figure
hold on
plot(squeeze(mean(maleLOOX,[1,2])),squeeze(mean(maleLOOY,[1,2])))
plot(mean(mPX),mean(mPY))
title(['male LOO; p = ',...
    num2str((sum(reshape(mP,1,numel(mP))>mean(maleLOOAUC,[1,2]))+1)/...
    (numel(mP)+1))])
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['male: ',num2str(round(mean(maleLOOAUC,[1,2]),2)),'\pm',...
    num2str(round(conf(reshape(maleLOOAUC,1,numel(maleLOOAUC)),0.95)...
    ,2))],['permuted: ',num2str(round(mean(mP,[1,2]),2)),'\pm',...
    num2str(round(conf(reshape(mP,1,numel(mP)),0.95),2))]})
box off

figure
hold on
plot(mean(mf80X),mean(mf80Y))
plot(mean(fm80X),mean(fm80Y))
plot(mean(f80PX),mean(f80PY))
plot(mean(m80PX),mean(m80PY))
legend({['m->f: ',num2str(round(mean(mf80A),2)),'\pm',...
    num2str(round(conf(mf80A,0.95),3)),' ',...
    num2str((sum(f80P>mean(mf80A))+1)/(numel(f80P)+1))],['f->m: ',...
    num2str(round(mean(fm80A),2)),'\pm',...
    num2str(round(conf(fm80A,0.95),3)),...
    ' ',num2str((sum(m80P>mean(fm80A))+1)/(numel(m80P)+1))],...
    ['permuted female: ',num2str(round(mean(f80P),2)),'\pm',...
    num2str(round(conf(f80P,0.95),2))],...
    ['permuted male: ',num2str(round(mean(m80P),2)),'\pm',...
    num2str(round(conf(m80P,0.95),2))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
title('80:20')
box off

figure
hold on
plot(squeeze(mean(mfX,[1,2])),squeeze(mean(mfY,[1,2])))
plot(squeeze(mean(fmX,[1,2])),squeeze(mean(fmY,[1,2])))
legend({['m->f: ',num2str(round(mean(mfA,[1,2,3]),2)),'\pm',...
    num2str(round(conf(reshape(mfA,1,numel(mfA)),0.95),3))],['f->m: ',...
    num2str(round(mean(fmA,[1,2,3]),2)),'\pm',...
    num2str(round(conf(reshape(fmA,1,numel(mfA)),0.95),3))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
title('LOO')
box off

%%
figure
hold on
% for ii = 1:5
%     plot(squeeze(mean(maleX(:,ii,:),1)),squeeze(mean(maleY(:,ii,:),1)),'k')
%     plot(squeeze(mean(femaleX(:,ii,:),1)),squeeze(mean(femaleY(:,ii,:),1)),'g')
% end
plot(squeeze(mean(mean(maleX))),squeeze(mean(mean(maleY))),'k')
plot(squeeze(mean(mean(femaleX))),squeeze(mean(mean(femaleY))),'g')
legend({['male: ',num2str(round(mean(mean(maleAUC)),2)),'\pm',...
    num2str(round(conf(reshape(maleAUC,1,400),0.95),2))],...
    ['female: ',num2str(round(mean(mean(femaleAUC)),2)),'\pm',...
    num2str(round(conf(reshape(femaleAUC,1,500),0.95),2))]},...
    'location','se')
title('within sex models')
xlabel('FPR'); ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)

figure
hold on
% plot(mean(mfXall,1),mean(mfYall,1))
% plot(mean(fmXall,1),mean(fmYall,1))
% legend({['m->f: ',num2str(round(mean(mfAUC),2)),'\pm',...
%     num2str(round(conf(mfAUC',0.95),3))],...
%     ['f->m: ',num2str(round(mean(fmAUC),2)),'\pm',...
%     num2str(round(conf(fmAUC',0.95),3))]},'location','se')
plot(squeeze(mean(mean(mean(mfXall,1),2),3)),...
    squeeze(mean(mean(mean(mfYall,1),2),3)))
plot(squeeze(mean(mean(mean(fmXall,1),2),3)),...
    squeeze(mean(mean(mean(fmYall,1),2),3)))
legend({['m->f: ',num2str(round(mean(mean(mean(mfAUC))),2)),'\pm',...
    num2str(round(conf(reshape(mfAUC,1,2000),0.95),3))],...
    ['f->m: ',num2str(round(mean(mean(mean(fmAUC))),2)),'\pm',...
    num2str(round(conf(reshape(fmAUC,1,2000),0.95),3))]},'location','se')
title('across sex models')
xlabel('FPR'); ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
%% single features
daFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'})';
maleAsign = squeeze(mean(maleA.*maleSign,[1,2]));
femaleAsign = squeeze(mean(femaleA.*femaleSign,[1,2]));
[~,maleInd] = sort(abs(maleAsign),'descend');
[~,femaleInd] = sort(abs(femaleAsign),'descend');
maleSort = maleAsign(maleInd);
femaleSort = femaleAsign(femaleInd);
maleFeat = daFeat(maleInd);
feamleFeat = daFeat(femaleInd);
figure
hold on
plot(maleAsign(1:24),femaleAsign(1:24),'ok')
plot(maleAsign(25:end),femaleAsign(25:end),'sk')
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
xlabel('male AUROC')
ylabel('female AUROC')
legend({'power','coherence'})
%%
figure
imagesc(mean(maleAUC,2))
colormap viridis
caxis([0.4 1])
title('male->male')
ylabel('male tested')
set(gca,'ytick',1:1:5,'xtick',[])

figure
imagesc(mean(mfAUC,3))
colormap viridis
caxis([0.4 1])
title('male->female')
ylabel('male left out')
xlabel('female tested')
set(gca,'xtick',1:5,'ytick',1:5)

figure
imagesc(mean(femaleAUC,2))
colormap viridis
caxis([0.4 1])
title('female->female')
ylabel('female left out')
set(gca,'xtick',[],'ytick',1:5)

figure
imagesc(mean(fmAUC,3))
colormap viridis
caxis([0.4 1])
title('female->male')
ylabel('female left out')
xlabel('male tested')
set(gca,'xtick',1:5,'ytick',1:5)
%% Pre-Drinking vs. Not Drinking
% files = fileSearch('D:\paper3\preDrinkCombined2\','.mat');
% fStuff = cellfun(@(x) strsplit(x,'_'),files,'UniformOutput',0);
% ids = cellfun(@(x) x{1},fStuff,'UniformOutput',0);
% allPre = []; preID = [];
% allNot = []; notID = [];
% nFiles = numel(files);
% [past,pastID] = deal(cell(nFiles,60));
% for ii = 1:nFiles
%     load(files{ii})
%     if ~isempty(trls{1,1})
%         [b,c,t] = size(psdTrls{1,1}.relPow);
%         thisPow = reshape(psdTrls{1,1}.relPow,b*c,t)';
%         [cmb,b,t] = size(coh{1,1}.rel);
%         thisCoh = reshape(permute(coh{1,1}.rel,[2,1,3]),cmb*b,t)';
%         allPre = [allPre;thisPow,thisCoh]; %#ok
%         preID = [preID;repmat(ids(ii),size(thisPow,1),1)]; %#ok
%         
%         [b,c,t] = size(psdTrls{1,end}.relPow);
%         thisPow = reshape(psdTrls{1,end}.relPow,b*c,t)';
%         [cmb,b,t] = size(coh{1,end}.rel);
%         thisCoh = reshape(permute(coh{1,end}.rel,[2,1,3]),cmb*b,t)';
%         allNot = [allNot;thisPow,thisCoh]; %#ok
%         notID = [notID;repmat(ids(ii),size(thisPow,1),1)]; %#ok
%         if numel(trls) > 2
%             count = 1;
%             for jj = 2:numel(trls)-1
%                 if ~isempty(trls{jj})
%                     [b,c,t] = size(psdTrls{1,jj}.relPow);
%                     thisPow = reshape(psdTrls{1,jj}.relPow,b*c,t)';
%                     [cmb,b,t] = size(coh{1,jj}.rel);
%                     thisCoh = reshape(permute(coh{1,jj}.rel,[2,1,3]),...
%                         cmb*b,t)';
%                 else
%                     thisPow = []; thisCoh = [];
%                 end
%                 past{ii,count} = [thisPow,thisCoh];
%                 pastID{ii,count} = ids(ii);
%                 count = count+1;
%             end
%         end
%     end
% end
% % Also load waterNot data
% load('D:\paper3\analyzed\final\waterAlcoholNotData.mat','watAlc','water'...
%     ,'notWater')
% % allPre = allPre;
% allPreID = preID;
% allNot = [allNot;table2array(water(:,1:60));...
%     table2array(watAlc(watAlc.group==0,1:60));...
%     table2array(notWater(:,1:60))];
% allNotID = [notID;table2array(water(:,61));...
%     table2array(watAlc(watAlc.group==0,61));...
%     table2array(notWater.ID)];
% save('D:\paper3\analyzed\final\preDrinkAllData.mat','allPre',...
%   'allPreID','allNot','allNotID','past','pastID')
%% Find overlapping windows with drinking
% First get subindices for animals with more channels
waterFeat = names({'ILL','CA1L','PL','SL','PR','CA1R','ILR','SR'},...
    {'d','t','a','b','lg','hg'});
alcFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
% Find overlap
waInds = zeros(1,60);
for ii = 1:60
    ind = logicFind(alcFeat{ii},waterFeat,'==');
    if ~isempty(ind)
        waInds(ii) = ind;
    end
end
% Doesn't find PR-SL coherence since in waterFeat it is SL-PR
waInds(43:48) = 157:162;
% Determine overlap between alcFeat and waterAlcFeat (same features,
% different order)
waterAlcFeat = names({'PL','SL','PR','SR'},{'d','t','a','b','lg','hg'});
inds2 = zeros(1,60);
for ii = 1:60
    ind = logicFind(alcFeat{ii},waterAlcFeat,'==');
    if ~isempty(ind)
        inds2(ii) = ind;
    end
end
% Doesn't find PR-SL coherence since in waterFeat it is SL-PR
inds2(43:48) = 43:48;
files = fileSearch('D:\paper3\periDrink\','.mat')';
fStuff = cellfun(@(x) strsplit(x,'_'),files,'UniformOutput',0);
nFiles = numel(files);
watInds = [25:36,52:61,67:69];
for ii = 1:nFiles
    % Get water, alcohol, and not data/samps
    if ~any(ismember(ii,watInds))
       load(['D:\paper3\waterAlcohol\notDrink\',fStuff{ii}{1},'_',fStuff{ii}{2},...
           '_drink_vs_~drink.mat'])
       if ~isempty(trls{1})
           alcSamps{ii} = trls{1}.sampleinfo;
       else
           alcSamps{ii} = [];
       end
       notSamps{ii} = trls{2}.sampleinfo;
       for jj = 1:numel(psdTrls)
           if ~isempty(trls{jj})
               [b,c,t] = size(psdTrls{1,jj}.relPow);
               thisPow = reshape(psdTrls{1,jj}.relPow,b*c,t)';
               [cmb,b,t] = size(coh{1,jj}.normBandCoh);
               thisCoh = reshape(permute(coh{1,jj}.normBandCoh,...
                   [2,1,3]),cmb*b,t)';
               if jj == 1
                  alcData{ii} = [thisPow,thisCoh]; 
               else
                  notData{ii} = [thisPow,thisCoh];
               end
           end
       end
    else
        load(['D:\paper3\waterAlcohol\notDrink\',fStuff{ii}{1},'_',...
            fStuff{ii}{2},'_',fStuff{ii}{3},'_~Both.mat'])
        neitherSamps{ii} = trls{1}.sampleinfo;
        [b,c,t] = size(psdTrls{1,1}.relPow);
        thisPow = reshape(psdTrls{1,1}.relPow,b*c,t)';
        [cmb,b,t] = size(coh{1,1}.normBandCoh);
        thisCoh = reshape(permute(coh{1,1}.normBandCoh,[2,1,3]),cmb*b,t)';
        neitherData{ii} = [thisPow,thisCoh];
        if exist(['D:\paper3\waterAlcohol\processed\',fStuff{ii}{1},'_',...
                fStuff{ii}{2},'_',fStuff{ii}{3},'_water_vs_alcohol.mat'])
            load(['D:\paper3\waterAlcohol\processed\',fStuff{ii}{1},'_',...
                fStuff{ii}{2},'_',fStuff{ii}{3},'_water_vs_alcohol.mat'])
            if ~isempty(trls{1})
                watSamps{ii} = trls{1}.sampleinfo;
                [b,c,t] = size(psdTrls{1,1}.relPow);
                thisPow = reshape(psdTrls{1,1}.relPow,b*c,t)';
                [cmb,b,t] = size(coh{1,1}.normBandCoh);
                thisCoh = reshape(permute(coh{1,1}.normBandCoh,...
                    [2,1,3]),cmb*b,t)';
                theseData = [thisPow,thisCoh];
                theseData = theseData(:,inds2);
                waterData{ii} = theseData;
            else
                waterData{ii} = [];
                watSamps{ii} = [];
            end
            if ~isempty(psdTrls{1,2})
                alcSamps{ii} = trls{1,2}.sampleinfo;
                [b,c,t] = size(psdTrls{1,2}.relPow);
                thisPow = reshape(psdTrls{1,2}.relPow,b*c,t)';
                [cmb,b,t] = size(coh{1,2}.normBandCoh);
                thisCoh = reshape(permute(coh{1,2}.normBandCoh,...
                    [2,1,3]),cmb*b,t)';
                theseData = [thisPow,thisCoh];
                theseData = theseData(:,inds2);
                alcData{ii} = theseData;
            else
                alcData{ii} = [];
                alcSamps{ii} = [];
            end
            % Combine neither and water into not
            notSamps{ii} = [neitherSamps{ii};watSamps{ii}];
            notData{ii} = [neitherData{ii};waterData{ii}];
        else
            notSamps{ii} = neitherSamps{ii};
            notData{ii} = neitherData{ii};
        end
    end
end
% Go through preSamps and remove corresponding notData
% Reset inds and go through files
inds = [];
pastWin = 50;
past = cell(69,pastWin);

allPre = []; allPreID = []; allNot = []; allNotID = [];
for k = 1:nFiles
    load(['D:\paper3\periDrink\',files{k}])
    %     for ii = 1:size(trls{1,31}.sampleinfo,1)
    if any(k == [1:24,37:51,62:66])
        alcInd = 6;
    else
        alcInd = 4;
    end
    % Get alc and water samples for this recordings
    theseAlcSamps = [];
    if ~isempty(alcSamps{k})
        for ii = 1:size(alcSamps{k},1)
            theseAlcSamps = [theseAlcSamps,...
                alcSamps{k}(ii,1):alcSamps{k}(ii,2)];
        end
    end
    theseWaterSamps = [];
    if ~isempty(watSamps{k})
        for ii = 1:size(watSamps{k},1)
            theseWaterSamps = [theseWaterSamps,...
                watSamps{k}(ii,1):watSamps{k}(ii,2)];
        end
    end
    % Go through each possible drink event
    for ii = 1:size(hist.eventTs.t{1,alcInd},1)
        % And go back through time, checking if any data came from that
        % epoch
        for jj = 1:pastWin
            if ~isempty(trls{1,pastWin+1-jj})
                if jj > 1 && overlapPre{k,pastWin+2-jj}(ii) == 1
                    overlapPre{k,pastWin+1-jj}(ii) = 1;
                else
                    % Get index of corresponding pre-window for each drink
                    % event going backwards through time
                    dummy = logicFind(...
                        nearest_idx3(hist.eventTs.t{1,alcInd}(ii),...
                        LFPTs.tvec)-(2000+400*(jj-1)),...
                        trls{1,pastWin+1-jj}.sampleinfo(:,1),'==');
                    if isempty(dummy)
                        inds{k}(ii,jj) = 0;
                    else
                        inds{k}(ii,jj) = dummy;
                        % Get samples for pre window
                        preSamps = trls{pastWin+1-jj}.sampleinfo(...
                            inds{k}(ii,jj),1):...
                            trls{pastWin+1-jj}.sampleinfo(inds{k}(ii,jj),2);
                    end
                    % Check for overlap with notSamps (remove pre-drinking
                    % data from notData)
                    for m = 1:size(notSamps{k},1)
                        if any(ismember(preSamps,notSamps{k}(m,1):...
                                notSamps{k}(m,2)))
                            overlap{k}(m,pastWin+1-jj,ii) = 1;
                        else
                            overlap{k}(m,pastWin+1-jj,ii) = 0;
                        end
                    end
                    % Check for overlap with other alcSamps (remove
                    % pre-drinking data that overlaps with previous
                    % drinking)
                    if any(ismember(preSamps,theseAlcSamps))
                        overlapPre{k,pastWin+1-jj}(ii) = 1;
                    else
                        if any(ismember(preSamps,theseWaterSamps))
                            overlapPre{k,pastWin+1-jj}(ii) = 2;
                        else
                            overlapPre{k,pastWin+1-jj}(ii) = 0;
                        end
                    end
                end
            else
                inds{k}(ii,jj) = 0;
                overlapPre{k,pastWin+1-jj}(ii) = 3;
            end
            % Grab pre-drinking data and store in past
            if overlapPre{k,pastWin+1-jj}(ii) == 0 && inds{k}(ii,jj) ~= 0
                [b,c,~] = size(psdTrls{1,pastWin+1-jj}.relPow);
                thisPow = reshape(psdTrls{1,pastWin+1-jj}.relPow(:,:,...
                    inds{k}(ii,jj)),b*c,1)';
                [cmb,b,~] = size(coh{1,pastWin+1-jj}.normBandCoh);
                thisCoh = reshape(permute(coh{1,pastWin+1-jj}.normBandCoh(:,:,...
                    inds{k}(ii,jj)),[2,1,3]),cmb*b,1)';
                theseData = [thisPow,thisCoh];
                % If a waterAlcohol recording, get rid of extra channels
                % and reorder features
                if alcInd == 4
                theseData = theseData(:,waInds);
                end
                if jj == 1
                    allPre = [allPre;theseData];
                    allPreID = [allPreID;repmat(fStuff{k}(1),...
                        size(theseData,1),1)];
                else
                    past{k,pastWin+1-jj} = [past{k,jj};theseData];
                end
            end
        end
    end
    % Add to pastID
    pastID{k} = fStuff{k}{1};
    % Add not data
    allNot = [allNot;...
        notData{k}(~any(any(overlap{k},3),2),:)];
    allNotID = [allNotID;repmat(fStuff{k}(1),...
        sum(~any(any(overlap{k},3),2)),1)];
end
% save('D:\paper3\analyzed\final\preDrinkAllData2.mat','allNot',...
%     'allNotID','allPre','allPreID','past','pastID')
%% Use everyone; either weight or downsample to 10 sample equivalent
% load('D:\paper3\analyzed\final\preDrinkAllData2.mat')
% Count samples per animal
uP = unique(allPreID);
uN = unique(allNotID);
% pSamps = cellfun(@(x) numel(logicFind(x,allPreID,'==')),uP);
% nSamps = cellfun(@(x) numel(logicFind(x,allNotID,'==')),uN);
samp = 10;
[trainX,trainY,testX,testY,trainWeight,testWeight,a,aP,aSub,allLeft,...
    trainAlcX,trainAlcY,testAlcX,testAlcY,aAlc,aAlcP,alcMdl,...
    trainAlcWeight] = deal([]);
[mdl,mdlNotInds,mdlNotSubInds] = deal(cell(100,19));
drinkShellInds = [7:12,19:24,49:54];
for jj = 1:100
%     disp(jj)
    mdlNotData = []; mdlNotWeight = [];
    mdlNotID = [];
    mdlNotSubData = []; mdlNotSubWeight = [];
    mdlPreData = []; mdlPreWeight = [];
    thisData = []; thisWeight = [];
    % Grab samps from each allNot animal
    for ii = 1:numel(uN)
        all = allNot(logicFind(uN{ii},allNotID,'=='),:);
        if size(all,1) < samp
            thisData = all;
            thisWeight = repmat(samp/size(all,1),size(all,1),1);
            thisNotInd = 1:size(all,1);
            thisID = repmat(uN(ii),size(all,1),1);
        else
            thisNotInd = randperm(size(all,1),samp);
            thisData = all(thisNotInd,:);
            thisWeight = ones(samp,1);
            thisID = repmat(uN(ii),samp,1);
        end
        mdlNotData = [mdlNotData;thisData]; %#ok
        mdlNotWeight = [mdlNotWeight;thisWeight]; %#ok
        mdlNotInds{jj,ii} = thisNotInd';
        mdlNotID = [mdlNotID;thisID]; %#ok
        % If an alcohol only animal, set aside
        if ~isempty(logicFind(uN{ii},uP,'=='))
            mdlNotSubData = [mdlNotSubData;thisData]; %#ok
            mdlNotSubWeight = [mdlNotSubWeight;thisWeight]; %#ok
            test(ii) = 1;
        end
        thisData = []; thisWeight = [];
    end
    % Grab samps from each preNot animal
    for ii = 1:numel(uP)
        all = allPre(logicFind(uP{ii},allPreID,'=='),:);
        if size(all,1) < samp
            thisData = all;
            thisWeight = repmat(samp/size(all,1),size(all,1),1);
        else
            thisData = all(randperm(size(all,1),samp),:);
            thisWeight = ones(samp,1);
        end
        mdlPreData = [mdlPreData;thisData]; %#ok
        mdlPreWeight = [mdlPreWeight;thisWeight]; %#ok
        thisData = []; thisWeight = [];
    end
    %% Combine - using all only alcohol animals (animals with pre data)
    alcData = [mdlPreData;mdlNotSubData];
    alcWeight = [mdlPreWeight;mdlNotSubWeight];
    alcY = [ones(numel(mdlPreWeight),1);zeros(numel(mdlNotSubWeight),1)];
    % Build train and test sets
    [trainAlcInds,~,testAlcInds,~] = trainTest((1:numel(alcWeight))',...
        (1:numel(alcWeight))',0.2);
    trainAlcX(jj,:,:) = alcData(trainAlcInds,drinkShellInds);
    testAlcX(jj,:,:) = alcData(testAlcInds,drinkShellInds);
    trainAlcY(jj,:,:) = alcY(trainAlcInds);
    testAlcY(jj,:) = alcY(testAlcInds);
    trainAlcWeight(jj,:) = alcWeight(trainAlcInds);
    % Build logistic model
    alcMdl{jj} = fitglm(squeeze(trainAlcX(jj,:,:)),trainAlcY(jj,:),...
        'distribution','binomial','binomialSize',size(trainAlcY,2),...
        'weights',trainAlcWeight(jj,:));
    prob = predict(alcMdl{jj},squeeze(testAlcX(jj,:,:)));
    [~,~,~,aAlc(jj)] = perfcurve(testAlcY(jj,:),prob,1);
    prob = predict(alcMdl{jj},squeeze(testAlcX(jj,...
        randperm(size(testAlcY,2),size(testAlcY,2)),:)));
    [~,~,~,aAlcP(jj)] = perfcurve(testAlcY(jj,:),prob,1);
    %% Combine - using all animals in notData
    postData = [mdlPreData;mdlNotData];
    allWeight = [mdlPreWeight;mdlNotWeight];
    allY = [ones(numel(mdlPreWeight),1);zeros(numel(mdlNotWeight),1)];
    % Build train and test sets
    [trainInds,~,testInds,~] = trainTest((1:numel(allWeight))',...
        (1:numel(allWeight))',0.2);
    trainX(jj,:,:) = postData(trainInds,drinkShellInds);
    testX(jj,:,:) = postData(testInds,drinkShellInds);
    trainY(jj,:,:) = allY(trainInds);
    testY(jj,:) = allY(testInds);
    trainWeight(jj,:) = allWeight(trainInds);
    testWeight(jj,:) = allWeight(testInds);
    % Build logistic model
    mdl{jj} = fitglm(squeeze(trainX(jj,:,:)),trainY(jj,:),'distribution'...
        ,'binomial','binomialSize',size(trainY,2),'weights',...
        trainWeight(jj,:));
    prob = predict(mdl{jj},squeeze(testX(jj,:,:)));
    [~,~,~,a(jj)] = perfcurve(testY(jj,:),prob,1);
    prob = predict(mdl{jj},squeeze(testX(jj,randperm(size(testY,2),...
        size(testY,2)),:)));
    [~,~,~,aP(jj)] = perfcurve(testY(jj,:),prob,1);
    %% Repeat model testing, on test set that replaces non-alcohol animal
    % data
    % Find indices in test set that need to be replaced
    rm = zeros(size(testInds));
    for ii = 1:numel(testInds)
       if testInds(ii) > size(mdlPreData,1)
            ind = logicFind(mdlNotID(testInds(ii)-size(mdlPreData,1)),...
                uP,'==');
            if isempty(ind)
               rm(ii) = 1;
            else
                rm(ii) = 0;
            end
       else
           rm(ii) = 0;
       end
    end
    % Collect notPre data that haven't been used yet
   for ii = 1:numel(uP)
        all = allNot(logicFind(uP{ii},allNotID,'=='),:);
        usedInds = mdlNotInds{jj,logicFind(uP{ii},uN,'==')};
        leftInds = ~ismember(1:size(all,1),usedInds);
        allLeft = [allLeft;all(leftInds,:)]; %#ok
   end
   % Replace values in rm with random assortment of allLeft
   testX(jj,logical(rm),:) = allLeft(randperm(size(allLeft,1),sum(rm)),:);
   % Test model
   prob = predict(mdl{jj},squeeze(testX(jj,:,:)));
   [~,~,~,aSub(jj)] = perfcurve(testY(jj,:),prob,1);
end
% save('D:\paper3\analyzed\final\preDrinkModelSubset.mat','trainX','trainY',...
%     'testX','testY','trainWeight','testWeight','a','aP','mdl','uN','uP',...
%     'trainAlcX','trainAlcY','testAlcX','testAlcY','aAlc','aAlcP',...
%     'alcMdl','trainAlcWeight','mdlNotInds','aSub')
%% Compare single feature preDrink to drinkNot
% load('F:\paper3\analyzed\final\preDrinkModel2.mat','trainAlcX',...
%     'trainAlcY','testAlcX','testAlcY','trainAlcWeight')
for ii = 1:100
    for jj = 1:size(trainAlcX,3)
        mdl = fitglm(squeeze(trainAlcX(ii,:,jj)),trainAlcY(ii,:),...
            'distribution','binomial','binomialSize',size(trainAlcY,2),...
            'weights',trainAlcWeight(ii,:));
        prob = predict(mdl,squeeze(testAlcX(ii,:,jj))');
        [~,~,~,aS(ii,jj)] = perfcurve(testAlcY(ii,:),prob,1);
        s(ii,jj) =  sign(table2array(mdl.Coefficients(2,1)));
    end
end
waterAlcFeat = names({'PL','SL','PR','SR'},{'d','t','a','b','lg','hg'});
[sortAS,sortASind] = sort(mean(aS,1),'descend');
waterAlcFeat = waterAlcFeat(sortASind);
sPDM = sign(mean(s,1));
% Get single features of drinkOther
[sF,fA,sD,dA] = deal([]);
for ii = 1:100
    load(['D:/paper3/analyzed/final/drinkOther/drinkOther',...
        num2str(ii),'.mat'],'s','sDA')
    sD(ii,:) = mean(s,1);
    dA(ii,:) = mean(sDA,1);
end
[sortDA,sortDAind] = sort(mean(dA,1),'descend');
sDM = sign(mean(sD,1));
daFeat = names({'PL','SL','PR','SR'},{'d','t','a','b','lg','hg'});
daShellFeat = daFeat(drinkShellInds);
daFeatSort = daShellFeat(sortDAind)';
figure
hold on
plot(mean(aS(:,1:12),1).*sPDM(1:12),mean(dA(:,1:12),1).*sDM(1:12),'o')
plot(mean(aS(:,13:end),1).*sPDM(13:end),...
    mean(dA(:,13:end),1).*sDM(13:end),'s')
plot([-1 1],[0.5 0.5],'k')
plot([-1 1],[-0.5 -0.5],'k')
plot([0.5 0.5],[-1 1],'k')
plot([-0.5 -0.5],[-1 1],'k')
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
xlabel('pre-drink auc')
ylabel('drink auc')
legend({'power','coherence'})
title('pre-drink vs. drink single features')
%% Apply pre models backwards
% Pre-drink
% Load pre-drink data
load('D:/paper3/analyzed/final/preDrinkAllData2.mat')
% Load uP
load('D:/paper3/analyzed/final/preDrinkModelSubset.mat','uN','mdlNotInds')
drinkShellInds = [7:12,19:24,49:54];
% drinkShellInds = 1:60;
% Unpack pastID
% for ii = 1:44
%     for jj = 1:60
%         if ~isempty(pastID{ii,jj})
%             pastID{ii,jj} = pastID{ii,jj}{1,1};
%         end
%     end
% end
% Count pre-samples
% pSamps = cell2mat(cellfun(@(x) size(x,1),past,'uniformoutput',0));
% ids = {'AH11','AH12','AH16','AH1','AH17','AH20','AH4','AH5','AH8'};
ids = unique(allPreID);
samp = 10;
% Use past data and subset of allNot data
% Get number of past windows, subtract one for empty immediate window
pastN = size(past,2)-1;
% Preallocate
[xPast,yPast,xPastP,yPastP] = deal(cell(pastN,100));
[aPast,aPastP] = deal(zeros(pastN,100));
[aPre,aPreP] = deal(zeros(1,100));
for ii = 1:100
    % Load lasso model
%     load(['f:/paper3/analyzed/final/drinkOtherAll/drinkOtherAll',num2str(ii),...
%         '.mat'],'acc','accR')
    load(['D:/paper3/analyzed/final/preDrinkSubset/preDrinkSubset',...
        num2str(ii),'.mat'],'acc','accR')
    % Preallocate
    thisX = cell(1,pastN); thisWeight = cell(1,pastN); 
    thisY = cell(1,pastN); thisXBal = cell(1,pastN); 
    thisYBal = cell(1,pastN);
    for jj = 1:pastN
        for n = 1:numel(ids)
            % Find indices of this animal
            inds{n,jj} = logicFind(ids{n},pastID,'==');
            % Concatenate corresponding past data
            thesePast{n,jj} = cat(1,past{inds{n,jj},jj}); %#ok
            % First build test set using the same procedure as training
            if ~isempty(thesePast{n,jj})
                if size(thesePast{n,jj},1) < samp
                    % Weight
                    thisX{jj} = [thisX{jj};thesePast{n,jj}];
                    thisWeight{jj} = [thisWeight{jj};...
                        repmat(samp/size(thesePast{n,jj},1),...
                        size(thesePast{n,jj},1),1)];
                    thisY{jj} = [thisY{jj};...
                        ones(size(thesePast{n,jj},1),1)];
                else
                    % Downsample
                    thisX{jj} = [thisX{jj};...
                        thesePast{n,jj}(randperm(size(thesePast{n,jj},1)...
                        ,samp),:)];
                    thisWeight{jj} = [thisWeight{jj};...
                        ones(size(thesePast{n,jj},1),1)];
                    thisY{jj} = [thisY{jj};ones(samp,1)];
                end
                % Find corresponding notDrink data that wasn't used in
                % training
                thisInd = logicFind(ids{n},uN,'==');
                theseNot = allNot(strcmp(ids{n},allNotID),:);
                % Remove used inds
                theseNot(mdlNotInds{ii,thisInd},:) = [];
                thisX{jj} = [thisX{jj};...
                    theseNot(randperm(size(theseNot,1),samp),:)];
                thisY{jj} = [thisY{jj};zeros(samp,1)];
                thisWeight{jj} = [thisWeight{jj};ones(samp,1)];
            end
        end
%         prob = predict(mdl{ii},thisX{jj});
        prob = cvglmnetPredict(acc{1}.mdl{1},thisX{jj}(:,drinkShellInds)...
            ,'lambda_1se','response');
        [xPast{jj,ii},yPast{jj,ii},~,aPast(jj,ii)] = ...
            perfcurve(thisY{jj},prob,1);
%         prob = predict(mdl{ii},thisX{jj}(randperm(size(thisX{jj},1),size(thisX{jj},1)),:));
        prob = cvglmnetPredict(accR{1}.mdl{1},thisX{jj}(:,drinkShellInds)...
            ,'lambda_1se','response');
        [xPastP{jj,ii},yPastP{jj,ii},~,aPastP(jj,ii)] = ...
            perfcurve(thisY{jj},prob,1);
        % Concatenated - use all pre and match with not
        thisPast = cat(1,past{:,jj});
        thisCatX = [thisPast;allNot(randperm(size(allNot,1),size(thisPast,1)),:)];
        thisCatY = [ones(size(thisPast,1),1);zeros(size(thisPast,1),1)];
%         prob = predict(mdl{ii},thisCatX);
        prob = cvglmnetPredict(acc{1}.mdl{1},thisCatX(:,drinkShellInds),...
            'lambda_1se','response');
        [~,~,~,aCat(ii,jj)] = perfcurve(thisCatY,prob,1);
    end
    aPre(ii) = acc{1}.acc;
    aPreP(ii) = accR{1}.acc;
    xPre(ii,:) = interp1(linspace(0,1,numel(acc{1}.x)),acc{1}.x,linspace(0,1,50));
    yPre(ii,:) = interp1(linspace(0,1,numel(acc{1}.y)),acc{1}.y,linspace(0,1,50));
    xPreP(ii,:) = interp1(linspace(0,1,numel(accR{1}.x)),accR{1}.x,linspace(0,1,50));
    yPreP(ii,:) = interp1(linspace(0,1,numel(accR{1}.y)),accR{1}.y,linspace(0,1,50));
end
% save('D:/paper3/analyzed/final/pastDrinkLassoSubset.mat','aPre'...
%     ,'aPreP','aPast','xPast','yPast','xPastP','yPastP','xPre','yPre',...
%     'xPreP','yPreP','aPastP')
%% Plot a couple time points
figure
hold on
plot(mean(xPre,1),mean(yPre,1))
thisX = []; thisY = [];
for ii = 1:100
    thisX(ii,:) = interp1(linspace(0,1,numel(xPast{9,ii})),...
        xPast{9,ii},linspace(0,1,50));
    thisY(ii,:) = interp1(linspace(0,1,numel(yPast{9,ii})),...
        yPast{9,ii},linspace(0,1,50));
end
plot(mean(thisX),mean(thisY))
thisX = []; thisY = [];
for ii = 1:100
    thisX(ii,:) = interp1(linspace(0,1,numel(xPastP{9,ii})),...
        xPastP{9,ii},linspace(0,1,50));
    thisY(ii,:) = interp1(linspace(0,1,numel(yPastP{9,ii})),...
        yPastP{9,ii},linspace(0,1,50));
end
plot(mean([thisX;xPreP]),mean([thisY;yPreP]))
legend(['pre: ',num2str(round(mean(aPre),2)),'\pm',...
    num2str(round(conf(aPre,0.95),2)),' ',...
    num2str((sum([aPreP,aPastP(8,:)]>mean(aPre))+1)/(200+1))],...
    ['past: ',num2str(round(mean(aPast(9,:)),2)),'\pm',...
    num2str(round(conf(aPast(9,:),0.95),3)),' ',...
    num2str((sum([aPreP,aPastP(9,:)]>mean(aPast(9,:)))+1)/(200+1))],...
    ['perm: ',num2str(round(mean([aPreP,aPastP(9,:)]),2)),'\pm',...
    num2str(round(conf([aPreP,aPastP(9,:)],0.95),2))])
%% Plot
% load('D:/paper3/analyzed/final/pastDrinkLasso2.mat')
load('D:/paper3/analyzed/final/preDrinkModelSubset.mat','mdl','a','aP','uN',...
    'uP')
preA = a;
% load('D:/paper3/pastDrinkMdls.mat')
figure
hold on
% Flip around past accuracies since they are stored going backwards in time
h = shadedErrorBar(1:26,[mean(aPre);fliplr(mean(aPast,2))],[std(aPre);...
    fliplr(std(aPast,[],2))],'r');
% h = shadedErrorBar(1:26,[mean(aPre),mean(aCat,1)],[std(aPre),...
%     std(aCat,[],1)],'r');
h2 = shadedErrorBar(1:26,[mean(aPreP);fliplr(mean(aPastP,2))],[std(aPreP);...
    fliplr(std(aPastP,[],2))],'k','1');

set(gca,'xtick',(2:5:61),'xticklabel',-2.5:-5:-62.5)
legend([h.mainLine,h2.mainLine],{'balanced','perm'})
xlabel('time before drinking')
ylabel('auc')
%% Single feature preDrink models (just to find best possible single 
% feature for plotting temporal dynamics - expected shitty accurarcy)
load('F:\paper3\analyzed\final\preDrinkModel.mat')
aS = zeros(100,60);
for ii = 1:100
    disp(ii)
    for jj = 1:60
        mdl = fitglm(trainX(ii,:,jj)',trainY(ii,:),...
            'distribution','binomial',...
            'binomialSize',numel(trainY(ii,:)),...
            'weights',trainWeight(ii,:));
        prob = predict(mdl,testX(ii,:,jj)');
        [~,~,~,aS(ii,jj)] = perfcurve(testY(ii,:),prob,1);
    end
end
[sortrMAS,sortInd] = sort(mean(aS,1),'descend'); %#ok
%% Prep data for preBinge model
[data,samp,preFiles] = collateData('F:\paper3\preBingeCombinedNew\',...
    {'dep24';'chow'},{'pow','coh'},'trl','rel');
% Only use animals that exist in both sets
both = [1:3,5:8,10,11];
sub{1,1} = data{1,1}(both,:);
sub{1,2} = data{1,2};
subSamp{1,1} = samp{1,1}(both,:);
subSamp{1,2} = samp{1,2};
% Set up 5 sec preBinge model with notBinge data
preData = sub{1}(:,2:end-1);
pre5 = sub{1}(:,1);
% Grab notBinge data (use notBinge data from sweet that excludes overlap
% with pre, and all notBinge from chow including pre - use notFeed data
% from 24sweetChowNotData to keep data windowing from analysis more
% consistent)
load('F:\paper3\analyzed\final\24sweetChowNotDataNew.mat','notFeed')
notBinge = cell(9,1);
for ii = 1:9
    notBinge{ii,1} = [sub{1}{ii,end};notFeed{ii,2}];
end
% save('F:/paper3/analyzed/final/preBingeAllData.mat','pre5','notBinge',...
%     'preData')
%% Build models
load('D:/paper3/analyzed/final/preBingeAllData.mat')
% Use 10 samples per animal 
samps = 10;
uP = cell2mat(cellfun(@(x) size(x,1),pre5,'uniformoutput',0));
uN = cell2mat(cellfun(@(x) size(x,1),notBinge,'uniformoutput',0));
feedShellInds = [1:12,25:30];
% Preallocate
[allX,allY,allA] = deal([]);
mdlNotInd = cell(12,100);
for ii = 1:100
    mdlPre = []; mdlPreWeight = [];
    for jj = 1:numel(uP)
        thisPre = []; thisWeight = []; %#ok
       if uP(jj) >= samps
          % Downsample
          thisPre = pre5{jj}(randperm(uP(jj),samps),:);
          thisWeight = ones(samps,1);
       else
          % Weight
          thisPre = pre5{jj};
          thisWeight = repmat(samps/uP(jj),uP(jj),1); 
       end
       mdlPre = [mdlPre;thisPre]; %#ok
       mdlPreWeight = [mdlPreWeight;thisWeight]; %#ok
    end
    mdlNot = []; mdlNotWeight = [];
    for jj = 1:numel(uP)
        % Find corresponding animal
        thisNot = []; thisWeight = []; %#ok
       if uN(jj) >= samps
          % Downsample
          thisNotInd = randperm(uN(jj),samps);
          thisNot = notBinge{jj}(thisNotInd,:);
          thisWeight = ones(samps,1);
       else
          % Weight
          thisNotInd = 1:uN{jj};
          thisNot = notBinge{jj};
          thisWeight = repmat(samps/uN(jj),uN(jj),1); 
       end
       mdlNot = [mdlNot;thisNot]; %#ok
       mdlNotWeight = [mdlNotWeight;thisWeight]; %#ok
       mdlNotInd{jj,ii} = thisNotInd;
    end
    % Combine
    data = [mdlPre;mdlNot];
    weight = [mdlPreWeight;mdlNotWeight];
    y = [ones(size(mdlPre,1),1);zeros(size(mdlNot,1),1)];
    [trainInds,~,testInds,~] = trainTest((1:numel(weight))',...
        (1:numel(weight))',0.2);
    trainX(ii,:,:) = data(trainInds,feedShellInds);
    testX(ii,:,:) = data(testInds,feedShellInds);
    trainY(ii,:,:) = y(trainInds);
    testY(ii,:) = y(testInds);
    trainWeight(ii,:) = weight(trainInds);
    testWeight(ii,:) = weight(testInds);
    % Build logistic model
    mdl = fitglm(squeeze(trainX(ii,:,:)),trainY(ii,:),'distribution',...
        'binomial','binomialSize',size(trainY,2),'weights',...
        trainWeight(ii,:));
    prob = predict(mdl,squeeze(testX(ii,:,:)));
    [thisX,thisY,~,allA(ii)] = perfcurve(testY(ii,:),prob,1);
    allX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,size(testY,2)));
    allY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,size(testY,2)));
    prob = predict(mdl,squeeze(testX(ii,randperm(size(testY,2),...
        size(testY,2)),:)));
    [thisX,thisY,~,aP(ii)] = perfcurve(testY(ii,:),prob,1);
    xP(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,...
        size(testY,2))); %#ok
    yP(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,...
        size(testY,2))); %#ok
end
% save('D:\paper3\analyzed\final\preBingeModelDataSubset.mat','trainX','trainY',...
%   'testX','testY','trainWeight','testWeight','allX','allY','allA','xP',...
%   'yP','aP')%,'mdlNotInds')
%% Pre-feed single features
load('D:\paper3\analyzed\final\preBingeModelDataSubset.mat')
singleA = zeros(size(trainX,1),size(trainX,3));
for ii = 1:size(trainX,1)
    for jj = 1:size(trainX,3)
        mdl = fitglm(squeeze(trainX(ii,:,jj)),trainY(ii,:),...
            'distribution','binomial','binomialSize',size(trainY,2),...
            'weights',trainWeight(ii,:));
        pred = predict(mdl,squeeze(testX(ii,:,jj))');
        [~,~,~,singleA(ii,jj)] = perfcurve(testY(ii,:),pred,1);
        signA(ii,jj) = sign(table2array(mdl.Coefficients(2,1)));
    end
end
feedShellInds = [1:12,25:30];
mA = mean(singleA,1);
[msA,sortInd] = sort(mA,'descend');
feat = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
featSub = feat(feedShellInds);
featSort = feat(sortInd);
% save('F:\paper3\analyzed\final\preBingeSingleFeature.mat','mA','msA',...
%     'sortInd','featSort','signA')
%% plot in comparison to feedNot single features
for ii = 1:100
    load(['D:\paper3\analyzed\final\feedOther\feedOther',...
        num2str(ii),'.mat'])
    sF(ii,:) = mean(s,1);
    sA(ii,:) = mean(sFA,1);
end
[sortFA,sortFAind] = sort(mean(sA,1),'descend');
sFM = sign(mean(sF,1));
feat = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
faFeat = featSub(sortFAind)';

figure
hold on
plot(mean(singleA(:,1:12),1).*mean(signA(:,1:12)),...
    mean(sA(:,1:12),1).*sFM(1:12),'o')
plot(mean(singleA(:,13:end),1).*mean(signA(:,13:end)),...
    mean(sA(:,13:end),1).*sFM(13:end),'s')
plot([-1 1],[0.5 0.5],'k')
plot([-1 1],[-0.5 -0.5],'k')
plot([0.5 0.5],[-1 1],'k')
plot([-0.5 -0.5],[-1 1],'k')
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
xlabel('pre-feed auc')
ylabel('feed auc')
legend({'power','coherence'})
title('pre-feed vs. feed single features')
%% Load preFeed models and apply backwards through time
load('D:\paper3\analyzed\final\preBingeModelData.mat','mdlNotInd')
load('D:\paper3\analyzed\final\preBingeAllData.mat')
feedShellInds = [1:12,25:30];
% Pre-allocate
[aPast,aPastP,aFeedOtherPast] = deal(zeros(120,100));
[a,aP,aFeedOther] = deal(zeros(1,100));
samps = 10;
for ii = 1:100
    disp(num2str(ii))
    cd('D:\paper3\analyzed\final\preFeedSubset\')
    load(['preFeedSubset',num2str(ii),'.mat'],'acc','accR','hist')
    a(ii) = acc{1}.acc;
    aP(ii) = accR{1}.acc;
    x(ii,:) = interp1(linspace(0,1,numel(acc{1}.x)),acc{1}.x,linspace(0,1,50));
    y(ii,:) = interp1(linspace(0,1,numel(acc{1}.y)),acc{1}.y,linspace(0,1,50));
    xP(ii,:) = interp1(linspace(0,1,numel(accR{1}.x)),accR{1}.x,linspace(0,1,50));
    yP(ii,:) = interp1(linspace(0,1,numel(accR{1}.y)),accR{1}.y,linspace(0,1,50));
    preMdl = acc{1}.mdl{1};
    preMdlP = accR{1}.mdl{1};
    cd('D:\paper3\analyzed\final\feedOther\')
    load(['feedOther',num2str(ii),'.mat'],'acc')
    feedMdl = acc{1}.mdl{1};
    thisPre = []; thisNot = [];
    for jj = 1:120%size(preData,2)
        for k = 1:size(preData,1)
            % Grab preData
            if size(preData{k,jj},1) < samps
                thisPre = [thisPre;preData{k,jj}]; %#ok
            else
                thisPre = [thisPre;...
                    preData{k,jj}(randperm(size(preData{k,jj},1),samps),...
                    :)]; %#ok
            end
            % Grab not data
            theseNot = notBinge{k};
            % Remove values already used in training
            theseNot(mdlNotInd{k,ii},:) = [];
            thisNot = [thisNot;...
                theseNot(randperm(size(theseNot,1),samps),:)]; %#ok
        end
        % Combine
        testX = [thisPre;thisNot];
        testY = [ones(size(thisPre,1),1);zeros(size(thisNot,1),1)];
        % Test
        pred = cvglmnetPredict(preMdl,testX(:,feedShellInds),'lambda_1se',...
            'response');
        [thisX,thisY,~,aPast(jj,ii)] = perfcurve(testY,pred,1);
        xPast(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,50));
        yPast(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,50));
        % Test using random model
        pred = cvglmnetPredict(preMdlP,testX(:,feedShellInds),'lambda_1se',...
            'response');
        [thisX,thisY,~,aPastP(jj,ii)] = perfcurve(testY,pred,1);
        xPastP(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,50));
        yPastP(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,50));
        % Also test feedNot model
        pred = cvglmnetPredict(feedMdl,testX(:,feedShellInds),'lambda_1se',...
            'response');
        [~,~,~,aFeedOtherPast(jj,ii)] = perfcurve(testY,pred,1);
    end
    pred = cvglmnetPredict(feedMdl,hist.cfg.naive.testX,'lambda_1se',...
            'response'); 
    [~,~,~,aFeedOther(ii)] = perfcurve(hist.cfg.naive.testY,pred,1);
end
% save('D:/paper3/analyzed/final/preFeedPastSubset.mat','a','aP','aPast',...
%     'aPastP','aFeedOther','xPast','yPast','xPastP','yPastP','x','y','xP','yP')
%% Plot couple of time points
figure
hold on
plot(mean(x),mean(y))
plot(squeeze(mean(xPast(:,9,:),1)),squeeze(mean(yPast(:,9,:),1)))
plot(mean([xP;squeeze(xPastP(:,9,:))]),mean([yP;squeeze(yPastP(:,9,:))]))
legend(['pre: ',num2str(round(mean(a),2)),'\pm',...
    num2str(round(conf(a,0.95),2)),' ',...
    num2str((sum([aP,aPast(9,:)]>mean(a))+1)/(200+1))],['past: ',...
    num2str(round(mean(aPast(9,:)),2)),'\pm',...
    num2str(round(conf(aPast(9,:),0.95),2)),' ',...
    num2str((sum([aP,aPastP(9,:)]>mean(aPast(9,:)))+1)/(200+1))],['perm: ',...
    num2str(round(mean([aP,aPast(9,:)]),2)),'\pm',...
    num2str(round(conf([aP,aPast(9,:)],0.95),2))])
%%
figure
hold on
h1 = shadedErrorBar(1:121,[mean(a);mean(aPast,2)],[std(a);...
    std(aPast,[],2)],'k');
% shadedErrorBar(1:121,[mean(aFeedOther);mean(aFeedOtherPast,2)],...
%     [std(aFeedOther);std(aFeedOtherPast,[],2)],'r');
h2 = shadedErrorBar(1:121,[mean(aP);mean(aPastP,2)],[std(aP);...
    std(aPastP,[],2)],'b');
set(gca,'xtick',(2:10:121),'xticklabel',-2.5:-10:-121.5)
xlabel('time before feeding (s)')
ylabel('auc')
legend([h1.mainLine,h2.mainLine],{'real','permuted'})
%% Go through postBinge files and remove trials that overlap with feeding
files = fileSearch('D:/paper3/postFeed/','dep','in');
postData = cell(1,20);
for fi = 1:numel(files)
    load(files{fi})
    % Concatenate all feeding times
    feedSamps = [];
    for ii = 1:size(hist.eventTs.t(7),1)
        % Convert to samples
        startSamp = nearest_idx3(hist.eventTs.t{1,7}(ii),LFPTs.tvec);
        endSamp = nearest_idx3(hist.eventTs.t{1,8}(ii),LFPTs.tvec);
        feedSamps = [feedSamps,startSamp:endSamp]; %#ok
    end
    for ii = 1:20
        if ~isempty(trls{ii})
            % Reshape data
            [b,c,t] = size(psdTrls{ii}.relPow);
            thisPow = reshape(psdTrls{ii}.relPow,b*c,t)';
            [cmb,b,t] = size(coh{ii}.normBandCoh);
            thisCoh = reshape(permute(coh{ii}.normBandCoh,[2,1,3]),...
                cmb*b,t)';
            % Default to keeping all data
            keep = ones(1,size(trls{ii}.sampleinfo,1));
            % Check that none of the post trls overlap with feeding samples
            if ii >= 10
                for k = 1:size(trls{ii}.sampleinfo,1)
                    overlap = ismember(trls{1,ii}.sampleinfo(k,1):...
                        trls{1,ii}.sampleinfo(k,2),feedSamps);
                    if any(overlap)
                        % If any overlap, set that index in keep to 0
                        keep(k) = 0;
                    end
                end
            end
            thisData = [thisPow(logical(keep),:),thisCoh(logical(keep),:)];
            postData{ii} = [postData{ii};thisData];
        end
    end
end
%% Then go through and load preBinge data
files = fileSearch('D:\paper3\preBingeCombinedNew','dep','in');
notFeed = []; preData = cell(1,120);
for ii = 1:size(files,2)
    load(files{ii})
    [b,c,t] = size(psdTrls{241}.relPow);
    notFeedPow = reshape(psdTrls{241}.relPow,b*c,t)';
    [cmb,b,t] = size(coh{241}.rel);
    notFeedCoh = reshape(permute(coh{241}.rel,[2,1,3]),cmb*b,t)';
    notFeed = [notFeed;notFeedPow,notFeedCoh]; %#ok
    for jj = 1:120
        [b,c,t] = size(psdTrls{jj}.relPow);
        thisPow = reshape(psdTrls{jj}.relPow,b*c,t)';
        [cmb,b,t] = size(coh{jj}.rel);
        thisCoh = reshape(permute(coh{jj}.rel,[2,1,3]),cmb*b,t)';
        thisData = [thisPow,thisCoh];
        preData{jj} = [preData{jj};thisData];
    end
end
%% Load feeding data
[feedData,~,~] = collateData('D:\paper3\feed\',{'dep'},{'pow','coh'},...
    'avg','rel');
feedDataCat = cell(1,10);
for ii = 1:10
    feedDataCat{ii} = cat(1,feedData{1}{:,ii});
end
%% Plot
% load('F:\paper3\analyzed\final\preBingeSingleFeature.mat','sortInd')
n = 3;
feat = sortInd(n);
thisPre = cellfun(@(x) mean(x(:,feat)),preData);
thisPreS = cellfun(@(x) std(x(:,feat)),preData);
thisFeed = cellfun(@(x) mean(x(:,feat)),feedDataCat);
thisFeedS = cellfun(@(x) std(x(:,feat)),feedDataCat);
thisPost = cellfun(@(x) mean(x(:,feat)),postData);
thisPostS = cellfun(@(x) std(x(:,feat)),postData);
% Load all binge data
% load('D:/paper3/analyzed/final/24sweetChowNotDataNew.mat','allFeedData')
% feedAverage = mean(cat(1,allFeedData{:,1}),1);
feedAverage = mean(thisFeed);
notAverage = mean(notFeed,1);
figure
hold on
shadedErrorBar(1:60,fliplr(thisPre(1:60)),fliplr(thisPreS(1:60)))
shadedErrorBar(61:64,thisFeed(1:4),thisFeedS(1:4),{'color',[0 0.45 0.74]})
shadedErrorBar(64:70,thisFeed(4:10),thisFeedS(4:10),{'color','b'})
shadedErrorBar(75:78,thisPost(1:4),thisPostS(1:4),{'color','b'})
shadedErrorBar(78:83,thisPost(4:9),thisPostS(4:9),{'color',[0 0.45 0.74]})
shadedErrorBar(84:94,thisPost(10:20),thisPostS(10:20))
plot(1:94,ones(1,94).*notAverage(feat),'k')
plot(1:94,ones(1,94).*feedAverage,'--b')
xlab = [-61.5:1:7.5,0,0,0,0,-6.5:1:12.5];
set(gca,'xtick',1:3:94,'xticklabel',xlab(1:3:94))
title(featSub(sortInd(n)))
%% Drinking feature through time
% Load data
[data,~,files] = collateData('D:\paper3\periDrink\',{'.mat'},{'pow',...
    'coh'},'avg','rel');
% Load an example cfg for eoi data
load(files{1}{1},'hist')
preDrinkInds = 1:26;
preDrinkOverlapInds = 27:30;
drinkStartInds = 31:35;
drinkEndInds = 36:38;
postDrinkOverlapInds = 39:42;
postDrinkInds = 43:50;
% Combine all indices
allInds = {preDrinkInds,preDrinkOverlapInds,drinkStartInds,...
    drinkEndInds,postDrinkOverlapInds,postDrinkInds};
% Try easy way first, use all loaded data (may contain some overlaping
% intervals)
allDataM = cell(1,6);
allDataS = cell(1,6);
for ii = 1:numel(allInds)
    c = 1;
    for jj = allInds{ii}
        allDataM{ii}(c,:) = mean(cat(1,data{1}{:,allInds{ii}(c)}),1);
        allDataS{ii}(c,:) = std(cat(1,data{1}{:,allInds{ii}(c)}),[],1);
        c = c+1;
    end
end
%%
% load('D:\paper3\analyzed\final\preDrinkAllData.mat')
% Set which feature to plot - use 12 for the third 'best' single feature
% from preDrink models which highlights the temporal dynamics
feat = 60;
figure
hold on
shadedErrorBar(1:26,allDataM{1}(:,feat),allDataS{1}(:,feat))
shadedErrorBar(27:30,allDataM{2}(:,feat),allDataS{2}(:,feat),{'color',...
    [0 0.45 0.74]})
shadedErrorBar(31:35,allDataM{3}(:,feat),allDataS{3}(:,feat),{'color','b'})
shadedErrorBar(40:42,allDataM{4}(:,feat),allDataS{4}(:,feat),...
    {'color','b'})
shadedErrorBar(43:46,allDataM{5}(:,feat),allDataS{5}(:,feat),{'color',...
    [0 0.45 0.74]})
shadedErrorBar(47:54,allDataM{6}(:,feat),allDataS{6}(:,feat))
xlab = [-27.5:1:6.5,0,0,0,0,-4.5:1:9.5];
xtickangle(45)
set(gca,'xtick',1:2:54,'xticklabel',xlab(1:2:end));
alcFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
title(alcFeat(feat))
plot(1:54,ones(1,54).*mean(allPre(:,feat)),'k')
plot(1:54,ones(1,54).*mean(table2array(alc(:,feat))),'--b')
%% Build model feed + pre5 vs. not feed
load('F:\paper3\analyzed\final\24sweetChowNotDataNew.mat')
load('F:\paper3\analyzed\final\preBingeAllData.mat')
% Use 30 samps per animal
samps = 30;
[thisFeed,thisNot] = deal(zeros(samps*9,60));
for ii = 1:9
    % Combine sweet with pre
    thisCat = [allFeedData{ii,1};pre5{ii}];
    thisFeed(1+30*(ii-1):30*ii,:) = thisCat(randperm(size(thisCat,1),...
        samps),:);
    thisNotCat = [notBinge{ii};allFeedData{ii,2}];
    thisNot(1+30*(ii-1):30*ii,:) = thisNotCat(...
    randperm(size(thisNotCat,1),samps),:);   
end
a = zeros(1,100);
aS = zeros(100,60);
for ii = 1:100
    [trainX,trainY,testX,testY] = trainTest([thisFeed;thisNot],...
        [ones(1,270);zeros(1,270)],0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial',...
        'binomialSize',numel(trainY));
    prob = predict(mdl,testX);
    [~,~,~,a(ii)] = perfcurve(testY,prob,1);
    for jj = 1:60
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial',...
            'binomialSize',numel(trainY));
        prob = predict(mdl,testX(:,jj));
        [~,~,~,aS(ii,jj)] = perfcurve(testY,prob,1);
    end
end
maS = mean(aS,1);
[sortAS,sortInd] = sort(maS,'descend');
feat = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'}); %#ok
%% OLD FEATURE PLOTTING CODE
% Power or Coh: PreBinge-Binge-PostBinge
% Channels: SL,SR,CL,CR; SLSR,SLCL,SLCR,SRCL,SRCR,CLCR
% Freq: delta,theta,alpha,beta,lgamma,hgamma
clear chan pair freq
chan = [];
pair = 1;
freq = 3;
feat = 'a';
loc = 'SL-SR';
% Get preBinge and notBinge data
files = fileSearch('D:\paper2\preBingeCombined','dep','in');
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
files = fileSearch(['D:\paper2\'...
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
files = fileSearch(['D:\paper2\'...
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
%%
% Plot
figure
hold on
shadedErrorBar(1:61,fliplr(mPre.*100),fliplr(sPre.*100))
shadedErrorBar(61:64,fliplr(mBinge(28:31).*100),fliplr(sBinge(28:31).*100),...
    {'color',[0 0.45 0.74]})
shadedErrorBar(65:92,fliplr(mBinge(1:28).*100),fliplr(sBinge(1:28).*100),...
    {'color','b'})
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
%% Get timeline of alcohol drinking dynamics and feeding dynamics
load('F:\paper3\analyzed\final\bottleTime.mat')
files = fileSearch('F:\paper3\mat\','.mat');
tvec = 0:0.0001:1;
masterDrinkT = zeros(numel(files),numel(tvec));
for ii = 1:numel(files)
    ind = logicFind(1,cellfun(@(x) isequal(files{ii}(1:end-4),x),...
        bottleTime(:,1)),'==');
    botSeconds = bottleTime{ind,2}*60+bottleTime{ind,3};
    load(files{ii},'eventTs','LFPTs')
    total = LFPTs.tvec(end)-botSeconds;
    starts = (eventTs.t{6}-botSeconds)./total;
    ends = (eventTs.t{7}-botSeconds)./total;
    for jj = 1:numel(starts)
        startInd = nearest_idx3(starts(jj),tvec);
        endInd = nearest_idx3(ends(jj),tvec);
        masterDrinkT(ii,startInd:endInd) = 1;
    end
end
files = fileSearch('F:\paper2\preBinge_240sec_redo','Dep');
tvec = 0:0.0001:1;
masterFeedT = zeros(numel(files),numel(tvec));
for ii = 1:numel(files)
    load(files{ii},'hist','LFPTs')
    foodTime = hist.eventTs.t{11};
    total = LFPTs.tvec(end)-foodTime;
    starts = (hist.eventTs.t{7}-foodTime)./total;
    ends = (hist.eventTs.t{8}-foodTime)./total;
    for jj = 1:numel(starts)
        startInd = nearest_idx3(starts(jj),tvec);
        endInd = nearest_idx3(ends(jj),tvec);
        masterFeedT(ii,startInd:endInd) = 1;
    end
end
%% Plot dynamics (both together and alone)
figure
hold on
plot(tvec,smoothdata(mean(masterDrinkT,1)),'b')
plot(tvec,smoothdata(mean(masterFeedT,1)),'r')
xlabel('normalized time (% of session)')
ylabel('% of animals consuming')
box off
title('consumption dynamics')
legend({'drinking','feeding'})

figure
plot(tvec,smoothdata(mean(masterDrinkT,1)),'b')
xlabel('normalized time (% of session)')
ylabel('% of animals drinking')
box off
title('drinking dynamics')
ylim([0 0.6])

figure
plot(tvec,smoothdata(mean(masterFeedT,1)),'r')
xlabel('normalized time (% of session)')
ylabel('% of animals feeding')
box off
title('feeding dynamics')
ylim([0 0.6])
%% Try and compare other sampling (inequal contributions from animals)
% % Use everyone all pre equal not
% load('D:\paper3\preDrinkAllData.mat')
% for ii = 1:100
%     allX = [allPre;allNot(randperm(size(allNot,1),size(allPre,1)),:)];
%     [trainX,trainY,testX,testY] = trainTest(allX,[ones(size(allPre,1),1)...
%         ;zeros(size(allPre,1),1)],0.2);
%     mdl = fitglm(trainX,trainY,'distribution','binomial','binomialSize',...
%         numel(trainY));
%     prob = predict(mdl,testX);
%     [pX(ii,:),pY(ii,:),~,pA(ii)] = perfcurve(testY,prob,1);
% end
% allX = [];
% % Use everyone all pre all not
% for ii = 1:100
%    [trainX,trainY,testX,testY]  = trainTest([allPre;allNot],...
%        [ones(size(allPre,1),1);zeros(size(allNot,1),1)],0.2);
%    mdl = fitglm(trainX,trainY,'distribution','binomial','binomialSize',...
%        numel(trainY));
%    prob = predict(mdl,testX);
%    [allX(ii,:),allY(ii,:),~,allA(ii)] = perfcurve(testY,prob,1);
% end
% % Get even contribution lasso
% cd 'D:/paper3/analyzed/preDrink/'
% for ii = 1:100
%     load(['preDrink',num2str(ii),'.mat'])
%     [thisX,thisY,~,lassoA(ii)] = perfcurve(hist.cfg.naive.testY,...
%         acc{1}.pred,1);
%     if isequal(thisX,thisY)
%        [lassoX(ii,:),lassoY(ii,:)] = deal(linspace(0,1,56));
%     else
%         lassoX(ii,:) = thisX; lassoY(ii,:) = thisY;
%     end
% %     lassoA(ii) = acc{1}.acc;
%     lassoAP(ii) = accR{1}.acc;
% end
% % Get even contribution logistics
% load('D:\paper3\preDrinkModel.mat','x','y','a','xP','yP','aP')
% % Plot
% figure
% hold on
% plot(mean(x,1),mean(y,1))
% plot(mean(pX,1),mean(pY,1))
% plot(mean(allX,1),mean(allY,1))
% plot(mean(lassoX,1),mean(lassoY,1))
% plot(mean(xP,1),mean(yP,1))
% legend({['balanced: ',num2str(round(mean(a),2))],['evenPre: ',...
%     num2str(round(mean(pA),2))],['all: ',num2str(round(mean(allA),2))],...
%     ['lasso balanced: ',num2str(round(mean(lassoA),2))],['perm: ',...
%     num2str(round(mean(aP),2))]})
% title('pre-drink log mdl comparisons')
% ylabel('TPR'); xlabel('FPR')
% % save('D:/paper3/preDrinkMdls.mat','x','y','a','xP','yP','aP','pX',...
% %   'pY','pA','allX','allY','allA','lassoX','lassoY','lassoA')