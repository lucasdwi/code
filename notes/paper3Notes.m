%% Grab all water and alcohol data
% Cross compare feature lists
waterFeat = names({'ILL','CA1L','PL','SL','PR','CA1R','ILR','SR'},...
    {'d','t','a','b','lg','hg'});
alcFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
% Find overlap
for ii = 1:60
    ind = logicFind(alcFeat{ii},waterFeat,'==');
    if ~isempty(ind)
        inds(ii) = ind;
    end
end
% Doesn't find PR-SL coherence since in waterFeat it is SL-PR
inds(43:48) = 157:162;
% Determine overlap between alcFeat and waterAlcFeat (same features,
% differnet order)
waterAlcFeat = names({'PL','SL','PR','SR'},{'d','t','a','b','lg','hg'});
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
[wData,samp,wFiles] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\waterProcessed\'],{'.mat'},{'pow','coh'},'trl','rel');
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),wData{1,1}(:,1),'uniformoutput',0));
animal = []; day = [];
for ii = 1:size(wData{1},1)
    parts = strsplit(wFiles{1}{ii},'_');
    animal{ii,1} = repmat({['AH',parts{2}]},n(ii,1),1);
    day{ii,1} = repmat({parts{3}},n(ii,1),1);
end
% Concatenate and z-score water data
allWater = zscore(cat(1,wData{1}{:,1}));
% Limit and tabulate water data
water = [array2table(allWater(:,inds)),cat(1,animal{:}),cat(1,day{:})];
% Add groups (1 = alcohol; 0 = water)
water(:,63) = array2table(ones(size(water,1),1));
water.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];
%__________________________________________________________________________
% Grab alcohol only data
[aData,samp,aFiles] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\drinkNot\'],{'.mat'},{'pow','coh'},'trl','rel');
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),aData{1,1}(:,1),'uniformoutput',0));
animal = []; day = [];
for ii = 1:size(aData{1},1)
    parts = strsplit(aFiles{1}{ii},'_');
    animal{ii,1} = repmat({parts{1}},n(ii,1),1);
    day{ii,1} = repmat({parts{2}},n(ii,1),1);
end
% Concenate and z-score alcohol data
allAlc = zscore(cat(1,aData{1}{:,1}));
% Limit and tabulate alcohol data
alc = [array2table(allAlc(:,1:60)),cat(1,animal{:}),cat(1,day{:})];
% Add groups (1 = alcohol; 0 = water)
alc(:,63) = array2table(ones(size(alc,1),1));
alc.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];
%__________________________________________________________________________
% Grab water/alcohol data
[awData,samp,awFiles] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\waterAlcohol\processed\'],{'.mat'},{'pow','coh'},'trl'...
    ,'rel');
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),awData{1,1},'uniformoutput',0));
animal = []; day = [];
for ii = 1:size(awData{1},1)
    parts = strsplit(awFiles{1}{ii},'_');
    for jj = 1:2
        animal{ii,jj} = repmat({['AH',parts{2}]},n(ii,jj),1);
        day{ii,jj} = repmat({parts{4}},n(ii,jj),1);
    end
end
% Concenate and z-score water/alcohol data
allWatAlc = [zscore(cat(1,awData{1}{:,1}));zscore(cat(1,awData{1}{:,2}))];
% Limit and tabulate water/alcohol data
watAlc = [array2table(allWatAlc(:,inds2)),[cat(1,animal{:,1});...
    cat(1,animal{:,2})],[cat(1,day{:,1});cat(1,day{:,2})]];
% Add groups (1 = alcohol; 0 = water)
watAlc(:,63) = array2table([zeros(sum(n(:,1)),1);ones(sum(n(:,2)),1)]);

watAlc.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];
%__________________________________________________________________________
% Combine all tables
allData = [water;alc;watAlc];
%% Get information on data
% Count number of trials for each group from each animal
ids = unique(allData.ID);
for ii = 1:size(ids,1)
    nSamps(ii,1) = sum(allData.group(strcmp(allData.ID,ids{ii}))==1);
    nSamps(ii,2) = sum(allData.group(strcmp(allData.ID,ids{ii}))==0);
end
% Go through and find the animals with at least 1 samp in each group
for ii = 1:size(nSamps)
    both(ii) = nSamps(ii,1)>0 && nSamps(ii,2)>0;
end
bothIDs = logicFind(1,both,'==');
minSamp = min(nSamps(bothIDs,:),[],2);
% Subset data to include only animals with >=5 samples
subData = [];
for ii = 1:size(bothIDs,2)
    if minSamp(ii)>5
    thisWat = allData(strcmp(allData.ID,ids(bothIDs(ii))) & ...
        allData.group==0,:);
    thisAlc = allData(strcmp(allData.ID,ids(bothIDs(ii))) & ...
        allData.group==1,:);
    subData = [subData;thisWat(randperm(size(thisWat,1),minSamp(ii)),:);...
        thisAlc(randperm(size(thisAlc,1),minSamp(ii)),:)];
    end
end
%% Using allData build full models
% Determine 80/20 split
nTrain = round(size(allData,1)*.8); nTest = size(allData,1)-nTrain;
% Set up formula to be used in linear mixed model
formula = [];
for ii = 1:60
    formula = [formula,watAlc.Properties.VariableNames{ii},'+'];
end
formula = ['group~',formula,'(1|ID)'];
% Build linear mixed model with ID as random effect
for ii = 1:100
    testInds = randperm(size(allData,1),nTest);
    trainInds = ~ismember(1:size(allData,1),testInds);
    mdl = fitglme(allData(trainInds,:),formula,'Distribution','binomial'...
        ,'binomialsize',nTrain);
    prob = predict(mdl,allData(testInds,:));
    [glmeFPR(ii,:),glmeTPR(ii,:),~,glmeA(ii)] = perfcurve(...
        allData.group(testInds),prob,1);
    [glmePrec(ii,:),glmeRecall(ii,:),~,glmeA2(ii)] = perfcurve(...
        allData.group(testInds),prob,1,'XCrit','prec','YCrit','reca');
    
    mdl = fitglm(allData(trainInds,:),formula(1:end-7),'Distribution','binomial'...
        ,'binomialsize',nTrain);
    prob = predict(mdl,allData(testInds,:));
    [glmFPR(ii,:),glmTPR(ii,:),~,glmA(ii)] = perfcurve(...
        allData.group(testInds),prob,1);
end
%% Just use animals with data from both groups (subData)
nTrain = round(size(subData,1)*.8); nTest = size(subData,1)-nTrain;
for ii = 1:100
    testInds = randperm(size(subData,1),nTest);
    trainInds = ~ismember(1:size(subData,1),testInds);
    
    mdl = fitglme(subData(trainInds,:),formula,'Distribution','binomial'...
        ,'binomialsize',nTrain);
    prob = predict(mdl,subData(testInds,:));
    [glmeFPR(ii,:),glmeTPR(ii,:),~,glmeA(ii)] = perfcurve(...
        subData.group(testInds),prob,1);
    [glmePrec(ii,:),glmeRecall(ii,:),~,glmeA2(ii)] = perfcurve(...
        subData.group(testInds),prob,1,'XCrit','prec','YCrit','reca');
    
    mdl = fitglm(subData(trainInds,:),formula(1:end-7),'Distribution','binomial'...
        ,'binomialsize',nTrain);
    prob = predict(mdl,subData(testInds,:));
    [glmFPR(ii,:),glmTPR(ii,:),~,glmA(ii)] = perfcurve(...
        subData.group(testInds),prob,1);
    
    mdl = fitglme(subData(trainInds,:),formula,'Distribution','binomial'...
        ,'binomialsize',nTrain);
    prob = predict(mdl,subData(testInds(randperm(nTest,nTest)),:));
    [glmerFPR(ii,:),glmerTPR(ii,:),~,glmerA(ii)] = perfcurve(...
        subData.group(testInds),prob,1);
end
%% Plot
figure
hold on
plot(mean(glmeFPR,1),mean(glmeTPR,1))
plot(mean(glmerFPR,1),mean(glmerTPR,1),'--k')
legend({['GLME: ',num2str(round(mean(glmeA),2)),'\pm',...
    num2str(round(conf(glmeA,0.95),2))],['GLMEp: ',...
    num2str(round(mean(glmerA),2)),'\pm',...
    num2str(round(conf(glmerA,0.95),2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1);
xlabel('FPR'); ylabel('TPR')
title('Water vs. Alcohol')
%% Run single feature models
for ii = 1:10
    disp(num2str(ii))
    testInds = randperm(size(subData,1),nTest);
    trainInds = ~ismember(1:size(subData,1),testInds);
    for jj = 1:size(subData,2)-3
%         mdl = fitglme(allData(trainInds,:),['Var218~Var',num2str(jj),...
%             '+(1|Var217)'],'distribution','binomial','binomialsize',nTrain);
%         mdl = fitglme(allData(trainInds,:),['Var62~Var',num2str(jj),...
%             '+(1|Var61)'],'distribution','binomial','binomialsize',nTrain);
%         prob = predict(mdl,allData(testInds,:));
        mdl = fitglme(subData(trainInds,:),...
            ['group~',alcFeat{jj},'+(1|ID)'],'distribution','binomial',...
            'binomialsize',nTrain);
        prob = predict(mdl,subData(testInds,:));
%         mdl = fitglm(table2array(allData(trainInds,jj)),table2array(allData(trainInds,end)),'distribution','binomial','binomialsize',nTrain);
%         prob = predict(mdl,table2array(allData(testInds,jj)));
%         [glmeLogFPR(ii,jj,:),glmeLogTPR(ii,jj,:),~,glmeLogA(ii,jj)] = ...
%             perfcurve(allData.Var218(testInds),prob,1);
%         [logFPR(ii,jj,:),logTPR(ii,jj,:),~,logA(ii,jj)] = ...
%             perfcurve(allData.Var218(testInds),prob,1);
        [glmeLogFPR(ii,jj,:),glmeLogTPR(ii,jj,:),~,glmeLogA(ii,jj)] = ...
            perfcurve(subData.group(testInds),prob,1,'TVals',0:1/100:1,...
            'UseNearest',0);
    end
end
%% Grab drinkNot data (taken from waterAlcohol cohort but set to look at 
% data corresponding to neither behavior)
[dataNot,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab'...
    '\data\paper3\waterAlcohol\notDrink\'],{'.mat'},{'pow','coh'},'trl',...
    'rel');
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),dataNot{1,1}(:,1),'uniformoutput',0));
animal = []; day = [];
for ii = 1:size(dataNot{1},1)
    parts = strsplit(files{1}{ii},'_');
    animal{ii,1} = repmat({[parts{1},parts{2}]},n(ii,1),1);
    day{ii,1} = repmat({parts{4}},n(ii,1),1);
end
% Concenate and z-score alcohol data
allDataNot = zscore(cat(1,dataNot{1}{:,1}));
% Limit and tabulate alcohol data
notD = [array2table(allDataNot(:,1:60)),cat(1,animal{:}),cat(1,day{:})];
% Add groups (0 = not)
notD(:,63) = array2table(zeros(size(notD,1),1));
notD.Properties.VariableNames = [alcFeat(1:60),'ID','day','group'];
%% Build models going between and within behavior
for jj = 1:100
    disp(jj)
    allAlc = []; allWater = [];
    for ii = 1:size(bothIDs,2)
        if minSamp(ii)>5
            thisWat = allData(strcmp(allData.ID,ids(bothIDs(ii))) & ...
                allData.group==0,:);
            thisWat.group(1:size(thisWat,1)) = ones(size(thisWat,1),1);
            thisNot = notD(strcmp(notD.ID,ids(bothIDs(ii))),:);
            % Random subset of thisNot
            allWater = [allWater;thisWat(randperm(size(thisWat,1),minSamp(ii)),:);...
                thisNot(randperm(size(thisNot,1),minSamp(ii)),:)];
            
            thisAlc = allData(strcmp(allData.ID,ids(bothIDs(ii))) & ...
                allData.group==1,:);
            allAlc = [allAlc;thisAlc(randperm(size(thisAlc,1),minSamp(ii)),:);...
                thisNot(randperm(size(thisNot,1),minSamp(ii)),:)];
        end
    end
    nTrain = round(size(allAlc,1)*0.8); nTest = size(allAlc,1)-nTrain;
    trainInds = randperm(size(allWater,1),nTrain);
    testInds = ~ismember(1:size(allWater,1),trainInds);
    % Wat -> wat
    mdl = fitglm(table2array(allWater(trainInds,1:60)),...
        table2array(allWater(trainInds,63)),'distribution','binomial',...
        'binomialsize',nTrain);   
    prob = predict(mdl,table2array(allWater(testInds,1:60)));
    [wwX(jj,:),wwY(jj,:),~,wwA(jj)] = perfcurve(...
        allWater.group(testInds),prob,1);%,'TVals',0:1/100:1,'UseNearest',0);
    % Alc -> alc
    mdl = fitglm(table2array(allAlc(trainInds,1:60)),...
        table2array(allAlc(trainInds,63)),'distribution','binomial',...
        'binomialsize',nTrain);
    prob = predict(mdl,table2array(allAlc(testInds,1:60)));
    [aaX(jj,:),aaY(jj,:),~,aaA(jj)] = perfcurve(...
        allAlc.group(testInds),prob,1);%,'TVals',0:1/100:1,'UseNearest',0);
    % Wat -> alc
    mdl = fitglm(table2array(allWater(trainInds,1:60)),...
        table2array(allWater(trainInds,63)),'distribution','binomial',...
        'binomialsize',nTrain);
    prob = predict(mdl,table2array(allAlc(testInds,1:60)));
    [waX(jj,:),waY(jj,:),~,waA(jj)] = perfcurve(...
        allAlc.group(testInds),prob,1);%,'TVals',0:1/100:1,'UseNearest',0);
    % Alc -> wat
    mdl = fitglm(table2array(allAlc(trainInds,1:60)),...
        table2array(allAlc(trainInds,63)),'distribution','binomial',...
        'binomialsize',nTrain);
    prob = predict(mdl,table2array(allWater(testInds,1:60)));
    [awX(jj,:),awY(jj,:),~,awA(jj)] = perfcurve(...
        allWater.group(testInds),prob,1);%,'TVals',0:1/100:1,'UseNearest',0);
    % Alc -> rand alc
    mdl = fitglm(table2array(allAlc(trainInds,1:60)),...
        table2array(allAlc(trainInds,63)),'distribution','binomial',...
        'binomialsize',nTrain);
    randTestAlc = table2array(allAlc(testInds,1:60));
    randTestAlc = randTestAlc(randperm(size(randTestAlc,1),...
        size(randTestAlc,1)),:);
    prob = predict(mdl,randTestAlc);
    [aarX(jj,:),aarY(jj,:),~,aarA(jj)] = perfcurve(...
        allAlc.group(testInds),prob,1);%,'TVals',0:1/100:1,'UseNearest',0);
end
%% Plot performances
figure
hold on
plot(mean(wwX,1),mean(wwY,1))
plot(mean(aaX,1),mean(aaY,1))
plot(mean(waX,1),mean(waY,1))
plot(mean(awX,1),mean(awY,1))
plot(mean(aarX,1),mean(aarY,1),'--k')

mWWA = round(mean(wwA),2);
mAAA = round(mean(aaA),2);
mWAA = round(mean(waA),2);
mAWA = round(mean(awA),2);
mAARA = round(mean(aarA),2);

cWWA = round(conf(wwA,0.95),2);
cAAA = round(conf(aaA,0.95),2);
cWAA = round(conf(waA,0.95),2);
cAWA = round(conf(awA,0.95),2);
cAARA = round(conf(aarA,0.95),2);

legend({['W->W: ',num2str(mWWA),'\pm',num2str(cWWA)],...
    ['A->A: ',num2str(mAAA),'\pm',num2str(cAAA)],...
    ['W->A: ',num2str(mWAA),'\pm',num2str(cWAA)],...
    ['A->W: ',num2str(mAWA),'\pm',num2str(cAWA)],...
    ['A->P: ',num2str(mAARA),'\pm',num2str(cAARA)]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1);
xlabel('FPR'); ylabel('TPR')
title('Water vs. Alcohol')
%% Predict drinking amounts from LFPs during drinking session; use average 
% of both drinking and notDrinking
[data,samps,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\drinkNot\'],{'.mat'},{'pow','coh'},'trl','rel');
for ii = 1:size(data{1,1},1)
    avgData(ii,:) = mean(cat(1,data{1,1}{ii,:}),1);
end
%%
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\drinkAmountData.mat')
mini = min(cell2mat(drink(:,1)));
maxi = max(cell2mat(drink(:,1)));
bins = mini:0.1:maxi;
for ii = 1:14
    if ii~= 14
        n(ii,:) = cell2mat(drink(:,1))>=bins(ii) & cell2mat(drink(:,1))<bins(ii+1); 
    else
        n(ii,:) = cell2mat(drink(:,1))>=bins(ii); 
    end
end
sumN = sum(n,2);
w = 1./sumN;

weights = n.*w;
weights(weights==0) = [];
save('C:\Users\Pythia\Documents\GreenLab\data\paper3\drinkAmountData.mat','avgData','drink','weights')
%% Run single feature models
nTest = 9;
for ii = 41%1:60
    for jj = 1:100
        testInds = randperm(size(avgData,1),nTest);
        trainInds = ~ismember(1:size(avgData,1),testInds);
        mdl = fitglm(avgData(trainInds,ii),cell2mat(drink(trainInds,1)));
        c(jj) = table2array(mdl.Coefficients(2,1));
        prob = predict(mdl,avgData(testInds,ii));
        err(ii,jj,:) = cell2mat(drink(testInds,1))-prob;
        mse(ii,jj) = mean(cell2mat(drink(testInds,1))-prob)^2;
        r(ii,jj) = mdl.Rsquared.Ordinary;
    end
end
%%
[sMSE,sInd] = sort(abs(mean(mse,2)),'ascend');
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
sName = nameVect(sInd)';
%%
for ii = 1:100
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\drinkingAmountsWeighted\',num2str(ii),'.mat'],'accArray','accArrayRand')
   a(ii,:) = accArray{1}.acc;
   aR(ii,:) = accArrayRand{1}.acc;
end
doubleHist(reshape(a,1,900),reshape(aR,1,900),'xlab','Error (g/kg)','main','Drinking Amount');
%% Only keep NAc shell and PL channels
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper3\waterAlcohol\mat\','.mat');
for ii = 1:size(files,2)
    disp(num2str(ii))
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\waterAlcohol\mat\',files{ii}])
    LFPTs.data = LFPTs.data([3,4,5,8],:);
    LFPTs.label = LFPTs.label([3,4,5,8]);
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\waterAlcohol\NAcPL\',files{ii}],'LFPTs','eventTs','pl2','adfreq')
end
%% Grab water/alcohol data
% All channels
% [data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
%     'data\paper3\waterAlcohol\processedOld\'],{'.mat'},{'pow','coh'},'trl',...
%     'rel');
% Just NAc Shell and PL
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\waterAlcohol\processed\'],{'.mat'},{'pow','coh'},'trl',...
    'rel');
oldData = data;
%% Create animal ID assignment
% First limit animal 5 to only 10 samples
data{1,1}{5,1} = data{1,1}{5,1}(1:10,:);
n = cell2mat(cellfun(@(x) size(x,1),data{1,1},'uniformoutput',0));
for ii = 1:size(data{1},1)
    animal{ii,1} = repmat({files{1,1}{ii}(4:5)},n(ii,1),1);
    animal{ii,2} = repmat({files{1,1}{ii}(4:5)},n(ii,2),1);
end
%% Collate data together into water and alcohol tables for ME models
water = array2table(zscore(cat(1,data{1}{:,1})));
water(:,end+1) = cell2table(cat(1,animal{:,1}));
water(:,end+1) = array2table(zeros(size(water,1),1));

alc = array2table(zscore(cat(1,data{1}{:,2})));
alc(:,end+1) = cell2table(cat(1,animal{:,2}));
alc(:,end+1) = array2table(ones(size(alc,1),1));

allData = [water;alc];
%% Try only data from animals with both water and alcohol
inds = [];
for ii = 1:24
    if ~isempty(oldData{1,1}{ii,1}) && ~isempty(oldData{1,1}{ii,2})
        data(ii,:) = oldData{1,1}(ii,:);
        inds = [inds,ii];
    end
end
n = cell2mat(cellfun(@(x) size(x,1),data,'uniformoutput',0));
for ii = 1:size(data,1)
    animal{ii,1} = repmat({files{1,1}{ii}(4:5)},n(ii,1),1);
    animal{ii,2} = repmat({files{1,1}{ii}(4:5)},n(ii,2),1);
end
% % Only keep 20 water of animal 4
% data{4,1} = data{4,1}(1:20,:);
% animal{4,1} = animal{4,1}(1:20,:);
% Only keep 20 water of animal 5
data{5,1} = data{5,1}(1:20,:);
animal{5,1} = animal{5,1}(1:20,:);
water = array2table(zscore(cat(1,data{inds,1})));
water(:,end+1) = cell2table(cat(1,animal{inds,1}));
water(:,end+1) = array2table(zeros(size(water,1),1));

alc = array2table(zscore(cat(1,data{inds,2})));
alc(:,end+1) = cell2table(cat(1,animal{inds,2}));
alc(:,end+1) = array2table(ones(size(alc,1),1));

allData = [water;alc];
%% Build linear mixed effects models with 20% samples withheld
nTest = round(size(allData,1)*.2); nTrain = size(allData,1)-nTest;
formula = [];
for ii = 1:size(allData,2)-2
   formula = [formula,'Var',num2str(ii),'+']; 
end
for ii = 1:100
    testInds = randperm(size(allData,1),nTest);
    trainInds = ~ismember(1:size(allData,1),testInds);
    
%     mdl = fitglme(allData(trainInds,:),['Var218~',formula,'(1|Var217)'],...
%         'Distribution','binomial','binomialsize',nTrain);
    mdl = fitglme(allData(trainInds,:),['Var62~',formula,'(1|Var61)'],...
        'Distribution','binomial','binomialsize',nTrain);
    prob = predict(mdl,allData(testInds,:));
%     [glmeFPR(ii,:),glmeTPR(ii,:),~,glmeA(ii)] = perfcurve(...
%         allData.Var218(testInds),prob,1);
%     [glmePrec(ii,:),glmeRecall(ii,:),~,glmeA2(ii)] = perfcurve(...
%         allData.Var218(testInds),prob,1,'XCrit','prec','YCrit','reca');
    [glmeFPR(ii,:),glmeTPR(ii,:),~,glmeA(ii)] = perfcurve(...
        allData.Var62(testInds),prob,1);
    [glmePrec(ii,:),glmeRecall(ii,:),~,glmeA2(ii)] = perfcurve(...
        allData.Var62(testInds),prob,1,'XCrit','prec','YCrit','reca');
end
figure
plot(mean(glmeFPR,1),mean(glmeTPR,1))
figure
plot(mean(glmeRecall,1),mean(glmePrec,1))
%% Build single feature models (linear mixed effects)
for ii = 1:10
    disp(num2str(ii))
    testInds = randperm(size(allData,1),nTest);
    trainInds = ~ismember(1:size(allData,1),testInds);
    for jj = 1:size(allData,2)-2
%         mdl = fitglme(allData(trainInds,:),['Var218~Var',num2str(jj),...
%             '+(1|Var217)'],'distribution','binomial','binomialsize',nTrain);
%         mdl = fitglme(allData(trainInds,:),['Var62~Var',num2str(jj),...
%             '+(1|Var61)'],'distribution','binomial','binomialsize',nTrain);
%         prob = predict(mdl,allData(testInds,:));
        mdl = fitglm(table2array(allData(trainInds,jj)),table2array(allData(trainInds,end)),'distribution','binomial','binomialsize',nTrain);
        prob = predict(mdl,table2array(allData(testInds,jj)));
%         [glmeLogFPR(ii,jj,:),glmeLogTPR(ii,jj,:),~,glmeLogA(ii,jj)] = ...
%             perfcurve(allData.Var218(testInds),prob,1);
%         [logFPR(ii,jj,:),logTPR(ii,jj,:),~,logA(ii,jj)] = ...
%             perfcurve(allData.Var218(testInds),prob,1);
        [glmeLogFPR(ii,jj,:),glmeLogTPR(ii,jj,:),~,glmeLogA(ii,jj)] = ...
            perfcurve(allData.Var62(testInds),prob,1,'TVals',0:1/100:1,'UseNearest',0);
    end
end
%% Sort mean logistic AUCs
glmeLogAM = mean(glmeLogA,1);
%% Grab data for drinkNot
[data,~,~] = collateData(['C:\Users\Pythia\Documents\GreenLab\data\'...
    'paper3\drinkNot\'],{'.mat'},{'pow','coh'},'trl','rel');
% Concat all together to avoid issues of low drinking samples 
catData{1,1} = cat(1,data{1,1}{:,1});
catData{1,2} = cat(1,data{1,1}{:,2});
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
%     'drinkNotRaw.mat'])
%% Run evenDataSplit: use 50-50 training, but force train/test into 80/20
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drinkNotRaw.mat'])
% Use all data
[all,~,rnd,~] = evenDataSplit(catData,18494,4624,'ADA',100); %#ok
save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'concatData100.mat'],'all','rnd')
% Or foce into same sample size as bingeNot data
[all,~,rnd,~] = evenDataSplit(catData,6000,1200,'ADA',100);
save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000_100.mat'],'all','rnd')
%% Build drinkNot models using drink6000.mat for comparison to bingeNot
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000_100.mat'])
% load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
%     'concatData.mat'])
% Preallocate
[x,y,xRand,yRand] = deal(cell(1,100));
[a,aRand] = deal(zeros(1,100));
for n = 1:100
    disp([num2str(n),' of 100'])
    % Set up training data
    trainX = all.trainX{n};
    trainY = all.trainY{n};
    % Set up testing data
    testX = all.testX{n};
    testY = all.testY{n};
    % Build and test model on concat data
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [x{n},y{n},~,a(n)] = perfcurve(testY,prob,1,'TVals',0:1/4623:1,...
        'UseNearest',0);
    % Store outcomes
    concatData.model{n}  = mdl;
    % Set up testing data
    trainY = rnd.allTrainY{n};
    testY = rnd.allTestY{n};
    % Build and test model on concat data
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [xRand{n},yRand{n},~,aRand(n)] = perfcurve(testY,prob,1,'Tvals',...
        0:1/4623:1,'UseNearest',0);
end
concatData.rocX = x;
concatData.rocY = y;
concatData.auc = a;
% Calculate average
mX = mean(cat(2,x{:}),2);
mY = mean(cat(2,y{:}),2);
mRandX = mean(cat(2,xRand{:}),2);
mRandY = mean(cat(2,yRand{:}),2);
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
%     'drinkNotModel100.mat'],'concatData','aRand','xRand','yRand')
%% Plot
figure
hold on
for ii = 1:20
    p = plot(x{ii},y{ii});
    col = get(p,'Color');
    hsv = rgb2hsv(col);
    newCol = hsv2rgb(hsv-[0 0.5 0]);
    set(p,'Color',newCol);
    plot(xRand{ii},yRand{ii},'--','Color',[0.4 0.4 0.4])
end
h(1) = plot(mX,mY,'-k');
h(2) = plot(mRandX,mRandY,'--k');
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
ylabel('True Positive Rate')
xlabel('False Positive Rate')
legend(h,{['Real: ',num2str(round(mean(a),2))],...
    ['Permuted: ',num2str(round(mean(aRand),2))]},'Location','se')
title('Drink vs. Other')
%% drinkNot monads (for Angela)
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000.mat'])
for ii = 1:20
    trainX = all.trainX{ii};
    trainY = all.trainY{ii};
    testX = all.testX{ii};
    testY = all.testY{ii};
    for jj = 1:60
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        prob = predict(mdl,testX(:,jj));
        [x{ii,jj},y{ii,jj},~,a(ii,jj)] = perfcurve(testY,prob,1,...
            'TVals',0:1/4623:1,'UseNearest',0);
        dir(ii,jj) = table2array(mdl.Coefficients(2,1));
    end
end
mDir = exp(mean(dir,1)')>1;
mA = mean(a,1)';
[smA,inds] = sort(mA,'descend');
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
nameSort = nameVect(inds)';
%% Build drinkNot models split by sex and tested across sex with leave one 
% animal out testing
[data,~,~] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\drinkNot\'],{'.mat'},{'pow','coh'},'trl','rel');
% Split
male = data{1}([1:18,24:25],:);
female = data{1}([21:23,26:45],:);
% Set up animal number vector
animal{1} = [1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5];
animal{2} = [1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5];
save('C:\Users\Pythia\Documents\GreenLab\data\paper3\maleFemale.mat',...
    'male','female','animal')
%% Male vs. Female
% Preallocate
[ffA,mmA] = deal(zeros(1,5));
[fmA,mfA] = deal(zeros(5,5));
for ii = 1:5
   load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
       'maleFemale\maleFemale',num2str(ii),'.mat']) 
   ffA(ii) = femaleMod.auc;
   fmA(ii,:) = femaleMod.mAuc;
   mmA(ii) = maleMod.auc;
   mfA(ii,:) = maleMod.fAuc;
end
% Plot
figure
subplot(1,2,1)
hold on
scatter(1:5,mmA,200,'.k')
scatter(reshape(repmat(1:5,5,1),1,25),reshape(fmA,1,25),200,'.r')
set(gca,'Xtick',1:5)
xlabel('Test Animal'); ylabel('AUC')
legend({'Male Model','Female Model'},'location','nw')
title('Male')

subplot(1,2,2)
hold on
scatter(1:5,ffA,200,'.k')
scatter(reshape(repmat(1:5,5,1),1,25),reshape(mfA,1,25),200,'.r')
set(gca,'Xtick',1:5)
xlabel('Test Animal'); ylabel('AUC')
legend({'Female Model','Male Model'},'location','nw')
title('Female')
%% Build dep24 vs. chow models
load(['E:\paper2\analyzed\finalNew\bingeNotPreData.mat'])
% Get number of samples from dep24 and chow
bingeSamp = [cellfun(@(x) size(x,1),data{1,2}([1:3,5:8,10:11],1)),...
    cellfun(@(x) size(x,1),data{1,4}(:,1))];
dep = data{1,2}([1:3,5:8,10:11],1);
chow = data{1,4}(:,1);
% Find min for each animal
[minSamp] = min(bingeSamp,[],2);
% Set up formula
formula = [];
for ii = 1:60
    formula = [formula,'trainX',num2str(ii),'+'];
end
formula = ['Y~',formula,'(1|ID)'];
for jj = 1:100
    % Build models between chow and dep24
    allDep = []; allChow = []; ID = [];
    for ii = 1:size(minSamp,1)
        thisDep = dep{ii}(randperm(size(dep{ii},1),minSamp(ii)),:);
        allDep = [allDep;thisDep];
        thisChow = chow{ii}(randperm(size(chow{ii},1),minSamp(ii)),:);
        allChow = [allChow;thisChow];
        ID = [ID;ones(minSamp(ii),1).*ii];
    end
    % Determine number of samples from chow and dep to use in training
    halfTrain = round(size(allChow,1)*.8);
    halfTest = size(allChow,1)-halfTrain;
    depTrainInds = randperm(size(allDep,1),halfTrain);
    depTestInds = ~ismember(1:size(allDep,1),depTrainInds);
    chowTrainInds = randperm(size(allChow,1),halfTrain);
    chowTestInds = ~ismember(1:size(allChow,1),chowTrainInds);
    trainX = [allDep(depTrainInds,:);allChow(chowTrainInds,:)];
    trainY = [zeros(halfTrain,1);ones(halfTrain,1)];
    testX = [allDep(depTestInds,:);allChow(chowTestInds,:)];
    testY = [zeros(halfTest,1);ones(halfTest,1)];
    
    train = array2table(trainX);
    train.ID = [ID(depTrainInds,:);ID(chowTrainInds,:)];
    train.Y = trainY;
    
    test = array2table(testX);
    test.ID = [ID(depTestInds,:);ID(chowTestInds,:)];
    test.Y = testY;
    test.Properties.VariableNames = train.Properties.VariableNames;
    
%     mdl = fitglm(trainX,trainY,'distribution','binomial','binomialsize',...
%         halfTrain*2);
    mdl = fitglme(train,formula,'distribution','binomial',...
        'binomialsize',halfTrain*2);
    prob = predict(mdl,test);
    [dcX(jj,:),dcY(jj,:),~,dcA(jj)] = perfcurve(testY,prob,1);
    prob = predict(mdl,test(randperm(size(test,1)),:));
    [dcrX(jj,:),dcrY(jj,:),~,dcrA(jj)] = perfcurve(testY,prob,1);
end
figure
plot(mean(dcX,1),mean(dcY,1))
hold on
plot(mean(dcrX,1),mean(dcrY,1),'--k')
box off
legend({['Real: ',num2str(round(mean(dcA),2)),'\pm',...
    num2str(round(conf(dcA,0.95),3))],['Permuted: ',...
    num2str(round(mean(dcrA),2)),'\pm',num2str(round(conf(dcrA,0.95),3))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
title('Dep24 vs. Chow')
%% Build models comparing dep24 to chow
load(['E:\paper2\analyzed\finalNew\bingeNotPreData.mat'])
bingeSamp = [cellfun(@(x) size(x,1),data{1,2}([1:3,5:8,10:11],1)),...
    cellfun(@(x) size(x,1),data{1,4}(:,1))];
dep = data{1,2}([1:3,5:8,10:11],1);
depNot = data{1,2}([1:3,5:8,10:11],2);
chow = data{1,4}(:,1);
chowNot = data{1,4}(:,2);
% Find min for each animal
[minSamp] = min(bingeSamp,[],2);
for ii = 1:100
    allDep = []; allDepNot = []; allChow = []; allChowNot = [];
    for jj = 1:size(minSamp,1)
        thisDep = dep{ii}(randperm(size(dep{ii},1),minSamp(ii)),:);
        allDep = [allDep;thisDep];
        thisDepNot = depNot{ii}(randperm(size(depNot{ii},1),minSamp(ii)),:);
        allDepNot = [allDepNot;thisDepNot];
        
        thisChow = chow{ii}(randperm(size(chow{ii},1),minSamp(ii)),:);
        allChow = [allChow;thisChow];
        thisChowNot = chowNot{ii}(randperm(size(chowNot{ii},1),minSamp(ii)),:);
        allChowNot = [allChowNot;thisChowNot];
    end
    % Determine train and test sizes
    halfTrain = round(size(allChow,1)*.8);
    halfTest = size(allChow,1)-halfTrain;
    % Generate indices
    
    % Dep -> Dep
   
    % Dep -> Chow
   
    % Chow -> Chow
   
    % Chow -> Dep
end
%% Build models: binge, drink, and across
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'baseline500Each6000All50-50_100.mat'])
binge = all;
% bingeRnd = rnd;
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000_100.mat'],'all','rnd')
drink = all;
% drinkRnd = rnd;
% Grab only NAc indices for comparison purposes
bingeInds = [1:12,25:30];
drinkInds = [13:24,55:60];
% Grab prefrontal indices
% drinkInds = [1:12,25:30];
% Set step for TVals
tStep = 1/2400;
% Preallocate
[bbX,bbY,bbXperm,bbYperm,bdX,bdY,bdXperm,bdYperm,ddX,ddY,ddXperm,...
    ddYperm,dbX,dbY,dbXperm,dbYperm] = deal(cell(1,100));
[bbA,bbAperm,bdA,bdAperm,ddA,ddAperm,dbA,dbAperm] = deal(zeros(1,100));
for ii = 1:100
    % Binge -> Binge
    trainX = binge.trainX{ii}(:,bingeInds);
    trainY = binge.trainY{ii};
    testX = binge.testX{ii}(:,bingeInds);
    testY = binge.testY{ii};
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [bbX{ii},bbY{ii},~,bbA(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/1200:1,'UseNearest',0);
    [bbXperm{ii},bbYperm{ii},~,bbAperm(ii)] = perfcurve(...
        testY(randperm(size(testY,1))),prob,1,'TVals',0:tStep:1,...
        'UseNearest',0);
    % Binge -> Drink
    testX = drink.testX{ii}(:,drinkInds);
    testY = drink.testY{ii};
    prob = predict(mdl,testX);
    [bdX{ii},bdY{ii},~,bdA(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:tStep:1,'UseNearest',0);
    [bdXperm{ii},bdYperm{ii},~,bdAperm(ii)] = perfcurve(...
    testY(randperm(size(testY,1))),prob,1,'TVals',0:tStep:1,...
    'UseNearest',0);
    % Drink -> Drink
    trainX = drink.trainX{ii}(:,drinkInds);
    trainY = drink.trainY{ii};
    testX = drink.testX{ii}(:,drinkInds);
    testY = drink.testY{ii};
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [ddX{ii},ddY{ii},~,ddA(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/1200:1,'UseNearest',0);
    [ddXperm{ii},ddYperm{ii},~,ddAperm(ii)] = perfcurve(...
        testY(randperm(size(testY,1))),prob,1,'TVals',0:tStep:1,...
    'UseNearest',0);
    % Drink -> Binge
    testX = binge.testX{ii}(:,bingeInds);
    testY = binge.testY{ii};
    prob = predict(mdl,testX);
    [dbX{ii},dbY{ii},~,dbA(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/1200:1,'UseNearest',0);
    [dbXperm{ii},dbYperm{ii},~,dbAperm(ii)] = perfcurve(...
        testY(randperm(size(testY,1))),prob,1,'TVals',0:tStep:1,...
        'UseNearest',0);
end
% Confidence intervals - effect sizes are meaningless
% bbES = distES(bbA,bbAperm);
bbConf = conf(bbA,0.95);
bbPermConf = conf(bbAperm,0.95);
% bdES = distES(bdA,bdAperm);
bdConf = conf(bdA,0.95);
bdPermConf = conf(bdAperm,0.95);
% ddES = distES(ddA,ddAperm);
ddConf = conf(ddA,0.95);
ddPermConf = conf(ddAperm,0.95);
% dbES = distES(dbA,dbAperm);
dbConf = conf(dbA,0.95);
dbPermConf = conf(dbAperm,0.95);
%% Build and test generalized models
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'baseline500Each6000All50-50_100.mat'])
binge = all;
bingeRnd = rnd; %#ok
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000_100.mat'],'all','rnd')
drink = all;
drinkRnd = rnd; %#ok
bingeInds = [1:12,25:30];
drinkInds = [13:24,55:60];
% Set step for TVals
tStep = 1/2400;
% Preallocate
[ggX,ggY,ggXperm,ggYperm,gbX,gbY,gbXperm,gbYperm,gdX,gdY,gdXperm,...
    gdYperm] = deal(cell(1,100));
[ggA,ggAperm,gbA,gbAperm,gdA,gdAperm] = deal(zeros(1,100));
betas = zeros(18,100);
for ii = 1:100
    % Train and test gen model
    trainX = [binge.trainX{ii}(:,bingeInds);drink.trainX{ii}(:,drinkInds)];
    trainY = [binge.trainY{ii};drink.trainY{ii}];
    testX = [binge.testX{ii}(:,bingeInds);drink.testX{ii}(:,drinkInds)];
    testY = [binge.testY{ii};drink.testY{ii}];
    % Normal
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    betas(:,ii) = table2array(mdl.Coefficients(2:end,1));
    prob = predict(mdl,testX);
    [ggX{ii},ggY{ii},~,ggA(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:tStep:1,'UseNearest',0);
    % Permuted
    [ggXperm{ii},ggYperm{ii},~,ggAperm(ii)] = perfcurve(...
        testY(randperm(size(testY,1))),prob,1,'TVals',0:tStep:1,...
        'UseNearest',0);
    % Apply gen model to binge
    testX = binge.testX{ii}(:,bingeInds);
    testY = binge.testY{ii};
    prob = predict(mdl,testX);
    % Normal
    [gbX{ii},gbY{ii},~,gbA(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:tStep:1,'UseNearest',0);
    % Permuted
    [gbXperm{ii},gbYperm{ii},~,gbAperm(ii)] = perfcurve(...
        testY(randperm(size(testY,1))),prob,1,'TVals',0:tStep:1,...
        'UseNearest',0);
    % Apply gen model to drink
    testX = drink.testX{ii}(:,drinkInds);
    testY = drink.testY{ii};
    prob = predict(mdl,testX);
    % Normal
    [gdX{ii},gdY{ii},~,gdA(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:tStep:1,'UseNearest',0);
    % Permuted
    [gdXperm{ii},gdYperm{ii},~,gdAperm(ii)] = perfcurve(...
        testY(randperm(size(testY,1))),prob,1,'TVals',0:tStep:1,...
        'UseNearest',0);
end
% Confidence intervals - effect sizes are meaningless
% ggES = distES(ggA,ggAperm);
ggConf = conf(ggA,0.95);
ggPermConf = conf(ggAperm,0.95);
% gbES = distES(gbA,gbAperm);
gbConf = conf(gbA,0.95);
gbPermConf = conf(gbAperm,0.95);
% gdES = distES(gdA,gdAperm);
gdConf = conf(gdA,0.95);
gdPermConf = conf(gdAperm,0.95);
%% Plot models
figure
% Drink>Drink
subplot(3,2,1); hold on
plot(mean(cat(2,ddX{:}),2),mean(cat(2,ddY{:}),2),'-k')
plot(mean(cat(2,ddXperm{:}),2),mean(cat(2,ddYperm{:}),2),'--k')
legend({['Actual: ',num2str(round(mean(ddA),2)),'\pm',...
    num2str(round(ddConf,2))],['Perm: ',num2str(round(mean(ddAperm),2)),...
    '\pm',num2str(round(ddPermConf,2))]},'location','se')
box off
title('Drink>Drink')
% Drink>Binge
subplot(3,2,2); hold on
plot(mean(cat(2,dbX{:}),2),mean(cat(2,dbY{:}),2),'-k')
plot(mean(cat(2,dbXperm{:}),2),mean(cat(2,dbYperm{:}),2),'--k')
legend({['Actual: ',num2str(round(mean(dbA),2)),'\pm',...
    num2str(round(dbConf,2))],['Perm: ',num2str(round(mean(dbAperm),2)),...
    '\pm',num2str(round(dbPermConf,2))]},'location','se')
box off
title('Drink>Binge')
% Binge>Binge
subplot(3,2,3); hold on
plot(mean(cat(2,bbX{:}),2),mean(cat(2,bbY{:}),2),'-k')
plot(mean(cat(2,bbXperm{:}),2),mean(cat(2,bbYperm{:}),2),'--k')
legend({['Actual: ',num2str(round(mean(bbA),2)),'\pm',...
    num2str(round(bbConf,2))],['Perm: ',num2str(round(mean(bbAperm),2)),...
    '\pm',num2str(round(bbPermConf,2))]},'location','se')
box off
title('Binge>Binge')
% Binge>Drink
subplot(3,2,4); hold on
plot(mean(cat(2,bdX{:}),2),mean(cat(2,bdY{:}),2),'-k')
plot(mean(cat(2,bdXperm{:}),2),mean(cat(2,bdYperm{:}),2),'--k')
legend({['Actual: ',num2str(round(mean(bdA),2)),'\pm',...
    num2str(round(bdConf,2))],['Perm: ',num2str(round(mean(bdAperm),2)),...
    '\pm',num2str(round(bdPermConf,2))]},'location','se')
box off
title('Binge>Drink')
% Gen>Drink
subplot(3,2,5); hold on
plot(mean(cat(2,gdX{:}),2),mean(cat(2,gdY{:}),2),'-k')
plot(mean(cat(2,gdXperm{:}),2),mean(cat(2,gdYperm{:}),2),'--k')
legend({['Actual: ',num2str(round(mean(gdA),2)),'\pm',...
    num2str(round(gdConf,2))],['Perm: ',num2str(round(mean(gdAperm),2)),...
    '\pm',num2str(round(gdPermConf,2))]},'location','se')
box off
title('Gen>Drink')
% Gen>Binge
subplot(3,2,6); hold on
plot(mean(cat(2,gbX{:}),2),mean(cat(2,gbY{:}),2),'-k')
plot(mean(cat(2,gbXperm{:}),2),mean(cat(2,gbYperm{:}),2),'--k')
legend({['Actual: ',num2str(round(mean(gbA),2)),'\pm',...
    num2str(round(gbConf,2))],['Perm: ',num2str(round(mean(gbAperm),2)),...
    '\pm',num2str(round(gbPermConf,2))]},'location','se')
box off
title('Gen>Binge')
%% Gen>Gen
figure; hold on
plot(mean(cat(2,ggX{:}),2),mean(cat(2,ggY{:}),2),'-k')
plot(mean(cat(2,ggXperm{:}),2),mean(cat(2,ggYperm{:}),2),'--k')
legend({['Actual: ',num2str(round(mean(ggA),2)),'\pm',...
    num2str(round(ggConf,2))],['Perm: ',num2str(round(mean(ggAperm),2)),...
    '\pm',num2str(round(ggPermConf,2))]},'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('Gen>Gen')
%% Monads
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'baseline500Each6000All50-50_100.mat'])
binge = all;
bingeRnd = rnd;
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000_100.mat'],'all','rnd')
drink = all;
drinkRnd = rnd;
bingeInds = [1:12,25:30];
drinkInds = [13:24,55:60];
% Preallocate
[bingeA,drinkA] = deal(cell(1,3));
for k = 1%:3
    cmbs = nchoosek(1:18,k);
    for ii = 1:100
        disp([num2str(k),': ',num2str(ii)])
        for jj = 1:size(cmbs,1)
            % Train and test gen model
            genTrainX = [binge.trainX{ii}(:,bingeInds);...
                drink.trainX{ii}(:,drinkInds)];
            genTrainY = [binge.trainY{ii};drink.trainY{ii}];
            genTestX = [binge.testX{ii}(:,bingeInds);...
                drink.testX{ii}(:,drinkInds)];
            genTestY = [binge.testY{ii};drink.testY{ii}];
            mdl = fitglm(genTrainX(:,cmbs(jj,:)),genTrainY,...
                'distribution','binomial');
            prob = predict(mdl,genTestX(:,cmbs(jj,:)));
            [~,~,~,a{k}(ii,jj)] = perfcurve(genTestY,prob,1);
            % Binge test
            testXbinge = binge.testX{ii}(:,bingeInds);
            testYbinge = [binge.testY{ii}];
            bingeProb = predict(mdl,testXbinge(:,cmbs(jj,:)));
            [~,~,~,bingeA{k}(ii,jj)] = perfcurve(testYbinge,bingeProb,1);
            % Drink test
            testXdrink = drink.testX{ii}(:,drinkInds);
            testYdrink = [drink.testY{ii}];
            drinkProb = predict(mdl,testXdrink(:,cmbs(jj,:)));
            [~,~,~,drinkA{k}(ii,jj)] = perfcurve(testYdrink,drinkProb,1);            
        end
    end
end
drinkAM = mean(drinkA{1},1)';
bingeAM = mean(bingeA{1},1)';
% genAM = mean(genA{1},1)';
mOdds = mean(exp(betas),2);
mBetas = mean(betas,2);
%%
figure
hold on
plot([0.5 0.5],[0 1],'--','col',[0.5 0.5 0.5])
plot([0 1],[0.5 0.5],'--','col',[0.5 0.5 0.5])
scatter(bingeAM,drinkAM,abs(mBetas)+1,'ok')
xlabel('Binge AUC'); ylabel('Drink AUC')
set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1)
title('Betas')
%% Sort vars
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
nameVect = nameVect(bingeInds);
% Preallocate
[aSort,sortInd,vars] = deal(cell(1,3));
for k = 1%:3
   [aSort{k},sortInd{k}] = sort(mean(a{k},1)','descend'); %#ok
   cmbs = nchoosek(1:18,k);
   for ii = 1:k
       vars{k}(:,ii) = nameVect(cmbs(sortInd{k},ii));
   end
end
%% Compare models built with different data subsets (all,pfc,nac)
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drinkNotSubsets.mat'])
figure
hold on
plot(mean(cat(2,ddXall{:}),2),mean(cat(2,ddYall{:}),2))
plot(mean(cat(2,ddXnac{:}),2),mean(cat(2,ddYnac{:}),2))
plot(mean(cat(2,ddXpfc{:}),2),mean(cat(2,ddYpfc{:}),2))
legend({['All: ',num2str(round(mean(ddAall),2)),'\pm',...
    num2str(round(conf(ddAall,0.95),2))],...
    ['NAc: ',num2str(round(mean(ddAnac),2)),'\pm',...
    num2str(round(conf(ddAnac,0.95),2))],...
    ['PFC: ',num2str(round(mean(ddApfc),2)),'\pm',...
    num2str(round(conf(ddApfc,0.95),2))]},'location','se')
title('DrinkNot: subsets')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
%% Pre-Drinking vs. Not Drinking
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper3\'...
    'preDrinkCombined2\'],'.mat');
preData = [];
notData = [];
for ii = 1:size(files,2)
   load(files{ii})
   if ~isempty(trls{1,1})
      [b,c,t] = size(psdTrls{1,1}.relPow);
      thisPow = reshape(psdTrls{1,1}.relPow,b*c,t)';
      [cmb,b,t] = size(coh{1,1}.rel);
      thisCoh = reshape(permute(coh{1,1}.rel,[2,1,3]),cmb*b,t)';
      preData = [preData;thisPow,thisCoh]; %#ok
      
      [b,c,t] = size(psdTrls{1,end}.relPow);
      thisPow = reshape(psdTrls{1,end}.relPow,b*c,t)';
      [cmb,b,t] = size(coh{1,end}.rel);
      thisCoh = reshape(permute(coh{1,end}.rel,[2,1,3]),cmb*b,t)';
      notData = [notData;thisPow,thisCoh]; %#ok
   end
end
% Collate and impute
catData{1,1} = preData;
catData{1,2} = notData;
% [all,each,rnd,~] = evenDataSplit(catData,17222,4306,'ADA',20);
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
%     'preDrinkData.mat'],'all','each','rnd')
[all,each,rnd,~] = evenDataSplit(catData,17222,4306,'ADA',100);
save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'preDrinkData100.mat'],'all','each','rnd')
%% Build full logistics
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'preDrinkData100.mat'])
% Preallocate
[rndX,rndY] = deal(cell(1,20));
rndA = zeros(1,20);
for n = 1:20
    disp([num2str(n),' of 20'])
    % Set up training data
    trainX = all.trainX{n};
    trainY = all.trainY{n};
    % Set up testing data
    testX = all.testX{n};
    testY = all.testY{n};
    % Build and test model on concat data
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [x{n},y{n},~,a(n)] = perfcurve(testY,prob,1);
    % Set up random testing data
    testY = rnd.allTestY{n};
    prob = predict(mdl,testX);
    [rndX{n},rndY{n},~,rndA(n)] = perfcurve(testY,prob,1);
end
% Test drinkNot models on preDrink data
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drinkNotModel100.mat'],'concatData')
% load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
%     'preDrinkData.mat'])
% Preallocate
[drinkX,drinkY] = deal(cell(1,20));
drinkA = zeros(1,20);
for ii = 1:20
    disp([num2str(ii),' of 20'])
    testX = all.testX{ii};
    testY = all.testY{ii};
    prob = predict(concatData.model{ii},testX);
    [drinkX{ii},drinkY{ii},~,drinkA(ii)] = perfcurve(testY,prob,1);
end
% Load preBinge files
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\'...
    'preBinge\'])
% Preallocate
[preA,preARand] = deal(zeros(20,8));
% Hardcoded first dimension from known size of roc curves in concatData
[preX,preY] = deal(zeros(21207,20));
[preRandX,preRandY] = deal(zeros(21222,20));
beta = zeros(20,58);
for ii = 1:20
   load([num2str(ii),'.mat'])
   preA(ii,:) = concatData{1,8}.auc;
   preX(:,ii) = concatData{1,1}.acc{1,1}.x;
   preY(:,ii) = concatData{1,1}.acc{1,1}.y;
end
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\'...
    'preBingeRand\'])
for ii = 1:20
   load([num2str(ii),'.mat'])
   preARand(ii,:) = concatData{1,8}.auc;
   preRandX(:,ii) = concatData{1,1}.acc{1,1}.x;
   preRandY(:,ii) = concatData{1,1}.acc{1,1}.y;
end
% Plot
figure
hold on
plot(mean(cat(2,x{:}),2),mean(cat(2,y{:}),2),'-k')
plot(mean(preX,2),mean(preY,2),'-.k')
plot(mean(cat(2,rndX{:}),2),mean(cat(2,rndY{:}),2),'--k')
plot(mean(preRandX,2),mean(preRandY,2),':k')
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend({['PreDrink: ',num2str(round(mean(a),2)),'\pm',...
    num2str(round(conf(a,0.95),2))],['PreBinge: ',...
    num2str(round(mean(preA(:,1)),2)),'\pm',...
    num2str(round(conf(preA(:,1)',0.95),2))],['Permuted Drink: ',...
    num2str(round(mean(rndA),2)),'\pm',...
    num2str(round(conf(rndA,0.95),2))],['Permuted Binge: ',...
    num2str(round(mean(preARand(:,1)),2)),'\pm',...
    num2str(round(conf(preARand(:,1)',0.95),2))]},'Location','se')
title('Pre-Drink vs. Other')
%% Set up config for analyzing preBinge data up to 240 seconds before the
% start of a binge moving in 1 second intervals
cfg.sdir = 'C:\Users\Pythia\Documents\GreenLab\data\paper2\toProcess\';
cfg.file = file;
cfg.nFilt = [57 63];
cfg.dsf = 5;
cfg.thresh = 1; 
cfg.onset = 1;
cfg.offset = 1;
cfg.foi = [1 1 100];
cfg.bands = {'delta',[1,4];
         'theta',[5,10];
         'alpha',[11,14];
         'beta',[15,30];
         'lgamma',[45,65];
         'hgamma',[70,90]};
cfg.overlap = 0.5;
cfg.ohMethod = 'mtm';
for ei = 1:240
    cfg.eoi(ei,:) = {'binge (s',[-4-ei 1-ei]};
end
cfg.vis = 'n';
cfg.saveParent = ['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    'preBinge_240sec\'];
%% Combine prebinge data with non-overlapping non-binge data
% Get fileNames of preBinge data
fNames = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    'preBinge_240sec\'],'.mat');
for ii = 1:size(fNames,2)
    names(ii,:) = strsplit(fNames{ii},'_'); %#ok<SAGROW>
end
% Cycle through files
% Preallocate
times = cell(1,size(names,1));
for ii = 1:size(names,1)
    disp(num2str(ii))
    load(fNames{ii})
    % Separate notBinge and preBinge
    notCoh = coh{end};
    notPow = psdTrls{end};
    notTrls = trls{end};
    
    preCoh = coh(1:end-1);
    prePow = psdTrls(1:end-1);
    preTrls = trls(1:end-1);
    % Grab times of binges
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\mat\',...
        names{ii,1},'_',names{ii,2},'.mat'],'eventTs')
    ind = logicFind(1,cell2mat(cellfun(@(x) strncmpi('binge (s',x,8),...
        eventTs.label,'UniformOutput',0)),'==');
    times{ii} = [eventTs.t{ind}*hist.adfreq,eventTs.t{ind+1}*hist.adfreq];
    % Clear generic variables so they can be used later
    clear coh psdTrls trls
    %% Check that pre trials do no overlap with binges
    % Preallocate
    nTrl = size(prePow,2); 
    [trls,psdTrls,coh] = deal(cell(1,nTrl));
    for k = 1:nTrl
        % Preallocate
        nSamp = size(preTrls{1,k}.sampleinfo,1);
        thisOverlap = zeros(nSamp,size(times{ii},1));
        overlap = zeros(1,nSamp);
        for jj = 1:nSamp
            for m = 1:size(times{ii},1)
                thisTrl = ismember(preTrls{1,k}.sampleinfo(jj,1):...
                    preTrls{1,k}.sampleinfo(jj,2),...
                    round(times{ii}(m,1)):round(times{ii}(m,2)));
                thisOverlap(jj,m) = any(thisTrl(:) == 1);
            end
        overlap(jj) = any(thisOverlap(jj,:)==1);
        end
        trls{1,k}.trial = preTrls{1,k}.trial(:,:,~overlap);
        trls{1,k}.time = preTrls{1,k}.time(~overlap);
        trls{1,k}.sampleinfo = preTrls{1,k}.sampleinfo(~overlap,:);
        
        psdTrls{1,k}.Pow = prePow{1,k}.Pow(:,:,~overlap);
        psdTrls{1,k}.f = prePow{1,k}.f;
        psdTrls{1,k}.hammSize = prePow{1,k}.hammSize;
        psdTrls{1,k}.bandPow = prePow{1,k}.bandPow(:,:,~overlap);
        psdTrls{1,k}.totPow = prePow{1,k}.totPow(:,:,~overlap);
        psdTrls{1,k}.relPow = prePow{1,k}.relPow(:,:,~overlap);
        psdTrls{1,k}.Pow = prePow{1,k}.Pow(:,:,~overlap);
        
        coh{1,k}.Cxy = preCoh{1,k}.Cxy(:,:,~overlap);
        coh{1,k}.rel = preCoh{1,k}.normBandCoh(:,:,~overlap);
        coh{1,k}.band = preCoh{1,k}.mBandCoh(:,:,~overlap);
        coh{1,k}.mRaw = preCoh{1,k}.mtCxy(:,:,~overlap);
        coh{1,k}.f = preCoh{1,k}.f;
    end
    % Find notTrls that do not overlap with preTrls
    samps = [];
    for jj = 1:size(trls,2)
        samps = [samps;trls{1,jj}.sampleinfo]; %#ok<AGROW>
    end
    overlap = zeros(1,nSamp);
    for jj = 1:size(notTrls.sampleinfo,1)
        thisTrl = ismember(samps,notTrls.sampleinfo(jj,1):...
            notTrls.sampleinfo(jj,2));
        overlap(jj) = any(thisTrl(:) == 1);
    end
    % Extract non-overlapping trials
    trls{1,k+1}.label = notTrls.label;
    trls{1,k+1}.fsample = notTrls.fsample;
    trls{1,k+1}.trial = notTrls.trial(:,:,~overlap);
    trls{1,k+1}.time = notTrls.time(~overlap);
    trls{1,k+1}.sampleinfo = notTrls.sampleinfo(~overlap,:);
    
    psdTrls{1,k+1}.Pow = notPow.Pow(:,:,~overlap);
    psdTrls{1,k+1}.f = notPow.f;
    psdTrls{1,k+1}.hammSize = notPow.hammSize;
    psdTrls{1,k+1}.bandPow = notPow.bandPow(:,:,~overlap);
    psdTrls{1,k+1}.totPow = notPow.totPow(:,:,~overlap);
    psdTrls{1,k+1}.relPow = notPow.relPow(:,:,~overlap);
    psdTrls{1,k+1}.Pow = notPow.Pow(:,:,~overlap);
    
    coh{1,k+1}.Cxy = notCoh.Cxy(:,:,~overlap);
    coh{1,k+1}.rel = notCoh.normBandCoh(:,:,~overlap);
    coh{1,k+1}.band = notCoh.mBandCoh(:,:,~overlap);
    coh{1,k+1}.mRaw = notCoh.mtCxy(:,:,~overlap);
    coh{1,k+1}.f = notCoh.f;
    overlap = [];
    
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\'...
        'preBingeCombined\',names{ii,1},'_',names{ii,2},'.mat'],'trls',...
        'psdTrls','coh')
end
%% Prep data for preBinge model
[data,samp,preFiles] = collateData(['C:\Users\Pythia\Documents\GreenLab'...
    '\data\paper3\preBingeCombined\'],{'base'},{'pow','coh'},'trl','rel');
% Preallocate
notBinge = cell(size(data{1},1),1);
% Grab bingeNot data that does not overlap with preBinge
for ii = 1:size(preFiles{1},2)
    disp(['Grabbin'' data from file ',num2str(ii),' of ',...
        num2str(size(preFiles{1},2))])
    % Grab all samples from this row's preBinge data
    theseSamp = cat(1,samp{1}{ii,1:end-1});
    sampCat = [];
    for jj = 1:size(theseSamp,1)
       sampCat = [sampCat,theseSamp(jj,1):theseSamp(jj,2)]; %#ok
    end
    % Take only unique values
    sampCat = unique(sampCat);
    % Check if notBinge data overlaps with sampCat
    overlap = sum(ismember(samp{1}{ii,end},sampCat),2);
    % If overlap remove those trials
    notBinge{ii,1} = data{1,1}{ii,end}(~overlap,:);
end
% Set up 5 sec preBinge model with notBinge data used for ADASYN
preData = data{1}(:,1:end-1);
pre5Cat = cat(1,preData{:,1});
% Grab notBinge data
notBinge = cat(1,data{1}{:,end});
% Determine ratio of preBinge to notBinge
nPre = sum(cellfun(@(x) size(x,1),samp{1}(:,1)));
nNot = size(notBinge,1);
r = nPre/(nPre+nNot);
% Use 21 samples per animal (252 total); thus, need 126 preBinge and 126
% notBinge. Use 1000 samples for testing due to extreme rarity of preBinge
% compared to notBinge.
nPreTest = round(r*1000);
% Preallocate
[thisTrainX,newPre,thisNotBinge] = deal(cell(1,20));
for ii = 1:100
    % Randomally generate indices corresponding to preBinge values to be
    % set aside for testing
    rPBI = randperm(size(pre5Cat,1),nPreTest);
    pre5.testX{ii} = pre5Cat(rPBI,:);
    % Take rest of values for training
    thisTrainX{ii} = pre5Cat(~ismember(1:size(pre5Cat,1),rPBI),:); 
    % Use 'notBinge' as majority case in ADASYN to impute preBinge up to
    % 126; uses 200 to ensure enough samples are generated
    [newPre{ii},~] = ADASYN([thisTrainX{ii};...
        notBinge(randperm(size(notBinge,1),200),:)],...
        [ones(size(thisTrainX{ii},1),1);zeros(200,1)]);
    % Concatenate real and imputed data together with notBinge
    thisPreCat = cat(1,thisTrainX{ii},newPre{ii});
    % Generate 1119 indices for notBinge; use first 126 for training and
    % last 993 for testing
    rNBI = randperm(size(notBinge,1),1119);
    % Grab 126 preBinge and 126 notBinge
    pre5.trainX{ii} = cat(1,thisPreCat(randperm(size(...
        thisPreCat,1),126),:),notBinge(rNBI(1:126),:));
    % Add notBinge to test set
    pre5.testX{ii} = [pre5.testX{ii};notBinge(rNBI(127:1119),:)];
    % Create Ys for training and testing
    pre5.trainY{ii} = [ones(126,1);zeros(126,1)];
    pre5.testY{ii} = [ones(nPreTest,1);zeros(993,1)];
    % Create test sets for pre240
    thisNotBinge{ii} = notBinge(~ismember(1:nNot,rNBI),:);
    % Use 50 preBinge trials (minimum number across all 240), which means
    % 7093 notBinge trials to approximate the ratio r
    % --CREATES HUGE STRUCTURE-- manually create each iterate on Discovery
    % instead
%     for jj = 2:240
%         thisPreCat = cat(1,preData{:,jj});
%         pre240.testX{ii,jj} = cat(1,...
%             thisPreCat(randperm(size(thisPreCat,1),50),:),...
%             thisNotBinge(randperm(size(thisNotBinge,1),7093),:));
%         pre240.testY{ii,jj} = cat(1,ones(50,1),zeros(7093,1));
%     end
end
% Save
save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'preBingeModelData.mat'],'pre5','thisNotBinge','preData')
%% Open binge240 files
[lassoA,logA,lassoRandA,logRandA] = deal(zeros(20,240));
for ii = 1:100
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
        'preBinge240\',num2str(ii),'.mat'])
    lassoA(ii,:) = cellfun(@(x) x.auc,preBingeLasso);
    logA(ii,:) = cellfun(@(x) x.auc,preBingeLog);
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
        'preBinge240\',num2str(ii),'_Rand.mat'])
    lassoRandA(ii,:) = cellfun(@(x) x.auc,preBingeLasso);
    logRandA(ii,:) = cellfun(@(x) x.auc,preBingeLog);
end
%%
% Get mean and std
lassoAM = mean(lassoA,1);
lassoAS = std(lassoA,[],1);
lassoAC = conf(lassoA',0.95);
logAM = mean(logA,1);
logAS = std(logA,[],1);

lassoRandAM = mean(lassoRandA,1);
lassoRandAS = std(lassoRandA,[],1);
lassoRandAC = conf(lassoRandA',0.95);
logRandAM = mean(logRandA,1);
logRandAS = std(logRandA,[],1);
% Run t-tests w/ correction
[~,~,pAdjLasso] = bulkT(lassoA,lassoRandA,1,'bc');
lassoInd = logicFind(0.05,pAdjLasso,'>=');
[~,~,pAdjLog] = bulkT(logA,logRandA,1,'bc');
logInd = logicFind(0.05,pAdjLog,'>=');
% Plot - scatterErr
% scatterErr(1:240,lassoAM,lassoAS,1)
% hold on
% scatterErr(1:240,lassoRandAM,lassoRandAS,0)
% 
% scatterErr(1:240,logAM,logAS,1)
% hold on
% scatterErr(1:240,logRandAM,logRandAS,0)
% Plot - shadedError
figure; hold on
shadedErrorBar(1:240,lassoAM,lassoAS,'b',1)
shadedErrorBar(1:240,lassoRandAM,lassoRandAS,'k',1)
plot(lassoInd,lassoRandAM(lassoInd),'.r')

figure; hold on
shadedErrorBar(1:240,logAM,logAS,'b',1)
shadedErrorBar(1:240,logRandAM,logRandAS,'k',1)
plot(logInd,logRandAM(logInd),'.r')
%% Prep bingeNot data for low samples -> ADASYN to compare w/ low sample 
% drinking
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\bingeNotData.mat'],'trlDat')
[all.trainX,all.trainY,all.testX,all.testY,rnd.trainX,rnd.trainY,...
    rnd.testX,rnd.testY] = deal(cell(1,20));
for ii = 1:20
    % Randomally grab 380 binge trials
    newBinge = trlDat{1,1}(randperm(size(trlDat{1},1),380),:);
    newData = {newBinge,trlDat{2,1}};
    [thisAll,~,thisRnd] = evenDataSplit(newData,6000,1200,'ADA',1);
    all.trainX{ii} = thisAll.trainX{1};
    all.trainY{ii} = thisAll.trainY{1};
    all.testX{ii} = thisAll.testX{1};
    all.testY{ii} = thisAll.testY{1};
    
    rnd.trainY{ii} = thisRnd.allTrainY{1};
    rnd.testY{ii} = thisRnd.allTestY{1};
end
% save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
%     'bingeNotLowSampleADA.mat'],'all','rnd');
%% Build and test low sample bingeNot models
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'bingeNotLowSampleADA.mat'])
% Preallocate
[x,y,t,xRand,yRand] = deal(cell(1,20));
[a,aRand] = deal(zeros(1,20));
for n = 1:20
    disp([num2str(n),' of 20'])
    % Set up training data
    trainX = all.trainX{n};
    trainY = all.trainY{n};
    % Set up testing data
    testX = all.testX{n};
    testY = all.testY{n};
    % Build and test model on concat data
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [x{n},y{n},~,a(n)] = perfcurve(testY,prob,1,'TVals',0:1/4623:1,...
        'UseNearest',0);
    % Store outcomes
    concatData.model{n}  = mdl;
%     concatData{n}.rocX = x;
%     concatData{n}.rocY = y;
%     concatData{n}.auc = a;
    % Set up testing data
    trainY = rnd.trainY{n};
    testY = rnd.testY{n};
    % Build and test model on concat data
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [xRand{n},yRand{n},~,aRand(n)] = perfcurve(testY,prob,1,'Tvals',...
        0:1/4623:1,'UseNearest',0);
end
concatData.rocX = x;
concatData.rocY = y;
concatData.auc = a;
% Get averages
mX = mean(cat(2,x{:}),2);
mY = mean(cat(2,y{:}),2);
mRandX = mean(cat(2,xRand{:}),2);
mRandY = mean(cat(2,yRand{:}),2);
% Plot
figure
hold on
for ii = 1:20
    p = plot(x{ii},y{ii});
    col = get(p,'Color');
    hsv = rgb2hsv(col);
    newCol = hsv2rgb(hsv-[0 0.5 0]);
    set(p,'Color',newCol);
    plot(xRand{ii},yRand{ii},'--','Color',[0.4 0.4 0.4])
end
h(1) = plot(mX,mY,'-k');
h(2) = plot(mRandX,mRandY,'--k');
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
ylabel('True Positive Rate')
xlabel('False Positive Rate')
legend(h,{['Real: ',num2str(round(mean(a),2))],...
    ['Permuted: ',num2str(round(mean(aRand),2))]},'Location','se')
title('Binge vs. Other')
%% PeriDrink
chan = [];
pair = []; %#ok
freq = 3;
loc = 'SLSR';
feat = 'a';
pow = zeros(44,32);
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper3\'...
    'periDrink\'],'.mat');
for ii = 1:size(files,2)
    load(files{ii})
    for jj = 1:size(psdTrls,2)
        if ~isempty(psdTrls{jj})
            pow(ii,jj) = mean(psdTrls{jj}.relPow(freq,chan,:));
        else
            pow(ii,jj) = NaN;
        end
    end
end
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drinkNotRaw.mat'])
mDrink = mean(catData{1,1}(:,(chan-1)*6+freq));
mNot = mean(catData{1,2}(:,(chan-1)*6+freq));
figure
hold on
% Pre
shadedErrorBar(15:16,fliplr(mean(pow(:,1:2),1,'omitnan').*100),...
    fliplr(std(pow(:,1:2),[],1,'omitnan').*100),{'color',[0 0.45 0.74]})
shadedErrorBar(12:15,fliplr(mean(pow(:,2:5),1,'omitnan').*100),...
    fliplr(std(pow(:,2:5),[],1,'omitnan').*100),'-b')
shadedErrorBar(1:12,fliplr(mean(pow(:,5:16),1,'omitnan').*100),...
    fliplr(std(pow(:,5:16),[],1,'omitnan').*100))
% Post
shadedErrorBar(17:18,mean(pow(:,17:18),1,'omitnan').*100,...
    std(pow(:,17:18),[],1,'omitnan').*100,{'color',[0 0.45 0.74]})
shadedErrorBar(18:21,mean(pow(:,18:21),1,'omitnan').*100,...
    std(pow(:,18:21),[],1,'omitnan').*100,'-b')
shadedErrorBar(21:32,mean(pow(:,21:32),1,'omitnan').*100,...
    std(pow(:,21:32),[],1,'omitnan').*100)
% Settings
set(gca,'XTick',[2:2:16,17:2:32],'XTickLabel',[-11.5:2:2.5,-2.5:2:12.5])
xlim([1 32])
xtickangle(90)
plot([1 32],[mDrink mDrink].*100,'--','color',[0 0.45 0.74])
plot([1 32],[mNot mNot].*100,'--k')
text(33,mDrink.*100,'Drink','color',[0 0.45 0.74])
text(33,mNot.*100,'Other')
xlabel('Time (s)')
ylabel(['% ',feat])
title([loc,feat])
%% PeriBinge
% Channels: SL,SR,CL,CR; SLSR,SLCL,SLCR,SRCL,SRCR,CLCR
% Freq: delta,theta,alpha,beta,lgamma,hgamma
clear chan pair freq
chan = [1]; %#ok
pair = [];
freq = 6;
feat = 'hg';
loc = 'SL';
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
   % Only grab up to [-56,-51]
   for t = 1:52
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
% Flip around
prePow = [fliplr(prePow(:,1:52)),prePow(:,53:end)];
preCoh = [fliplr(preCoh(:,1:52)),preCoh(:,53:end)];
% Get preBingeInfluence data
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    'preBingeInfluence'],'.mat');
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'psdTrls')
    else
        load(files{ii},'coh')
    end
    for t = 1:4
        if isempty(pair)
            prePow(ii,52+t) = mean(psdTrls{t}.relPow(freq,chan,:),...
            'omitnan');
        else
            preCoh(ii,52+t) = mean(coh{t}.normBandCoh(pair,freq,:),...
                'omitnan');
        end
    end
end
% Get Binge data
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper2\'...
    'binge'],'base','in');
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'psdTrls')
    else
        load(files{ii},'coh')
    end
    % Only grab up to [4,9]
    for t = 1:5
        if isempty(pair)
            prePow(ii,56+t) = mean(psdTrls{t}.relPow(freq,chan,:),...
                'omitnan');
        else
            preCoh(ii,56+t) = mean(coh{t}.rel(pair,freq,:),'omitnan');
        end
    end
end
if isempty(pair)
    mPre = mean(prePow,1,'omitnan');
    sPre = std(prePow,[],1,'omitnan');
else
    mPre = mean(preCoh,1,'omitnan');
    sPre = std(preCoh,[],1,'omitnan');
end
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\'...
    'bingeNotBingeTrial.mat'])
nameVect = names({'SL','SR','CL','CR'},{'d','t','a','b','lg','hg'});
featI = logicFind([loc,feat],nameVect,'==');
bingeCat = cat(1,data{1}{:,1});
notCat = cat(1,data{1}{:,2});
mNot = mean(notCat(:,featI));
mBinge = mean(bingeCat(:,featI));
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
% Pre and beginning of binge
shadedErrorBar(1:52,mPre(1:52).*100,sPre(1:52).*100)
shadedErrorBar(1:52,mPre(1:52).*100,sPre(1:52).*100)
shadedErrorBar(52:55,mPre(52:55).*100,sPre(52:55).*100,'-b')
shadedErrorBar(55:61,mPre(55:61).*100,sPre(55:61).*100,...
    {'color',[0 0.45 0.74]})
% End and post
shadedErrorBar(63:69,mPost(1:7).*100,sPost(1:7).*100,...
    {'color',[0 0.45 0.74]})
shadedErrorBar(69:72,mPost(7:10).*100,sPost(7:10).*100,'-b')
shadedErrorBar(72:123,mPost(10:61).*100,sPost(10:61).*100)
plot(1:123,ones(1,123).*mNot.*100,'--k')
plot(1:123,ones(1,123).*mBinge.*100,'--','color',[0 0.45 0.74])
% Settings
xlim([1 121])
set(gca,'XTick',[1:5:61,63:5:123],...
    'XTickLabel',[-53.5:5:6.5,-6.5:5:53.5])
xtickangle(90)
title([loc,' ',feat]);
ylabel(['% ',feat])
text(123,mean(mBinge)*100,'Binge','color',[0 0.45 0.74])
text(123,mNot*100,'Other')
xlabel('Time (s)')
box off
%% Find files with >20 drink samples
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\drinkNot\'],{'.mat'},{'pow','coh'},'trl','rel');
nDrink = cellfun(@(x) size(x,1),data{1}(:,1));
ind = logicFind(15,nDrink,'>=');
%% Calculate CFC for those animals
for ii = 1:size(ind,2)
    load(files{1}{ind(ii)},'trls')
    thisDrink = reshape(trls{1,1}.trial,4,2000*size(trls{1,1}.trial,3));
    thisNotDrink = reshape(trls{1,2}.trial(:,:,...
        randperm(size(trls{1,2}.trial,3),size(trls{1,1}.trial,3))),4,...
        2000*size(trls{1,1}.trial,3));
    [MIs,fSpace,Hvals,probVects,statsData] = gmwMI(thisDrink,400,1,100,...
        18,size(trls{1,1}.trial,3),2000,3,5,4);
    [notMIs,~,notHvals,notProbVects,notStatsData] = gmwMI(thisNotDrink,...
        400,1,100,18,size(trls{1,1}.trial,3),2000,3,5,4);
end
%%
for ii = 1:4
    figure
    subplot(1,2,1)
    imagesc(MIs(:,:,ii,ii)')
    fSpaceH = (fSpace)*(400/(2*pi));
    set(gca,'ytick',1:4:length(fSpace)),set(gca,'xtick',1:4:length(fSpace))
    set(gca,'xticklabel',round(fliplr(fSpaceH(1:4:length(fSpaceH)))',1))
    set(gca,'yticklabel',round((fSpaceH(1:4:length(fSpaceH)))',1))
    set(gca,'xdir','reverse')
    xlabel('Phase (Hz)')
    ylabel('Amplitude (Hz)')
    colormap('viridis')
    title('Drink')
    
    subplot(1,2,2)
    imagesc(notMIs(:,:,ii,ii)')
    fSpaceH = (fSpace)*(400/(2*pi));
    set(gca,'ytick',1:4:length(fSpace)),set(gca,'xtick',1:4:length(fSpace))
    set(gca,'xticklabel',round(fliplr(fSpaceH(1:4:length(fSpaceH)))',1))
    set(gca,'yticklabel',round((fSpaceH(1:4:length(fSpaceH)))',1))
    set(gca,'xdir','reverse')
    xlabel('Phase (Hz)')
    ylabel('Amplitude (Hz)')
    colormap('viridis')
    title('Not')
end