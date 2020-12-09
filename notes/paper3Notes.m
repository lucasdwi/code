
%% Feeding
load('D:/paper3/preBingeAllData.mat')
allNot = cat(1,notBingeCat{:});
for ii = 1:100
    load(['D:/paper3/analyzed/preFeed/preFeed',num2str(ii),'.mat'],'acc','accR')
    for jj = 1:size(preData,2)
        thisPast = cat(1,preData{:,jj});
        thisCatX = [thisPast;allNot(randperm(size(allNot,1),size(thisPast,1)),:)];
        thisCatY = [ones(size(thisPast,1),1);zeros(size(thisPast,1),1)];
        prob = cvglmnetPredict(acc{1}.mdl{1},thisCatX,'lambda_1se','response');
        [~,~,~,aLassoFeed(ii,jj)] = perfcurve(thisCatY,prob,1);
        prob = cvglmnetPredict(acc{1}.mdl{1},thisCatX(randperm(size(thisCatX,1),size(thisCatX,1)),:),'lambda_1se','response');
        [~,~,~,aLassoFeedP(ii,jj)] = perfcurve(thisCatY,prob,1);
    end
    a5(ii) = acc{1}.acc;
    a5p(ii) = accR{1}.acc;
end
%% Plot
figure
hold on
h2 = shadedErrorBar(1:241,[mean(a5p),mean(aLassoFeedP,1)],[std(a5p),std(aLassoFeedP)],'k');
h1 = shadedErrorBar(1:241,[mean(a5),mean(aLassoFeed,1)],[std(a5),std(aLassoFeed)],'r');
set(gca,'xtick',2:20:241,'xticklabel',-2.5:-20:-242.5)
legend([h2.mainLine,h1.mainLine],{'all','perm'})
xlabel('time before feeding')
ylabel('auc')


%% Male vs. Female
% Preallocate
[ffA,mmA] = deal(zeros(1,5));
[fmA,mfA] = deal(zeros(5,5));
for ii = 1:5
   load(['D:\paper3\analyzed\maleFemale\maleFemale',num2str(ii),'.mat']) 
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
%%

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
% Build linear mixed model with group as random effect
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
    
    mdl = fitglm(allData(trainInds,:),formula(1:end-7),'Distribution',...
        'binomial','binomialsize',nTrain);
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
    
    mdl = fitglm(subData(trainInds,:),formula(1:end-7),'Distribution',...
        'binomial','binomialsize',nTrain);
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
[dataNot,samp,files] = collateData(['D:\paper3\'...
    'waterAlcohol\notDrink\'],{'.mat'},{'pow','coh'},'trl','rel');
% Get animal ID and day ID
n = cell2mat(cellfun(@(x) size(x,1),dataNot{1,1}(:,1),'uniformoutput',0));
animal = []; day = [];
for ii = 1:size(dataNot{1},1)
    parts = strsplit(files{1}{ii},'_');
    animal{ii,1} = repmat({[parts{1},parts{2}]},n(ii,1),1);
    day{ii,1} = repmat({parts{4}},n(ii,1),1);
end
% Concatenate and z-score alcohol data
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
%% Build model alcohol vs. water/notDrik
load('D:\paper3\waterAlcoholNotData.mat')
% Count number of trials for each group from each animal
ids = unique(allData.ID);
for ii = 1:size(ids,1)
    nSamps(ii,1) = sum(allData.group(strcmp(allData.ID,ids{ii}))==1);
    nSamps(ii,2) = sum(allData.group(strcmp(allData.ID,ids{ii}))==0);
end
% Double water samples to account for water and not data being combined
nSamps(:,2) = nSamps(:,2)*2;
% Go through and find the animals with at least 1 samp in each group
for ii = 1:size(nSamps)
    both(ii) = nSamps(ii,1)>0 && nSamps(ii,2)>0;
end
bothIDs = logicFind(1,both,'==');
minSamp = min(nSamps(bothIDs,:),[],2);
minSamp(isodd(minSamp)) = minSamp(isodd(minSamp))-1;
% Subset data to include only animals with >=5 samples
subData = [];
for ii = 1:size(bothIDs,2)
    if minSamp(ii)>5
    thisWat = allData(strcmp(allData.ID,ids(bothIDs(ii))) & ...
        allData.group==0,:);
    thisNot = notD(strcmp(notD.ID,ids(bothIDs(ii))),:);
    thisAlc = allData(strcmp(allData.ID,ids(bothIDs(ii))) & ...
        allData.group==1,:);
    % Grab 
    subData = [subData;thisWat(randperm(size(thisWat,1),minSamp(ii)/2),:);...
        thisNot(randperm(size(thisNot,1),minSamp(ii)/2),:);...
        thisAlc(randperm(size(thisAlc,1),minSamp(ii)),:)];
    end
end

for ii = 1:100
    load(['D:\paper3\analyzed\alcoholOther\alcoholOther',num2str(ii),'.mat'])
    [alcRocX(ii,:),alcRocY(ii,:),~,alcA(ii)] = perfcurve(hist.cfg.naive.testY,acc{1}.pred,1);
    alcOtherA(ii) = acc{1}.acc;
    [thisX,thisY,~,alcRA(ii)] = perfcurve(histR.cfg.naive.testY,accR{1}.pred,1);
    if numel(thisX) == 2
        alcRocXR(ii,:) = linspace(0,1,107);
        alcRocYR(ii,:) = linspace(0,1,107);
    else
        alcRocXR(ii,:) = thisX;
        alcRocYR(ii,:) = thisY;
    end
    alcOtherR(ii) = accR{1}.acc;
end
figure
hold on
plot(mean(alcRocX,1),mean(alcRocY,1),'-k')
plot(mean(alcRocXR,1),mean(alcRocYR,1),'--k')
legend({['Real: ',num2str(round(mean(alcA),2)),'\pm',...
    num2str(round(conf(alcA,0.95),2))],...
    ['Permuted: ',num2str(round(mean(alcRA),2)),'\pm',...
    num2str(round(conf(alcRA,0.95),2))]},'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR');
title('Alcohol vs. Water/notDrink')
%% Build model sweet-fat vs. chow/notBinge
load('D:/paper2/analyzed/finalNew/bingeNotPreData.mat')
samps(:,1) = cellfun(@(x) size(x,1),data{1,1}([1:3,5:8,10,11],1));
samps(:,2) = cellfun(@(x) size(x,1),data{1,4}(:,1))*2;
minSamp = min(samps,[],2);
minSamp(isodd(minSamp)) = minSamp(isodd(minSamp))-1;
x = []; y = [];
bingeInds = [1:3,5:8,10,11];
for ii = 1:9
    thisBinge = data{1,1}{bingeInds(ii),1};
    thisChow = data{1,4}{ii,1};
    thisNot = [data{1,1}{bingeInds(ii),2};data{1,4}{ii,1}];
    % Grab
    x = [x;thisBinge(randperm(size(thisBinge,1),minSamp(ii)),:);...
        thisChow(randperm(size(thisChow,1),minSamp(ii)/2),:);...
        thisNot(randperm(size(thisNot,1),minSamp(ii)/2),:)];
    y = [y;ones(minSamp(ii),1);zeros(minSamp(ii),1)];
end
for ii = 1:100
    load(['D:\paper3\analyzed\bingeOther\bingeOther',num2str(ii),'.mat'])
    [bingeRocX(ii,:),bingeRocY(ii,:),~,bingeA(ii)] = perfcurve(hist.cfg.naive.testY,acc{1}.pred,1,'Tvals',linspace(0,1,273),'usenearest',0);
    bingeOtherA(ii) = acc{1}.acc;
    [thisX,thisY,~,bingeRA(ii)] = perfcurve(histR.cfg.naive.testY,accR{1}.pred,1,'Tvals',linspace(0,1,273),'usenearest',0);
    if numel(thisX) == 2
        bingeRocXR(ii,:) = linspace(0,1,273);
        bingeRocYR(ii,:) = linspace(0,1,273);
    else
        bingeRocXR(ii,:) = thisX;
        bingeRocYR(ii,:) = thisY;
    end
    bingeOtherR(ii) = accR{1}.acc;
end
figure
hold on
plot(mean(bingeRocX,1),mean(bingeRocY,1),'-k')
plot(mean(bingeRocXR,1),mean(bingeRocYR,1),'--k')
legend({['Real: ',num2str(round(mean(bingeA),2)),'\pm',...
    num2str(round(conf(bingeA,0.95),4))],...
    ['Permuted: ',num2str(round(mean(bingeRA),2)),'\pm',...
    num2str(round(conf(bingeRA,0.95),2))]},'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR');
title('Sweet-Fat vs. Chow/notEat')
%% Find trade-off of number of animals and minimum of samples per behavior
load('waterAlcoholNotData.mat')
ids = unique(allData.ID);
for ii = 1:size(ids,1)
    nSamps(ii,1) = sum(allData.group(strcmp(allData.ID,ids{ii}))==1);
    nSamps(ii,2) = sum(allData.group(strcmp(allData.ID,ids{ii}))==0)*2;
end
for ii = 1:size(nSamps)
    both(ii) = nSamps(ii,1)>0 && nSamps(ii,2)>0;
end
nSamps = nSamps(both,:);
load('D:\paper2\analyzed\finalNew\bingeNotPreData.mat')
sweetSamps = cellfun(@(x) size(x,1),data{2}([1:3,5:8,10,11],:));
chowSamps = cellfun(@(x) size(x,1),data{4});
bingeSamps(:,1) = sweetSamps(:,1);
bingeSamps(:,2) = chowSamps(:,1)*2;
% bingeSamps(:,3) = sweetSamps(:,2)+chowSamps(:,2);
all = [nSamps;bingeSamps];
minBinge = min(bingeSamps,[],2);
minAlc = min(nSamps,[],2);
minSamps = sort(min(all,[],2),'ascend');
for ii = 1:size(minSamps,1)
    bingeN(ii) = sum(minBinge>=minSamps(ii));
    alcN(ii) = sum(minAlc>=minSamps(ii));
    total(ii) = minSamps(ii)*sum(minSamps>=minSamps(ii));
end
figure
hold on
yyaxis left
plot(minSamps,bingeN,'-o')
plot(minSamps,alcN,'-o')
yyaxis right
plot(minSamps,total,'-o')
%% Build best possible drinking model
load('D:/paper3/waterAlcoholNotData2.mat')
ids = unique([unique(notD.ID);unique(water.ID);unique(alc.ID);...
    unique(watAlc.ID)]);
samps = [];
for ii = 1:size(ids,1)
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
for ii = 1:size(theseIDs,1)
    thisAlc{ii,1} = [watAlc(strcmp(watAlc.ID,theseIDs{ii}) & ...
        watAlc.group == 1,:); alc(strcmp(alc.ID,theseIDs{ii}) & ...
        alc.group == 1,:)];
    thisWat{ii,1} = [watAlc(strcmp(watAlc.ID,theseIDs{ii}) & ...
        watAlc.group == 0,:); water(strcmp(water.ID,theseIDs{ii}) ...
        & water.group == 1,:)];
    thisNeither{ii,1} = notD(strcmp(notD.ID,theseIDs{ii}) & ...
        notD.group == 0,:);
end
allData = [thisAlc,thisWat,thisNeither];
% Set target number of samples for alcohol (water and neither will be 1/2)
n = 20;
target = [n;n/2;n/2];
these = []; theseWeights = [];
for ii = 1:size(allData,1)
    for jj = 1:3
        % If fewer samples exist than target, then pull all and compute
        % weights
        if theseSamps(ii,jj)<target(jj)
           these = [these;allData{ii,jj}];
           theseWeights = [theseWeights;repmat(target(jj)/theseSamps(ii,jj),...
               theseSamps(ii,jj),1)];
        % Otherwise, grab random subsample
        else
            these = [these;allData{ii,jj}(randperm(theseSamps(ii,jj),...
                target(jj)),:)];
            theseWeights = [theseWeights;ones(target(jj),1)];
        end
    end
end
% Set up training and testing set indices (80/20 random split)
[trainInd,~,testInd,~] = trainTest((1:numel(theseWeights))',...
    (1:numel(theseWeights))',.20);
% Use indices to create training and testing sets, and split weights
trainX = these(trainInd,1:60);
trainY = these.group(trainInd);
trainWeights = theseWeights(trainInd);
testX = these(testInd,1:60);
testY = these.group(testInd); 
testWeights = theseWeights(testInd);
% Build models - only use shell features (to match feeding features)
cfg = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se',trainWeights);
[~,lam,beta,fits,acc,hist] = lassoNet(trainX,trainY,'binomial','class',...
    1,10,1,cfg);
cfg = lassoNetCfg({testX,testY},[],'y','y','n',100,'1se',trainWeights);
[~,lamR,betaR,fitsR,accR,histR] = lassoNet(trainX,trainY,'binomial',...
    'class',1,10,1,cfg);
%% Load drinking model data
for ii = 1:100
%     % 40 sample model
%     cd D:\paper3\analyzed\alcoholOther40\
%     load(['alcoholOther',num2str(ii),'.mat'])
%     acc40(ii) = acc{1}.acc;
%     accR40(ii) = accR{1}.acc;
%     % 80 sample model
%     cd D:\paper3\analyzed\alcoholOther80\
%     load(['alcoholOther',num2str(ii),'.mat'])
%     acc80(ii) = acc{1}.acc;
%     accR80(ii) = accR{1}.acc;
    % 120 sample model
    cd D:\paper3\analyzed\alcoholOther120\
    load(['alcoholOther',num2str(ii),'.mat'])
    acc120(ii) = acc{1}.acc;
    accR120(ii) = accR{1}.acc;
%     % 160 sample model
%     cd D:\paper3\analyzed\alcoholOther160\
%     load(['alcoholOther',num2str(ii),'.mat'])
%     acc160(ii) = acc{1}.acc;
%     accR160(ii) = accR{1}.acc;
%     % 200 sample model
%     cd D:\paper3\analyzed\alcoholOther200\
%     load(['alcoholOther',num2str(ii),'.mat'])
%     acc200(ii) = acc{1}.acc;
%     accR200(ii) = accR{1}.acc;
end
%% Plot performance curve
figure
hold on
plot([40,80,120,160,200],[mean(acc40),mean(acc80),mean(acc120),mean(acc160),mean(acc200)],'-o')
plot([40,80,120,160,200],[mean(accR40),mean(accR80),mean(accR120),mean(accR160),mean(accR200)],'-o')
xlabel('Samples per animal')
ylabel('Mean AUC')
legend({'Real','Permuted'})
%% Load alcohol models by sex
cd D:\paper3\analyzed\alcoholOtherSex
for ii = 1:100
    load(['alcoholOtherMale',num2str(ii),'.mat'],'accMale','accMaleR',...
        'histMale','histMaleR')
    predBinary = round(accMale{1}.pred);
    maleA(ii) = sum(predBinary==histMale.cfg.naive.testY')/numel(predBinary);
    predBinaryR = round(accMaleR{1}.pred);
    maleAR(ii) = sum(predBinaryR==histMaleR.cfg.naive.testY')/numel(predBinaryR);
    [maleX(ii,:),maleY(ii,:),~,maleAcc(ii)] = perfcurve(histMale.cfg.naive.testY,...
        accMale{1}.pred,1,'TVals',linspace(0,1,53),'UseNearest',0);
    if numel(thisX)~=54
        maleX(ii,:) = interp1(thisX.*(1+eps),thisY.*(1+eps),linspace(0,1,54));       
        maleY(ii,:) = interp1(thisY,thisX,linspace(0,1,54));
    else
        maleX(ii,:) = thisX;
        maleY(ii,:) = thisY;
    end
    [maleXR(ii,:),maleYR(ii,:),~,maleAccR(ii)] = perfcurve(...
        histMaleR.cfg.naive.testY,accMaleR{1}.pred,1,'TVals',...
        linspace(0,1,53),'UseNearest',0);
    if numel(thisXR)~=54
        maleXR(ii,:) = interp1(thisXR,thisYR,linspace(0,1,54));       
        maleYR(ii,:) = interp1(thisYR,thisXR,linspace(0,1,54));
    else
        maleXR(ii,:) = thisXR;
        maleYR(ii,:) = thisYR;
    end
    
    load(['alcoholOtherFemale',num2str(ii),'.mat'],'accFemale',...
        'accFemaleR','histFemale','histFemaleR')
    predBinary = round(accFemale{1}.pred);
    femaleA(ii) = sum(predBinary==histFemale.cfg.naive.testY')/numel(predBinary);
    predBinaryR = round(accFemaleR{1}.pred);
    femaleAR(ii) = sum(predBinaryR==histFemaleR.cfg.naive.testY')/numel(predBinaryR);
    [femaleX(ii,:),femaleY(ii,:),~,femaleAcc(ii)] = perfcurve(...
        histFemale.cfg.naive.testY,accFemale{1}.pred,1,'TVals',...
        linspace(0,1,53),'UseNearest',0);
    if numel(thisX)~=54
        femaleX(ii,:) = interp1(thisX,thisY,linspace(0,1,54));      
        femaleY(ii,:) = interp1(thisY,thisX,linspace(0,1,54));
    else
        femaleX(ii,:) = thisX;
        femaleY(ii,:) = thisY;
    end
    [femaleXR(ii,:),femaleYR(ii,:),~,femaleAccR(ii)] = perfcurve(...
        histFemaleR.cfg.naive.testY,accFemaleR{1}.pred,1,'TVals',...
        linspace(0,1,53),'UseNearest',0);
     if numel(thisXR)~=54
        femaleXR(ii,:) = interp1(thisXR,thisYR,linspace(0,1,54));       
        femaleYR(ii,:) = interp1(thisYR,thisXR,linspace(0,1,54));
     else
         femaleXR(ii,:) = thisXR;
         femaleYR(ii,:) = thisYR;
    end
end
%
figure
hold on
plot(mean(maleX,1),mean(maleY,1),'-k')
plot(mean(maleXR,1),mean(maleYR,1),'--k')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('Male (n=4): Alcohol vs. Water + Other')
legend({['Real: ',num2str(round(mean(maleAcc),2)),'\pm',...
    num2str(round(conf(maleAcc,0.95),2))],['Permuted: ',...
    num2str(round(mean(maleAccR),2)),'\pm',...
    num2str(round(conf(maleAccR,0.95),2))]},'location','nw')
figure
hold on
plot(mean(femaleX,1),mean(femaleY,1),'-k')
plot(mean(femaleXR,1),mean(femaleYR,1),'--k')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('Female (n=5): Alcohol vs. Water + Other')
legend({['Real: ',num2str(round(mean(femaleAcc),2)),'\pm',...
    num2str(round(conf(femaleAcc,0.95),2))],['Permuted: ',...
    num2str(round(mean(femaleAccR),2)),'\pm',...
    num2str(round(conf(femaleAccR,0.95),2))]},'location','nw')
%% Build best possible feeding model
load('D:\paper2\analyzed\finalNew\bingeNotPreData.mat')
% Set up indices of the nine animals found in both 24 sweet dep and 24 chow
nineInds = [1:3,5:8,10,11];
sweetSamps = cellfun(@(x) size(x,1),data{2}(nineInds,:));
chowSamps = cellfun(@(x) size(x,1),data{4});
bingeSamps(:,1) = sweetSamps(:,1);
bingeSamps(:,2) = chowSamps(:,1);
bingeSamps(:,3) = sweetSamps(:,2)+chowSamps(:,2);
% Combine data into similar structure as drinking data
allFeedData(:,1) = data{2}(nineInds,1);
allFeedData(:,2) = data{4}(:,1);
for ii = 1:9
    allFeedData{ii,3} = [data{2}{nineInds(ii),2};data{4}{ii,2}];
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
% save('D:\paper3\24sweetChowNotData.mat','allFeedData','allY','bingeSamps')
% Use same model building process as drinking models
% Set target number of samples for 24 hour sweet dep (chow and neither will
% be 1/2)
n = 20;
target = [n;n/2;n/2];
these = []; theseWeights = []; theseY = [];
for ii = 1:size(allFeedData,1)
    for jj = 1:3
        % If fewer samples exist than target, then pull all and compute
        % weights
        if bingeSamps(ii,jj)<target(jj)
           these = [these;allFeedData{ii,jj}];
           theseWeights = [theseWeights;repmat(target(jj)/...
               bingeSamps(ii,jj),bingeSamps(ii,jj),1)];
           % Add to y vector
           theseY = [theseY;allY{ii,jj}];
        % Otherwise, grab random subsample
        else
            % Generate indices
            theseInds = randperm(bingeSamps(ii,jj),target(jj));
            these = [these;allFeedData{ii,jj}(theseInds,:)];
            theseWeights = [theseWeights;ones(target(jj),1)];
            % Add to y vector
            theseY = [theseY;allY{ii,jj}(theseInds)];
        end
    end
end
% Set up training and testing set indices (80/20 random split)
[trainInd,~,testInd,~] = trainTest((1:numel(theseWeights))',...
    (1:numel(theseWeights))',.20);
% Use indices to create training and testing sets, and split weights
trainX = these(trainInd,1:60);
trainY = theseY(trainInd);
trainWeights = theseWeights(trainInd);
testX = these(testInd,1:60);
testY = theseY(testInd); 
testWeights = theseWeights(testInd);
% Build models - only use shell features (to match feeding features)
cfg = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se',trainWeights);
[~,lam,beta,fits,acc,hist] = lassoNet(trainX,trainY,'binomial','class',...
    1,10,1,cfg);
cfg = lassoNetCfg({testX,testY},[],'y','y','n',100,'1se',trainWeights);
[~,lamR,betaR,fitsR,accR,histR] = lassoNet(trainX,trainY,'binomial',...
    'class',1,10,1,cfg);
%% Load feeding model data
allN = 20:40:380;
for ii = 1:10
    cd(['D:/paper3/analyzed/24bingeOther/24bingeOther',num2str(allN(ii))])
    for jj = 1:100
        load(['24bingeOther',num2str(jj),'.mat'])
        feedAcc(ii,jj) = acc{1}.acc;
        feedAccR(ii,jj) = accR{1}.acc;
    end
end
figure
hold on
plot(allN,mean(feedAcc,2),'-o')
plot(allN,mean(feedAccR,2),'-o')

%% Plot features from each data set
% Get sign data
feedSign = double(mean(feedBeta,2)>=0);
feedSign(feedSign==0) = -1;
drinkSign = double(mean(drinkBeta,2)>=0);
drinkSign(drinkSign==0) = -1;
% Load gen model to get gen single feature AUCs
% theseGenA = [];
% for ii = 1:100
%     load(['D:/paper3/analyzed/genModel/1feedWeight/feedDrinkOther',num2str(ii),'.mat'],'sA')
%     theseGenA(ii,:) = sA;
% end
% meanGenA = mean(theseGenA,1);
% [sortMeanGenA,sortInd] = sort(meanGenA,'descend'); 
figure
hold on
scatter(mean(feedA(1:12,:),2).*feedSign(1:12),mean(drinkA(1:12,:),2).*drinkSign(1:12),'ks')
scatter(mean(feedA(13:18,:),2).*feedSign(13:18),mean(drinkA(13:18,:),2).*drinkSign(13:18),'ro')
xlabel('Feed AUC'); ylabel('Drink AUC')
xlim([-0.8 0.8]); ylim([-0.8 0.8])
title('single feature AUC, built within dataset')
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
% Only keep 10 water of animal 4
data{4,1} = data{4,1}(1:20,:);
animal{4,1} = animal{4,1}(1:20,:);

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
[all,~,rnd,~] = evenDataSplit(catData,18494,4624,'ADA',20); %#ok
save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'concatData.mat'],'all','rnd')
% Or force into same sample size as bingeNot data
[all,~,rnd,~] = evenDataSplit(catData,6000,1200,'ADA',20);
save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000.mat'],'all','rnd')
%% Build drinkNot models using drink6000.mat for comparison to bingeNot
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000.mat'])
% load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
%     'concatData.mat'])
% Preallocate
[x,y,xRand,yRand] = deal(cell(1,20));
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
%% Build models: binge, drink, and across
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'])
binge = all;
% bingeRnd = rnd;
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000.mat'],'all','rnd')
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
    ddYperm,dbX,dbY,dbXperm,dbYperm] = deal(cell(1,20));
[bbA,bbAperm,bdA,bdAperm,ddA,ddAperm,dbA,dbAperm] = deal(zeros(1,20));
for ii = 1:20
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
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50.mat'])
binge = all;
bingeRnd = rnd; %#ok
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000.mat'],'all','rnd')
drink = all;
drinkRnd = rnd; %#ok
bingeInds = [1:12,25:30];
drinkInds = [13:24,55:60];
% Set step for TVals
tStep = 1/2400;
% Preallocate
[ggX,ggY,ggXperm,ggYperm,gbX,gbY,gbXperm,gbYperm,gdX,gdY,gdXperm,...
    gdYperm] = deal(cell(1,20));
[ggA,ggAperm,gbA,gbAperm,gdA,gdAperm] = deal(zeros(1,20));
betas = zeros(18,20);
for ii = 1:20
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
load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\baseline500Each6000All50-50v2.mat'])
binge = all;
bingeRnd = rnd;
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'drink6000.mat'],'all','rnd')
drink = all;
drinkRnd = rnd;
bingeInds = [1:12,25:30];
drinkInds = [13:24,55:60];
% Preallocate
[bingeA,drinkA] = deal(cell(1,3));
for k = 1%:3
    cmbs = nchoosek(1:18,k);
    for ii = 1:20
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
genAM = mean(genA{1},1)';
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