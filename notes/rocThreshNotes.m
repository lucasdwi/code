%% Go through datasets and recalculate ROC curve with thresholds determined 
% from all data - forces ROCs into same space
%% concatLog
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50.mat')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\concatLog\')
% Determine threshold - can just use 'concat' data since it includes all
% 'each' data
catThresh = [];
for ii = 1:20
    load([num2str(ii),'.mat'])
    catThresh = cat(1,catThresh,predict(concatData.model,all.testX{ii}(:,inds)));
end
catThresh = unique(catThresh);
% Check if 'catThresh' starts with zero or ends with 1, if not add them
if catThresh(1) ~= 0
    catThresh = [0;catThresh];
end
if catThresh(end) ~= 1
    catThresh = [catThresh;1];
end
% Apply threshold to data and predictions from model; save over old ROC
% data
for ii = 1:20
    load([num2str(ii),'.mat'])
    predY = predict(concatData.model,all.testX{ii}(:,inds));
    [x,y,~,a] = perfcurve(all.testY{ii},predY,1,'TVals',catThresh,'UseNearest',0);
    for jj = 1:12
        predY = predict(concatData.model,each.testX{ii,jj}(:,inds));
        [eachX{jj},eachY{jj},~,eachA(jj)] = perfcurve(each.testY{ii,jj},predY,1,'TVals',catThresh,'UseNearest',0);
    end
    concatData.rocX = x;
    concatData.rocY = y;
    concatData.auc = a;
    eachData.rocX = eachX;
    eachData.rocY = eachY;
    eachData.auc = eachA;
    save([num2str(ii),'.mat'],'eachData','concatData')
end
%% concatLogRand
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50.mat')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\concatLogRand\')
% Determine threshold - can just use 'concat' data since it includes all
% 'each' data
catThresh = [];
for ii = 1:20
    load([num2str(ii),'.mat'])
    catThresh = cat(1,catThresh,predict(concatData.model,all.testX{ii}(:,inds)));
end
catThresh = unique(catThresh);
% Check if 'catThresh' starts with zero or ends with 1, if not add them
if catThresh(1) ~= 0
    catThresh = [0;catThresh];
end
if catThresh(end) ~= 1
    catThresh = [catThresh;1];
end
% Apply threshold to data and predictions from model; save over old ROC
% data
for ii = 1:20
    load([num2str(ii),'.mat'])
    predY = predict(concatData.model,all.testX{ii}(:,inds));
    [x,y,~,a] = perfcurve(all.testY{ii},predY,1,'TVals',catThresh,'UseNearest',0);
    for jj = 1:12
        predY = predict(concatData.model,each.testX{ii,jj}(:,inds));
        [eachX{jj},eachY{jj},~,eachA(jj)] = perfcurve(each.testY{ii,jj},predY,1,'TVals',catThresh,'UseNearest',0);
    end
    concatData.rocX = x;
    concatData.rocY = y;
    concatData.auc = a;
    eachData.rocX = eachX;
    eachData.rocY = eachY;
    eachData.auc = eachA;
    save([num2str(ii),'.mat'],'eachData','concatData')
end
%% eachLog
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50.mat')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\eachLog\')
% Determine threshold - can just use 'concat' data since it includes all
% 'each' data
catThresh = [];
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    catThresh = cat(1,catThresh,predict(selfData.model,all.testX{iter}(:,inds)));
end
catThresh = unique(catThresh);
% Check if 'catThresh' starts with zero or ends with 1, if not add them
if catThresh(1) ~= 0
    catThresh = [0;catThresh];
end
if catThresh(end) ~= 1
    catThresh = [catThresh;1];
end
% Apply threshold to data and predictions from model; save over old ROC
% data
for ii = 1:240
    disp(ii)
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    for jj = 1:12
        if jj ~= animal
            predY = predict(selfData.model,each.testX{iter,jj}(:,inds));
            [eachX{jj},eachY{jj},~,eachA(jj)] = perfcurve(each.testY{iter,jj},predY,1,'TVals',catThresh,'UseNearest',0);
        else
            predY = predict(selfData.model,each.testX{iter,jj}(:,inds));
            [x,y,~,a] = perfcurve(each.testY{iter,jj},predY,1,'TVals',catThresh,'UseNearest',0);
        end
    end
    predY = predict(selfData.model,all.testX{iter}(:,inds));
    [concatX,concatY,~,concatA] = perfcurve(all.testY{iter},predY,1,'TVals',catThresh,'UseNearest',0);
    selfData.rocX = x;
    selfData.rocY = y;
    selfData.auc = a;
    concatData.rocX = concatX;
    concatData.rocY = concatY;
    concatData.auc = concatA;
    eachData.rocX = eachX;
    eachData.rocY = eachY;
    eachData.auc = eachA;
    save([num2str(ii),'.mat'],'selfData','eachData','concatData')
end
%% eachLogRand
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50.mat')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\eachLogRand\')
% Determine threshold - can just use 'concat' data since it includes all
% 'each' data
catThresh = [];
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    catThresh = cat(1,catThresh,predict(selfData.model,all.testX{iter}(:,inds)));
end
catThresh = unique(catThresh);
% Check if 'catThresh' starts with zero or ends with 1, if not add them
if catThresh(1) ~= 0
    catThresh = [0;catThresh];
end
if catThresh(end) ~= 1
    catThresh = [catThresh;1];
end
% Apply threshold to data and predictions from model; save over old ROC
% data
for ii = 1:240
    disp(ii)
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    for jj = 1:12
        if jj ~= animal
            predY = predict(selfData.model,each.testX{iter,jj}(:,inds));
            [eachX{jj},eachY{jj},~,eachA(jj)] = perfcurve(each.testY{iter,jj},predY,1,'TVals',catThresh,'UseNearest',0);
        else
            predY = predict(selfData.model,each.testX{iter,jj}(:,inds));
            [x,y,~,a] = perfcurve(each.testY{iter,jj},predY,1,'TVals',catThresh,'UseNearest',0);
        end
    end
    predY = predict(selfData.model,all.testX{iter}(:,inds));
    [concatX,concatY,~,concatA] = perfcurve(all.testY{iter},predY,1,'TVals',catThresh,'UseNearest',0);
    selfData.rocX = x;
    selfData.rocY = y;
    selfData.auc = a;
    concatData.rocX = concatX;
    concatData.rocY = concatY;
    concatData.auc = concatA;
    eachData.rocX = eachX;
    eachData.rocY = eachY;
    eachData.auc = eachA;
    save([num2str(ii),'.mat'],'selfData','eachData','concatData')
end
%% each
load('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50.mat')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
cd('C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\each\')
% Determine threshold - can just use 'concat' data since it includes all
% 'each' data
catThresh = [];
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    catThresh = cat(1,catThresh,cvglmnetPredict(selfData.model,all.testX{iter}(:,inds)));
end
catThresh = unique(catThresh);
% Check if 'catThresh' starts with zero or ends with 1, if not add them
if catThresh(1) ~= 0
    catThresh = [0;catThresh];
end
if catThresh(end) ~= 1
    catThresh = [catThresh;1];
end
% Apply threshold to data and predictions from model; save over old ROC
% data
for ii = 1:240
    disp(ii)
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    for jj = 1:12
        if jj ~= animal
            predY = predict(selfData.model,each.testX{iter,jj}(:,inds));
            [eachX{jj},eachY{jj},~,eachA(jj)] = perfcurve(each.testY{iter,jj},predY,1,'TVals',catThresh,'UseNearest',0);
        else
            predY = predict(selfData.model,each.testX{iter,jj}(:,inds));
            [x,y,~,a] = perfcurve(each.testY{iter,jj},predY,1,'TVals',catThresh,'UseNearest',0);
        end
    end
    predY = predict(selfData.model,all.testX{iter}(:,inds));
    [concatX,concatY,~,concatA] = perfcurve(all.testY{iter},predY,1,'TVals',catThresh,'UseNearest',0);
    selfData.rocX = x;
    selfData.rocY = y;
    selfData.auc = a;
    concatData.rocX = concatX;
    concatData.rocY = concatY;
    concatData.auc = concatA;
    eachData.rocX = eachX;
    eachData.rocY = eachY;
    eachData.auc = eachA;
    save([num2str(ii),'.mat'],'selfData','eachData','concatData')
end
%% pre
catThresh = [];
for ii = 1:20
    load(['C:\Users\Pythia\documents\GreenLab\data\paper2\analyzed\finalNew\preBinge\',num2str(ii),'.mat'])
    catThresh = cat(1,catThresh,concatData{1,1}.acc{1,1}.pred');
end
catThresh = unique(catThresh);
% Check if 'catThresh' starts with zero or ends with 1, if not add them
if catThresh(1) ~= 0
    catThresh = [0;catThresh];
end
if catThresh(end) ~= 1
    catThresh = [catThresh;1];
end

load('C:\Users\Pythia\documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50.mat','pInds','all')
load('C:\Users\Pythia\documents\GreenLab\data\paper2\analyzed\finalNew\preBingeTrainTest.mat')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
for ii = 1:20
    load(['C:\Users\Pythia\documents\GreenLab\data\paper2\analyzed\finalNew\preBinge\',num2str(ii),'.mat'])
    notFeedTest = all.testX{ii}(all.testY{ii}==0,inds);
    numPre = round((size(notFeedTest,1)-(1-0.006).*size(notFeedTest,1))./(1-0.006));
    testY = [ones(numPre,1);zeros(size(notFeedTest,1),1)];
    [concatData{1,1}.acc{1}.x,concatData{1,1}.acc{1,1}.y,~,~] = perfcurve(testY,concatData{1,1}.acc{1,1}.pred,1,'Tvals',catThresh,'UseNearest',0);
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\preBinge\',num2str(ii),'.mat'],'concatData')
end
%% preRand
catThresh = [];
for ii = 1:20
    load(['C:\Users\Pythia\documents\GreenLab\data\paper2\analyzed\finalNew\preBingeRand\',num2str(ii),'.mat'])
    catThresh = cat(1,catThresh,concatData{1,1}.acc{1,1}.pred');
end
catThresh = unique(catThresh);
% Check if 'catThresh' starts with zero or ends with 1, if not add them
if catThresh(1) ~= 0
    catThresh = [0;catThresh];
end
if catThresh(end) ~= 1
    catThresh = [catThresh;1];
end

load('C:\Users\Pythia\documents\GreenLab\data\paper2\analyzed\finalNew\baseline500Each6000All50-50.mat','pInds','all')
load('C:\Users\Pythia\documents\GreenLab\data\paper2\analyzed\finalNew\preBingeTrainTest.mat')
inds = 1:60;
inds = inds(~ismember(inds,pInds));
for ii = 1:20
    load(['C:\Users\Pythia\documents\GreenLab\data\paper2\analyzed\finalNew\preBingeRand\',num2str(ii),'.mat'])
    notFeedTest = all.testX{ii}(all.testY{ii}==0,inds);
    numPre = round((size(notFeedTest,1)-(1-0.006).*size(notFeedTest,1))./(1-0.006));
    testY = [ones(numPre,1);zeros(size(notFeedTest,1),1)];
    [concatData{1,1}.acc{1}.x,concatData{1,1}.acc{1,1}.y,~,~] = perfcurve(testY,concatData{1,1}.acc{1,1}.pred,1,'Tvals',catThresh,'UseNearest',0);
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\preBingeRand\',num2str(ii),'.mat'],'concatData')
end