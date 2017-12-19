%% Concat
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat50-50_500TrainNewRand.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\pInds.mat','pInds')
% Get usable indices
inds = 1:60;
inds = inds(~ismember(inds,pInds));
for ii = 1:20
    tic
    trainX = allTrainX{1,ii}(:,inds);
    trainY = allTrainY{1,ii};
    mdl = fitcsvm(trainX,trainY,'OptimizeHyperparameters','auto','HyperparameterOptimization',struct('kfold',10,'Verbose',0,'ShowPlots',0));
    scoreMdl{ii} = fitPosterior(mdl);
    testX = allTestX{1,ii}(:,inds);
    testY = allTestY{1,ii};
    [~,prob] = predict(scoreMdl{ii},testX);
    [concatX{ii},concatY{ii},~,concatA(ii)] = perfcurve(testY,prob(:,2),1);
    for jj = 1:12
        testX = [];
        testX = eachTestX{ii,jj}(:,inds);
        testY = [];
        testY = eachTestY{ii,jj};
        [~,prob] = predict(scoreMdl{ii},testX);
        [eachX{ii,jj},eachY{ii,jj},~,eachA(ii,jj)] = perfcurve(testY,prob(:,2),1);
    end
    toc
end
%% Regularized - LASSO
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat50-50_500TrainNew.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeNotData.mat','pInds')
% Build and test models from concat to concat
lambda = logspace(-5,-1,15);
% Get usable indices
inds = 1:60;
inds = inds(~ismember(inds,pInds));
for ii = 1:20
    tic
    disp(num2str(ii))
    trainX = allTrainX{1,ii}(:,inds);
    trainY = allTrainY{1,ii};
    CVmdl = fitrlinear(trainX,trainY,'kfold',10,'lambda',lambda,'Regularization','lasso');
    mse = kfoldLoss(CVmdl);
    ind = logicFind(min(mse),mse,'==');
    lam = lambda(ind);
    mdl{ii} = fitrlinear(trainX,trainY,'lambda',lam,'Regularization','lasso');
    testX = allTestX{1,ii}(:,inds);
    testY = allTestY{1,ii};
    pred = predict(mdl{ii},testX);
    [concatX{ii},concatY{ii},~,concatA(ii)] = perfcurve(testY,pred,1);
    % Test concat models on each data
    for jj = 1:12
        testX = [];
        testX = eachTestX{ii,jj}(:,inds);
        testY = [];
        testY = eachTestY{ii,jj};
        pred = predict(mdl{ii},testX);
        [eachX{ii,jj},eachY{ii,jj},~,eachA(ii,jj)] = perfcurve(testY,pred,1);
    end
    toc
end
%%
cMX = mean(cat(2,concatX{:}),2);
cMY = mean(cat(2,concatY{:}),2);
cSX = std(cat(2,concatX{:}),[],2);
cSY = std(cat(2,concatY{:}),[],2);
cEXp = cMX+cSX;
cEXm = cMX-cSX;
cEYp = cMY+cSY;
cEYm = cMY-cSY;
cXFill = [cEXp;flipud(cEXm)];
cYFill = [cEYp;flipud(cEYm)];
figure
hold on
fill(cXFill,cYFill,[.5 .5 .5]);
plot(cMX,cMY,'-k','LineWidth',2)
xlim([0 1])
ylim([0 1])
for ii = 1:12
    eXm(ii,:) = mean(cat(2,eachX{:,ii}),2);
    eYm(ii,:) = mean(cat(2,eachY{:,ii}),2);
    plot(eXm(ii,:),eYm(ii,:))
    hold on
end

%%
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\paper1data.mat')
allProb = [];
for ii = 1:2:24
    trainX = x(~ismember(1:24,ii:ii+1),:);
    trainY = y(~ismember(1:24,ii:ii+1),1);
    testX = x(ii:ii+1,:);
    testY = y(ii:ii+1,1);
    mdl = fitcsvm(trainX,trainY,'OptimizeHyperparameters','auto','HyperparameterOptimization',struct('kfold',3,'Verbose',0,'ShowPlots',0));
    scoreMdl{ii} = fitPosterior(mdl);    
    [prob] = predict(scoreMdl{ii},testX);
    allProb = [allProb;prob];
end
