%% Collate raw data: LSD vs. Saline
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\lsd\processed\'],{'lsd';'sal'},{'pow','coh'},...
    'trl','');
% For each animal, average that day's baseline and use to normalize post
allData = [data{1,1};data{1,2}];
mBase = cellfun(@(x) mean(x,1),allData(:,1),'UniformOutput',0);
for ii = 1:4
    dif{ii} = allData{ii,2}-repmat(mBase{ii},size(allData{ii,2},1),1);
end
postLSD = cat(1,dif{3:4});
postSal = cat(1,dif{1:2});
save('C:\Users\Pythia\Documents\GreenLab\data\lsd\normPost.mat',...
    'postLSD','postSal')
%% Collate data
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\lsd\processed\'],{'48hr';'lsd';'sal'},{'pow','coh'},...
    'trl','rel');
[waterData,waterSamp,waterFiles] = collateData(['C:\Users\Pythia\'...
    'Documents\GreenLab\data\paper3\waterProcessed\'],{'26';'29';'30';...
    '37';'38';'39';'80';'81';'82';'88';'90'},{'pow','coh'},'trl','rel');
waterData = cat(1,waterData{:});
waterFiles = cat(2,waterFiles{:})';
nameVect = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
%% Create animal variable for mixed effects models
group = cell(1,3);
for ii = 1:size(data,2)
    for jj = 1:size(data{ii},1)
        for k = 1:size(data{ii},2)
            group{ii}{jj,k} = repmat({files{ii}{jj}(4:5)},size(data{ii}...
                {jj,k},1),1);
        end
    end
end
n = cellfun(@(x) size(x,1),waterData,'UniformOutput',0);
waterGroup = cell(size(waterData));
for ii = 1:size(waterData,1)
    for jj = 1:size(waterData,2)
        waterGroup{ii,jj} = repmat({waterFiles{ii}(7:8)},n{ii,jj},1);
    end
end
%% Combine data across animals into two groups with base and post
% Get baseline saline and lsd
baseSal = cat(1,data{1,3}{:,1},waterData{[1:2,5:11,16:17],2});
baseSalG = cat(1,group{1,3}{:,1},waterGroup{[1:2,5:11,16:17],2});
baseLSD = cat(1,data{1,2}{:,1},waterData{[3:4,12:15,18:25],2});
baseLSDG = cat(1,group{1,2}{:,1},waterGroup{[3:4,12:15,18:25],2});
% Get 48hour saline and lsd
sal24 = cat(1,data{1,1}{[1,3:5,8],1});
sal24G = cat(1,group{1,1}{[1,3:5,8],1});
lsd24 = cat(1,data{1,1}{[2,6,7,9:11],1});
lsd24G = cat(1,group{1,1}{[2,6,7,9:11],1});
% Get post injection saline and lsd
postSal = cat(1,data{1,3}{:,2});
postSalG = cat(1,group{1,3}{:,2});
postLSD = cat(1,data{1,2}{:,2});
postLSDG = cat(1,group{1,2}{:,2});
%% Compare saline baselines to post injection
nTrain = 1300;
nTest = 325;
for ii = 1:100
    baseInd = randperm(size(baseSal,1),nTrain+nTest);
    postInd = randperm(size(postSal,1),nTrain+nTest);
    
    trainX = [baseSal(baseInd(1:nTrain),:);postSal(postInd(1:nTrain),:)];
    trainY = [zeros(nTrain,1);ones(nTrain,1)];
    testX = [baseSal(baseInd(nTrain+1:end),:);...
        postSal(postInd(nTrain+1:end),:)];
    testY = [zeros(nTest,1);ones(nTest,1)];
    mdl = fitglm(zscore(trainX),trainY,'distribution','binomial');
    beta(ii,:) = mdl.Coefficients.Estimate;
    prob = predict(mdl,zscore(testX));
    [fpr(:,ii),tpr(:,ii),~,a(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/100:1,'UseNearest',0);
    [prec(:,ii),recall(:,ii),~,a(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/100:1,'UseNearest',0,'XCrit','prec');
end
%% Compare saline baselines (+ post injection) to 24 hour
nTrain = 1300;
nTest = 325;
newBaseSal = [baseSal;postSal];
for ii = 1:100
    baseInd = randperm(size(newBaseSal,1),nTrain+nTest);
    postInd = randperm(size(sal24,1),nTrain+nTest);
    
    trainX = [newBaseSal(baseInd(1:nTrain),:);sal24(postInd(1:nTrain),:)];
    trainY = [zeros(nTrain,1);ones(nTrain,1)];
    testX = [newBaseSal(baseInd(nTrain+1:end),:);...
        sal24(postInd(nTrain+1:end),:)];
    testY = [zeros(nTest,1);ones(nTest,1)];
    mdl = fitglm(zscore(trainX),trainY,'distribution','binomial');
    beta24(ii,:) = mdl.Coefficients.Estimate;
    prob = predict(mdl,zscore(testX));
    [fpr24(:,ii),tpr24(:,ii),~,a24(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/100:1,'UseNearest',0);
    [prec24(:,ii),recall24(:,ii),~,a(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/100:1,'UseNearest',0,'XCrit','prec');
end
%% Compare base to post injection LSD
nTrain = 1300;
nTest = 325;
for ii = 1:100
    baseInd = randperm(size(baseLSD,1),nTrain+nTest);
    postInd = randperm(size(postLSD,1),nTrain+nTest);
    
    trainX = [baseLSD(baseInd(1:nTrain),:);postLSD(postInd(1:nTrain),:)];
    trainY = [zeros(nTrain,1);ones(nTrain,1)];
    testX = [baseLSD(baseInd(nTrain+1:end),:);...
        postLSD(postInd(nTrain+1:end),:)];
    testY = [zeros(nTest,1);ones(nTest,1)];
    mdl = fitglm(zscore(trainX),trainY,'distribution','binomial');
    lsdBeta(ii,:) = mdl.Coefficients.Estimate;
    prob = predict(mdl,zscore(testX));
    [lsdFPR(:,ii),lsdTPR(:,ii),~,lsdA(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/100:1,'UseNearest',0);
    [lsdPrec(:,ii),lsdRecall(:,ii),~,a(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/100:1,'UseNearest',0,'XCrit','prec');
end
%% Compare base to 24 hour LSD
nTrain = 1300;
nTest = 325;
for ii = 1:100
    baseInd = randperm(size(baseLSD,1),nTrain+nTest);
    postInd = randperm(size(lsd24,1),nTrain+nTest);
    
    trainX = [baseLSD(baseInd(1:nTrain),:);lsd24(postInd(1:nTrain),:)];
    trainY = [zeros(nTrain,1);ones(nTrain,1)];
    testX = [baseLSD(baseInd(nTrain+1:end),:);...
        lsd24(postInd(nTrain+1:end),:)];
    testY = [zeros(nTest,1);ones(nTest,1)];
    mdl = fitglm(zscore(trainX),trainY,'distribution','binomial');
    lsd24Beta(ii,:) = mdl.Coefficients.Estimate;
    prob = predict(mdl,zscore(testX));
    [lsd24FPR(:,ii),lsd24TPR(:,ii),~,lsd24A(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/100:1,'UseNearest',0);
    [lsd24Prec(:,ii),lsd24Recall(:,ii),~,a(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/100:1,'UseNearest',0,'XCrit','prec');
end
%% Compare LSD to saline acute
nTrain = 2600;
nTest = 650;
for ii = 1:100
    baseInd = randperm(size(postLSD,1),nTrain+nTest);
    postInd = randperm(size(postSal,1),nTrain+nTest);
    
    trainX = [postLSD(baseInd(1:nTrain),:);postSal(postInd(1:nTrain),:)];
    trainY = [zeros(nTrain,1);ones(nTrain,1)];
    testX = [postLSD(baseInd(nTrain+1:end),:);...
        postSal(postInd(nTrain+1:end),:)];
    testY = [zeros(nTest,1);ones(nTest,1)];
    mdl = fitglm(zscore(trainX),trainY,'distribution','binomial');
    lsdSalBeta(ii,:) = mdl.Coefficients.Estimate;
    prob = predict(mdl,zscore(testX));
    [lsdSalFPR(:,ii),lsdSalTPR(:,ii),~,lsdSalA(ii)] = perfcurve(testY,prob,...
        1,'TVals',0:1/100:1,'UseNearest',0);
    [lsdSalPrec(:,ii),lsdSalRecall(:,ii),~,a(ii)] = perfcurve(testY,prob,1,...
        'TVals',0:1/100:1,'UseNearest',0,'XCrit','prec');
end
%% Plot all models
figure
hold on
plot(mean(fpr,2),mean(tpr,2))
plot(mean(fpr24,2),mean(tpr24,2))
plot(mean(lsdFPR,2),mean(lsdTPR,2))
plot(mean(lsd24FPR,2),mean(lsd24TPR,2))
plot(mean(lsdSalFPR,2),mean(lsdSalTPR,2))
legend({'Sal Base vs. Sal Post','Sal Base + Post vs. Sal 24',...
    'LSD Base vs. LSD Post','LSD Base vs. LSD 24',...
    'LSD Post vs. Sal Post'},'location','se')
box off
xlabel('FPR'); ylabel('TPR')

figure
hold on
plot(mean(prec,2),mean(recall,2))
plot(mean(prec24,2),mean(recall24,2))
plot(mean(lsdPrec,2),mean(lsdRecall,2))
plot(mean(lsd24Prec,2),mean(lsd24Recall,2))
plot(mean(lsdSalPrec,2),mean(lsdSalRecall,2))
legend({'Sal Base vs. Sal Post','Sal Base + Post vs. Sal 24',...
    'LSD Base vs. LSD Post','LSD Base vs. LSD 24',...
    'LSD Post vs. Sal Post'},'location','sw')
box off
xlabel('Precision'); ylabel('Recall (TPR)')
%% Build LSD pre vs. acute model (Linear Mixed Effects)
nTrain = 2200;
nTest = 550;
% Combine data into table
t = array2table([zscore(baseLSD);zscore(postLSD)]);
t(:,end+1) = cell2table([baseLSDG;postLSDG]);
t(:,end+1) = array2table([ones(size(baseLSD,1),1);zeros(size(postLSD,1),...
    1)]);
parfor ii = 1:100
      baseLSDInd = randperm(size(baseLSD,1),nTrain+nTest);
      % Add number of baseLSD to get correct indices after concatenating
      % data
      postLSDInd = randperm(size(postLSD,1),nTrain+nTest)+size(baseLSD,1);
      % Separate trian and test data
      train = t([baseLSDInd(1:nTrain),postLSDInd(1:nTrain)],:);
      test = t([baseLSDInd(1:nTest),postLSDInd(1:nTest)],:);
      % Build linear mixed effects model with animal ID as random effect
      mdl = 
end
%% Single feature models
parfor ii = 1:100
    nTrain = 800;
    nTest = 200;
    % Build LSD model
    baseLSDInd = randperm(size(baseLSD,1),nTrain+nTest);
    lsd24Ind = randperm(size(lsd24,1),nTrain+nTest);
%     postLSDInd = randperm(size(postLSD,1),nTrain+nTest);
%     postSalInd = randperm(size(postSal,1),nTrain+nTest);
    trainX = [baseLSD(baseLSDInd(1:nTrain),:);...
        lsd24(lsd24Ind(1:nTrain),:)];
%     trainX = [baseLSD(baseLSDInd(1:nTrain),:);...
%         postLSD(postLSDInd(1:nTrain),:)];
%     trainX = [postSal(postSalInd(1:nTrain),:);...
%         postLSD(postLSDInd(1:nTrain),:)];

    trainY = [zeros(nTrain,1);ones(nTrain,1)];
    testX = [baseLSD(baseLSDInd(nTrain+1:end),:);...
        lsd24(lsd24Ind(nTrain+1:end),:)];
%     testX = [baseLSD(baseLSDInd(nTrain+1:end),:);...
%         postLSD(postLSDInd(nTrain+1:end),:)];
%     testX = [postSal(postSalInd(nTrain+1:end),:);...
%         postLSD(postLSDInd(nTrain+1:end),:)];
    testY = [zeros(nTest,1);ones(nTest,1)];
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        prob = predict(mdl,testX(:,jj));
        [~,~,~,a(ii,jj)] = perfcurve(testY,prob,1);
    end
end
%%
mA = mean(a,1);
for ii = 1:216
   c(ii) = conf(a(:,ii)',0.95); 
end
[sortA,sortI] = sort(mA,'descend');
figure
scatterErr(1:216,sortA,c(sortI),0)
%% Open Lasso results
lsdA = []; salA = [];
for ii = 1:100
%    load(['C:\Users\Pythia\Documents\GreenLab\data\lsd\lassoBaseAcute\',...
%        num2str(ii),'.mat']) 
%    load(['C:\Users\Pythia\Documents\GreenLab\data\lsd\lassoBase24\',...
%        num2str(ii),'.mat'])
   lsdA(ii) = aLSD;
   lsdX(ii,:) = xLSD;
   lsdY(ii,:) = yLSD;
   
   salA(ii) = aSal;
   salX(ii,:) = xSal;
   salY(ii,:) = ySal;
end
figure
plot(mean(salX,1),mean(salY,1))
hold on
plot(mean(lsdX,1),mean(lsdY,1))
legend({'Saline','LSD'},'location','se')
title('Pre vs. Acute')
xlabel('FPR')
ylabel('TPR')
box off
%% Saline vs. LSD: Acute and 24 Hour
for ii = 1:100
   load(['C:\Users\Pythia\Documents\GreenLab\data\lsd\lassoLSDvSal\',...
       num2str(ii),'.mat']) 
   acuteA(ii) = a;
   acuteX(ii,:) = x;
   acuteY(ii,:) = y;
   
   load(['C:\Users\Pythia\Documents\GreenLab\data\lsd\lasso24LSDvSal\',...
       num2str(ii),'.mat']) 
   a24(ii) = a;
   x24(ii,:) = x;
   y24(ii,:) = y;
end
figure
plot(mean(acuteX,1),mean(acuteY,1))
hold on
plot(mean(x24,1),mean(y24,1))
legend({'Acute','24 Hour'},'location','se')
title('LSD vs. Sal: Acute and 24 Hour')
xlabel('FPR')
ylabel('TPR')
box off
%% Go through LSD samples and index every x minutes
x = 10;
for ii = 1:size(samp{3},1)
    test = samp{3}{ii,2}(:,2)./400;
    c = 1;
    inds = [];
    while c < size(test,1)
        newC = logicFind(test(c)+(60*x),test,'>','first');
        if isempty(newC)
            newC = size(test,1);
        end
        inds = [inds;c,newC-1];
        c = newC;
    end
    for jj = 1:size(inds,1)
        post{ii}(jj,:) = mean(data{3}{ii,2}(inds(jj,1):inds(jj,2),:),1);
    end
end
%% Plot post feature y vs. baseline average of that feature
y = 134;
figure
hold on
for ii = 1:size(post,2)
    plot(1:size(inds,1),post{ii}(:,y))
end
shadedErrorBar([1 size(inds,1)],[mean(baseLSD(:,y)) mean(baseLSD(:,y))],...
    [std(baseLSD(:,y)) std(baseLSD(:,y))],':k',1)
% plot([1 size(inds,1)],[mean(baseLSD(:,y)) mean(baseLSD(:,y))],'k:')
set(gca,'xticklabel',x:x:x*size(inds,1))
xlabel('Minutes Post Injection')
title([nameVect{y}])
%% GLM
for jj = 1:2
    for ii = 1:100
        preInd = randperm(size(data{1,2}{jj,1},1),400);
        postInd = randperm(size(data{1,2}{jj,2},1),400);
        trainX = [data{1,2}{jj,1}(preInd(1:350),:);...
            data{1,2}{jj,2}(postInd(1:350),:)];
        trainY = [zeros(350,1);ones(350,1)];
        testX = [data{1,2}{jj,1}(preInd(351:400),:);...
            data{1,2}{jj,2}(postInd(351:400),:)];
        testY = [zeros(50,1);ones(50,1)];
        
        mdl = fitglm(trainX,trainY,'distribution','binomial');
        prob = predict(mdl,testX);
        [~,~,~,a(jj,ii)] = perfcurve(testY,prob,1);
    end
end
%% Lasso
for ii = 1:2
    for jj = 1:2
        for k = 1:100
            preInd = randperm(size(data{1,ii}{jj,1},1),400);
            postInd = randperm(size(data{1,ii}{jj,2},1),400);
            trainX = [data{1,ii}{jj,1}(preInd(1:350),:);...
                data{1,ii}{jj,2}(postInd(1:350),:)];
            trainY = [zeros(350,1);ones(350,1)];
            testX = [data{1,ii}{jj,1}(preInd(351:400),:);...
                data{1,ii}{jj,2}(postInd(351:400),:)];
            testY = [zeros(50,1);ones(50,1)];
            [cfg] = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se',[]);
            [allAlpha,allLambda,allBeta,cvFitsArray,accArray{ii,jj}(k),hist] = lassoNet(trainX,trainY,'binomial','auc',1,10,1,cfg);
        end
    end
end