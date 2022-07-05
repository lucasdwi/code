%% cohort 1
[c1Data,c1Samps,c1Files] = collateData(['F:\angelaAlcoholDual\'...
    'processedCohort1\'],{'Dual';'Con';'MIA';'AE'},{'pow','coh'},...
    'trl','rel');
% split by sex

% load cohort 2 pre data
load('F:/angelaAlcoholDual/cohort2.mat','c2PreData','c2PreFiles')
% add the four female dual rats
inds = 7:10;
c1c2Data = c1Data;
c1c2Files = c1Files;

%% Dual vs. con 
% LOO from each group
cmbs = [];
c = 1;
for ii = 1:8
   for jj = 1:15
      cmbs(c,:) = [ii*2-1,ii*2,jj*2-1,jj*2];
      c = c+1;
   end
end
% Go through all iterations of different left out pairs
for ii = 1:size(cmbs,1)
    dualTestInds = cmbs(ii,1:2);
    conTestInds = cmbs(ii,3:4);
    dualTrainInds = logicFind(1,~ismember(1:size(c1Data{1},1),...
        dualTestInds),'==');
    conTrainInds =  logicFind(1,~ismember(1:size(c1Data{2},1),...
        conTestInds),'==');
    trainX = [];
    for jj = dualTrainInds
        trainX = [trainX;c1Data{1}{jj}(randperm(size(c1Data{1}{jj},1),...
            1200),:)];
    end
    for jj = conTrainInds
        trainX = [trainX;c1Data{2}{jj}(randperm(size(c1Data{2}{jj},1),...
            1200),:)];
    end
    trainY = [ones(14*1200,1);zeros(28*1200,1)];
    
    testX = [];
    for jj = dualTestInds
        testX = [testX;c1Data{1}{jj}(randperm(size(c1Data{1}{jj},1),...
            1200),:)];
    end
    for jj = conTestInds
        testX = [testX;c1Data{2}{jj}(randperm(size(c1Data{2}{jj},1),...
            1200),:)];
    end
    testY = [ones(2400,1);zeros(2400,1)];
    % Single feature
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        pred = predict(mdl,testX(:,jj));
        [~,~,~,dcSingle(ii,jj)] = perfcurve(testY,pred,1);
        dcSign(ii,jj) = table2array(mdl.Coefficients(2,1));
    end
    % Real
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,c1DCA(ii)] = perfcurve(testY,pred,1);
    c1DCX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,4801));
    c1DCY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,4801));
    % Permuted
    mdl = fitglm(trainX(randperm(size(trainX,1),size(trainX,1)),:),...
        trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisXP,thisYP,~,c1DCAp(ii)] = perfcurve(testY,pred,1);
    c1DCXp(ii,:) = interp1(linspace(0,1,numel(thisXP)),thisXP,...
        linspace(0,1,4801));
    c1DCYp(ii,:) = interp1(linspace(0,1,numel(thisYP)),thisYP,...
        linspace(0,1,4801));
    disp(ii)
end
%% Con vs. MIA
cmbs = [];
c = 1;
for ii = 1:15
   for jj = 1:17
      cmbs(c,:) = [ii*2-1,ii*2,jj*2-1,jj*2];
      c = c+1;
   end
end
% Go through all iterations of different left out pairs
for ii = 1:size(cmbs,1)
    conTestInds = cmbs(ii,1:2);
    miaTestInds = cmbs(ii,3:4);
    conTrainInds = logicFind(1,~ismember(1:size(c1Data{2},1),...
        conTestInds),'==');
    miaTrainInds =  logicFind(1,~ismember(1:size(c1Data{3},1),...
        miaTestInds),'==');
    trainX = [];
    for jj = conTrainInds
        trainX = [trainX;c1Data{2}{jj}(randperm(size(c1Data{2}{jj},1),...
            1200),:)];
    end
    for jj = miaTrainInds
        trainX = [trainX;c1Data{3}{jj}(randperm(size(c1Data{3}{jj},1),...
            1200),:)];
    end
    trainY = [ones(28*1200,1);zeros(32*1200,1)];
    
    testX = [];
    for jj = conTestInds
        testX = [testX;c1Data{2}{jj}(randperm(size(c1Data{2}{jj},1),...
            1200),:)];
    end
    for jj = miaTestInds
        testX = [testX;c1Data{3}{jj}(randperm(size(c1Data{3}{jj},1),...
            1200),:)];
    end
    testY = [ones(2400,1);zeros(2400,1)];
    % Single feature
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        pred = predict(mdl,testX(:,jj));
        [~,~,~,cmSingle(ii,jj)] = perfcurve(testY,pred,1);
        cmSign(ii,jj) = table2array(mdl.Coefficients(2,1));
    end
    % Real
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,c1CMA(ii)] = perfcurve(testY,pred,1);
    c1CMX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,4801));
    c1CMY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,4801));
    % Permuted
    mdl = fitglm(trainX(randperm(size(trainX,1),size(trainX,1)),:),...
        trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisXP,thisYP,~,c1CMAp(ii)] = perfcurve(testY,pred,1);
    c1CMXp(ii,:) = interp1(linspace(0,1,numel(thisXP)),thisXP,...
        linspace(0,1,4801));
    c1CMYp(ii,:) = interp1(linspace(0,1,numel(thisYP)),thisYP,...
        linspace(0,1,4801));
    disp(ii)
end
%% Con vs. AE
cmbs = [];
c = 1;
for ii = 1:15
   for jj = 1:8
      cmbs(c,:) = [ii*2-1,ii*2,jj*2-1,jj*2];
      c = c+1;
   end
end
% Go through all iterations of different left out pairs
for ii = 1:size(cmbs,1)
    conTestInds = cmbs(ii,1:2);
    aeTestInds = cmbs(ii,3:4);
    conTrainInds = logicFind(1,~ismember(1:size(c1Data{2},1),...
        conTestInds),'==');
    aeTrainInds =  logicFind(1,~ismember(1:size(c1Data{4},1),...
        aeTestInds),'==');
    trainX = [];
    for jj = conTrainInds
        trainX = [trainX;c1Data{2}{jj}(randperm(size(c1Data{2}{jj},1),...
            1200),:)];
    end
    for jj = aeTrainInds
        trainX = [trainX;c1Data{4}{jj}(randperm(size(c1Data{4}{jj},1),...
            1200),:)];
    end
    trainY = [ones(28*1200,1);zeros(14*1200,1)];
    
    testX = [];
    for jj = conTestInds
        testX = [testX;c1Data{2}{jj}(randperm(size(c1Data{2}{jj},1),...
            1200),:)];
    end
    for jj = aeTestInds
        testX = [testX;c1Data{4}{jj}(randperm(size(c1Data{4}{jj},1),...
            1200),:)];
    end
    testY = [ones(2400,1);zeros(2400,1)];
    % Single feature
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        pred = predict(mdl,testX(:,jj));
        [~,~,~,caSingle(ii,jj)] = perfcurve(testY,pred,1);
        caSign(ii,jj) = table2array(mdl.Coefficients(2,1));
    end
    % Real
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,c1CAA(ii)] = perfcurve(testY,pred,1);
    c1CAX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,4801));
    c1CAY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,4801));
    % Permuted
    mdl = fitglm(trainX(randperm(size(trainX,1),size(trainX,1)),:),...
        trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisXP,thisYP,~,c1CAAp(ii)] = perfcurve(testY,pred,1);
    c1CAXp(ii,:) = interp1(linspace(0,1,numel(thisXP)),thisXP,...
        linspace(0,1,4801));
    c1CAYp(ii,:) = interp1(linspace(0,1,numel(thisYP)),thisYP,...
        linspace(0,1,4801));
    disp(ii)
end
save('c1mdls.mat','-append','caSingle','caSign','c1CAA','c1CAX','c1CAY',...
    'c1CAAp','c1CAXp','c1CAYp')
%%
for ii = 1:100
%     % leave 2 dual
%     these = randperm(8,2)';
%     dualTestInds = [these*2-1;these*2];
%     dualTrainInds = logicFind(1,~ismember(1:size(c1Data{1},1),...
%         dualTestInds),'==');
%     % and 1 from each other
%     this = randi(15,1);
%     conTestInds = [this*2-1;this*2];
%     conTrainInds = logicFind(1,~ismember(1:size(c1Data{2},1),...
%         conTestInds),'==');
%     this = randi(17,1);
%     miaTestInds = [this*2-1;this*2];
%     miaTrainInds = logicFind(1,~ismember(1:size(c1Data{3},1),...
%         miaTestInds),'==');
%     this = randi(8,1);
%     aeTestInds = [this*2-1;this*2];
%     aeTrainInds = logicFind(1,~ismember(1:size(c1Data{4},1),...
%         aeTestInds),'==');
%     % train
%     trainX = [];
%     for jj = dualTrainInds
%         trainX = [trainX;c1Data{1}{jj}(randperm(size(c1Data{1}{jj},1),...
%             1200),:)];
%     end
%     for jj = conTrainInds
%         trainX = [trainX;c1Data{2}{jj}(randperm(size(c1Data{2}{jj},1),...
%             1200),:)];
%     end
%     for jj = miaTrainInds
%         trainX = [trainX;c1Data{3}{jj}(randperm(size(c1Data{3}{jj},1),...
%             1200),:)];
%     end
%     for jj = aeTrainInds
%         trainX = [trainX;c1Data{4}{jj}(randperm(size(c1Data{4}{jj},1),...
%             1200),:)];
%     end
%     trainY = [ones(1200*12,1);zeros(1200*74,1)];
%     % test
%     testX = [];
%     for jj = dualTestInds'
%         testX = [testX;c1Data{1}{jj}(randperm(size(c1Data{1}{jj},1),...
%             1200),:)];
%     end
%     for jj = conTestInds'
%         testX = [testX;c1Data{2}{jj}(randperm(size(c1Data{2}{jj},1),...
%             1200),:)];
%     end
%     for jj = miaTestInds'
%         testX = [testX;c1Data{3}{jj}(randperm(size(c1Data{3}{jj},1),...
%             1200),:)];
%     end
%     for jj = aeTestInds'
%         testX = [testX;c1Data{4}{jj}(randperm(size(c1Data{4}{jj},1),...
%             1200),:)];
%     end
%     testY = [ones(1200*4,1);zeros(1200*6,1)];
%     % Single feature
%     for jj = 1:216
%         mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
%         pred = predict(mdl,testX(:,jj));
%         [~,~,~,dAllSingle(ii,jj)] = perfcurve(testY,pred,1);
%         dAllSign(ii,jj) = table2array(mdl.Coefficients(2,1));
%     end
%     % LOO
%     % Real
%     mdl = fitglm(trainX,trainY,'distribution','binomial');
%     pred = predict(mdl,testX);
%     [thisX,thisY,~,c1DAllA(ii)] = perfcurve(testY,pred,1);
%     c1DAllX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
%         linspace(0,1,4801));
%     c1DAllY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
%         linspace(0,1,4801));
%     % Permuted
%     mdl = fitglm(trainX(randperm(size(trainX,1),size(trainX,1)),:),...
%         trainY,'distribution','binomial');
%     pred = predict(mdl,testX);
%     [thisXP,thisYP,~,c1DAllAp(ii)] = perfcurve(testY,pred,1);
%     c1DAllXp(ii,:) = interp1(linspace(0,1,numel(thisXP)),thisXP,...
%         linspace(0,1,4801));
%     c1DAllYp(ii,:) = interp1(linspace(0,1,numel(thisYP)),thisYP,...
%         linspace(0,1,4801));
    % 80/20
    xAll = [];
    for jj = 1:4
        for k = 1:size(c1Data{jj},1)
            xAll = [xAll;...
                c1Data{jj}{k}(randperm(size(c1Data{jj}{k},1),1200),:)];
        end
    end
    % Hard coded y
    yAll = [ones(16*1200,1);zeros(80*1200,1)];
    [trainXall,trainYall,testXall,testYall] = trainTest(xAll,yAll,0.2);
    mdl = fitglm(trainXall,trainYall,'distribution','binomial');
    pred = predict(mdl,testXall);
    [thisX,thisY,~,c1DAllA80(ii)] = perfcurve(testYall,pred,1);
    c1DAllX80(ii,:) = interp1(linspace(0,1,numel(thisXP)),thisXP,...
        linspace(0,1,23041));
    c1DAllY80(ii,:) = interp1(linspace(0,1,numel(thisYP)),thisYP,...
        linspace(0,1,23041));
    disp(ii)
end
% save('c1mdls.mat','-append','dAllSingle','dAllSign','c1DAllA','c1DAllX',...
%     'c1DAllY','c1DAllAp','c1DAllXp','c1DAllYp')
%% Dual vs. con
figure
hold on
plot(mean(c1DCX,1),mean(c1DCY,1),'-k')
plot(mean(c1DCXp,1),mean(c1DCYp,1),'--k')
legend({['LOO: ',num2str(round(mean(c1DCA),2)),'\pm',...
    num2str(round(conf(c1DCA,0.95),2))],['permuted: ',...
    num2str(round(mean(c1DCAp),2)),'\pm',num2str(round(conf(c1DCAp,0.95)...
    ,2))]},'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('cohort 1; dual vs. con; GLM; LOO')
box off
this = sign(mean(dcSign,1));
[sortDC,dcInd] = sort(mean(dcSingle,1),'descend');
sortDC = sortDC'.*this(dcInd)';
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'})';
sortFeat = feat(dcInd);
%% MIA vs. con
figure
hold on
plot(mean(c1CMX,1),mean(c1CMY,1),'-k')
plot(mean(c1CMXp,1),mean(c1CMYp,1),'--k')
legend({['LOO: ',num2str(round(mean(c1CMA),2)),'\pm',...
    num2str(round(conf(c1CMA,0.95),2))],['permuted: ',...
    num2str(round(mean(c1CMAp),2)),'\pm',num2str(round(conf(c1CMAp,0.95)...
    ,2))]},'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('cohort 1; mia vs. con; GLM; LOO')
box off
this = sign(mean(cmSign,1));
[sortCM,cmInd] = sort(mean(cmSingle,1),'descend');
sortCM = sortCM'.*this(cmInd)';
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'})';
sortFeat = feat(cmInd);
%% AE vs. con
figure
hold on
plot(mean(c1CAX,1),mean(c1CAY,1),'-k')
plot(mean(c1CAXp,1),mean(c1CAYp,1),'--k')
legend({['LOO: ',num2str(round(mean(c1CAA),2)),'\pm',...
    num2str(round(conf(c1CAA,0.95),2))],['permuted: ',...
    num2str(round(mean(c1CAAp),2)),'\pm',num2str(round(conf(c1CAAp,0.95)...
    ,2))]},'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('cohort 1; ae vs. con; GLM; LOO')
box off
this = sign(mean(caSign,1));
[sortCA,caInd] = sort(mean(caSingle,1),'descend');
sortCA = sortCA'.*this(caInd)';
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'})';
sortFeat = feat(caInd);
%% Dual vs. all
figure
hold on
plot(mean(c1DAllX,1),mean(c1DAllY,1),'-k')
plot(mean(c1DAllXp,1),mean(c1DAllYp,1),'--k')
plot(mean(c1DAllX80,1),mean(c1DAllY80,1),':k')
legend({['LOO: ',num2str(round(mean(c1DAllA),2)),'\pm',...
    num2str(round(conf(c1DAllA,0.95),2))],['permuted: ',...
    num2str(round(mean(c1DAllAp),2)),'\pm',...
    num2str(round(conf(c1DAllAp,0.95),2))],['80/20: ',...
    num2str(round(mean(c1DAllA80),2)),'\pm',...
    num2str(round(conf(c1DAllA80,0.95),2))]},'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('cohort 1; dual vs. all; GLM; LOO')
box off
this = sign(mean(dAllSign,1));
[sortDAll,dAllInd] = sort(mean(dAllSingle,1),'descend');
sortDAll = sortDAll'.*this(dAllInd)';
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'})';
sortFeat = feat(dAllInd);
%% Load pre and post
[c2PreData,c2PreSamps,c2PreFiles] = collateData(...
    'F:\angelaAlcoholDual\processedCohort2\Pre\',{'Dual0_';'Dual05_';...
    'Dual10_';'Con0_';'Con05_';'Con10_'},{'pow','coh'},'trl','rel');
[c2PostData,c2PostSamps,c2PostFiles] = collateData(...
    'F:\angelaAlcoholDual\processedCohort2\Post\',{'Dual0_';'Dual05_';...
    'Dual10_';'Con0_';'Con05_';'Con10_'},{'pow','coh'},'trl','rel');
% split post data into first 30 minutes and post 45 minutes
alc30 = cell(6,12);
post45 = cell(6,12);
for ii = 1:numel(c2PostSamps)
    for jj = 1:numel(c2PostSamps{ii})
        ind30 = logicFind(1,c2PostSamps{ii}{jj}>=30*60,'==','first');
        ind45 = logicFind(1,c2PostSamps{ii}{jj}>=45*60,'==','first');
        alc30{ii,jj} = c2PostData{ii}{jj}(1:ind30,:);
        post45{ii,jj} = c2PostData{ii}{jj}(ind45:end,:);
    end
end
% find out number of samples per animal
thisID = cell(6,12);
thisGroup = cell(6,12);
for ii = 1:numel(c2PreFiles)
    for jj = 1:numel(c2PreFiles{ii})
        parts = strsplit(c2PreFiles{ii}{jj},'_');
        thisID{ii,jj} = parts{2};
        thisGroup{ii,jj} = parts{1};       
    end
end
uID = unique([unique(thisID(1,1:10)),unique(thisID(2,1:11)),...
    unique(thisID(3,1:12)),unique(thisID(4,1:9)),unique(thisID(5,1:8)),...
    unique(thisID(6,1:10))]);
last = [10,11,12,9,8,10];
for ii = 1:numel(uID)
    for jj = 1:6
        inds = logicFind(1,cellfun(@(x) ~isempty(x),...
            strfind(thisID(jj,1:last(jj)),uID{ii})),'==');
        onSamps{ii,jj} = size(cat(1,alc30{jj,inds}),1);
    end
end
%% Get collated feats for base, 0-30 min, and 31-60 min
load('F:/angelaAlcoholDual/cohort2.mat')
% Set up indices of interest
% inds = [1,6,24,42,48]; % old features
% Look at all features
inds = 1:216;
% Preallocate three cell arrays
[dataBase,data30,data60] = deal(cell(1,6));
% Get baseline data
for ii = 1:6
    % Combine all data per group
    dataBase{ii} = cat(1,c2PreData{1,ii}{:});
    % Pull out indices of interest
    dataBase{ii} = dataBase{ii}(:,inds);
end
data = cell(6,4);
for ii = 1:6
    for jj = 1:numel(c2PostSamps{ii})
        % 0-30 min and 31-60 min
        % Find index for 30 minutes
        samp30 = nearest_idx3(30*60,c2PostSamps{ii}{jj});
        % Find index for 60 minutes
        samp60 = nearest_idx3(60*60,c2PostSamps{ii}{jj});
        % Grab data from each recording and concatenate with others
        data30{ii} = [data30{ii};c2PostData{ii}{jj}(1:samp30,inds)];
        % Use samp30+1 to start with first index after 30 minutes
        data60{ii} = [data60{ii};c2PostData{ii}{jj}(samp30+1:samp60,inds)];
        % Grab every 15 minutes
        for h = 1:4
            start = nearest_idx3(1+(h-1)*15*60,c2PostSamps{ii}{jj});
            stop = nearest_idx3(h*15*60,c2PostSamps{ii}{jj});
            data{ii,h} = [data{ii,h};c2PostData{ii}{jj}(start:stop,:)];
        end
    end
end
% save('F:/angelaAlcoholDual/featuresOfInterest.mat','dataBase','data30',...
%     'data60')
%%
for ii = 1:216
    % ttest between combined dual 0.5 & 1.0 mg/kg and con 0 mg/kg at
    % baseline
%     [~,p(ii,1)] = ttest2([dataBase{3}(:,ii);dataBase{2}(:,ii)],...
%         dataBase{4}(:,ii));
    % ttest between combined dual 1.0 mg/kg and con 0 mg/kg at baseline
    [~,p(ii,1)] = ttest2(dataBase{3}(:,ii),dataBase{4}(:,ii));
    % ttest between combined dual 1.0 mg/kg and con 0 mg/kg, first 30
    % minutes of alcohol 
    [~,p(ii,2)] = ttest2([data30{3}(:,ii);data30{2}(:,ii)],...
        data30{4}(:,ii));
end
% Bonferroni correction
pAdj = p.*216;
% Find features that are significantly different at baseline and not in
% first 30 minutes
pSig = pAdj(:,1)<=0.05 & pAdj(:,2)>=0.05;

feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'})';
col = {'-r','--r',':r','k'};
% Plot significant features
for ii = 1%logicFind(1,pSig,'==')
    figure
    hold on
    for jj = 1:4
        % Plot base, 30, and 60 minutes
        errorbar([1,2,3],[mean(dataBase{jj}(:,ii)),mean(data30{jj}(:,ii)),...
            mean(data60{jj}(:,ii))],[std(dataBase{jj}(:,ii)),...
            std(data30{jj}(:,ii)),std(data60{jj}(:,ii))],col{jj})
        % Plot base, and every 15 minutes
%         errorbar([1:5],[mean(dataBase{jj}(:,ii)),mean(data{jj,1}(:,ii)),...
%             mean(data{jj,2}(:,ii)),mean(data{jj,3}(:,ii)),...
%             mean(data{jj,4}(:,ii))],[std(dataBase{jj}(:,ii)),...
%             std(data{jj,1}(:,ii)),std(data{jj,2}(:,ii)),...
%             std(data{jj,3}(:,ii)),std(data{jj,4}(:,ii))],col{jj})
    end
    legend({'D0','D05','D10','C0'})
    title(feat(ii))
end
%% dual0 on vs. off - GLM
% combine offs
for ii = 1:numel(c2PreData{1})
    allOff{ii} = [c2PreData{1}{ii};post45{1,ii}];
end
% combine data from same animals
thisOff = {allOff{1},cat(1,allOff{2:3}),allOff{4:10}};
thisAlc30 = {alc30{1},cat(1,alc30{2:3}),alc30{4:10}};
for k = 1:20
    % LOO
    for ii = 1:numel(thisOff)
        testX = [thisOff{ii}(randperm(size(thisOff{ii},1),1200),:);...
            thisAlc30{1,ii}(randperm(size(thisAlc30{1,ii},1),1200),:)];
        testY = [zeros(1200,1);ones(1200,1)];
        trainInds = logicFind(1,~ismember(1:numel(thisOff),ii),'==');
        trainX = [];
        trainY = [];
        for jj = trainInds
            trainX = [trainX;thisOff{jj}(randperm(size(thisOff{jj},1),...
                1200),:);thisAlc30{1,jj}(randperm(size(thisAlc30{1,jj},...
                1),1200),:)];
            trainY = [trainY;zeros(1200,1);ones(1200,1)];
        end
        % GLM
        mdlLOO_0 = fitglm(trainX,trainY,'distribution','binomial');
        pred = predict(mdlLOO_0,testX);
        [glmLOOX_0{k,ii},glmLOOY_0{k,ii},~,glmLOOA_0(k,ii)] = perfcurve(...
            testY,pred,1);
    end
    % 80/20
    allX = [];
    allY = [];
    for ii = 1:numel(thisOff)
        allX = [allX;thisOff{ii}(randperm(size(thisOff{ii},1),...
            1200),:);thisAlc30{1,ii}(randperm(size(thisAlc30{1,ii},1),...
            1200),:)];
        allY = [allY;zeros(1200,1);ones(1200,1)];
    end
    [trainX,trainY,testX,testY] = trainTest(allX,allY,0.2);
    % GLM
    mdl80_0 = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl80_0,testX);
    [glm80X_0{k},glm80Y_0{k},~,glm80A_0(k)] = perfcurve(testY,pred,1);
    disp(k)
end
save('c2DualOnOff_0_glm.mat','glmLOOX_0','glmLOOY_0','glmLOOA_0',...
    'glm80X_0','glm80Y_0','glm80A_0')
% dual05 on vs. off - GLM
allOff = cell(1,numel(c2PreData{2}));
% combine offs
for ii = 1:numel(c2PreData{2})
    allOff{ii} = [c2PreData{2}{ii};post45{2,ii}];
end
% combine data from same animals
thisOff = {allOff{1},cat(1,allOff{2:3}),allOff{4},cat(1,allOff{5:6}),...
    allOff{7:8},cat(1,allOff{9:10}),allOff{11}};
thisAlc30 = {alc30{1},cat(1,alc30{2:3}),alc30{4},cat(1,alc30{5:6}),...
    alc30{7:8},cat(1,alc30{9:10}),alc30{11}};
for k = 1:20
    % LOO
    for ii = 1:numel(thisOff)
        testX = [thisOff{ii}(randperm(size(thisOff{ii},1),1200),:);...
            thisAlc30{1,ii}(randperm(size(thisAlc30{1,ii},1),1200),:)];
        testY = [zeros(1200,1);ones(1200,1)];
        trainInds = logicFind(1,~ismember(1:numel(thisOff),ii),'==');
        trainX = [];
        trainY = [];
        for jj = trainInds
            trainX = [trainX;thisOff{jj}(randperm(size(thisOff{jj},1),...
                1200),:);thisAlc30{1,jj}(randperm(size(thisAlc30{1,jj},...
                1),1200),:)];
            trainY = [trainY;zeros(1200,1);ones(1200,1)];
        end
        % GLM
        mdlLOO_05 = fitglm(trainX,trainY,'distribution','binomial');
        pred = predict(mdlLOO_05,testX);
        [glmLOOX_05{k,ii},glmLOOY_05{k,ii},~,glmLOOA_05(k,ii)] = ...
            perfcurve(testY,pred,1);
        for jj = 1:216
           mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
           pred = predict(mdl,testX(:,jj));
           [~,~,~,sA05(k,ii,jj)] = perfcurve(testY,pred,1);
           sA05Sign(k,ii,jj) = sign(table2array(mdl.Coefficients(2,1))); 
        end
    end
    % 80/20
    allX = [];
    allY = [];
    for ii = 1:numel(thisOff)
        allX = [allX;thisOff{ii}(randperm(size(thisOff{ii},1),...
            1200),:);thisAlc30{1,ii}(randperm(size(thisAlc30{1,ii},1),...
            1200),:)];
        allY = [allY;zeros(1200,1);ones(1200,1)];
    end
    [trainX,trainY,testX,testY] = trainTest(allX,allY,0.2);
    % GLM
    mdl80_05 = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl80_05,testX);
    [glm80X_05{k},glm80Y_05{k},~,glm80A_05(k)] = perfcurve(testY,pred,1);
    disp(k)
end
save('c2DualOnOff_05_glm.mat','glmLOOX_05','glmLOOY_05','glmLOOA_05',...
    'glm80X_05','glm80Y_05','glm80A_05','sA05','sA05Sign')
% dual10 on vs. off - GLM
allOff = cell(1,numel(c2PreData{3}));
% combine offs
for ii = 1:numel(c2PreData{3})
    allOff{ii} = [c2PreData{3}{ii};post45{3,ii}];
end
% combine data from same animals
thisOff = {cat(1,allOff{1:2}),allOff{3:4},cat(1,allOff{5:6}),allOff{7},...
    cat(1,allOff{8:9}),allOff{10},cat(1,allOff{11:12})};
thisAlc30 = {cat(1,alc30{1:2}),alc30{3:4},cat(1,alc30{5:6}),alc30{7},...
    cat(1,alc30{8:9}),alc30{10},cat(1,alc30{11:12})};
for k = 1:20
    % LOO
    for ii = 1:numel(thisOff)
        testX = [thisOff{ii}(randperm(size(thisOff{ii},1),1200),:);...
            thisAlc30{1,ii}(randperm(size(thisAlc30{1,ii},1),1200),:)];
        testY = [zeros(1200,1);ones(1200,1)];
        trainInds = logicFind(1,~ismember(1:numel(thisOff),ii),'==');
        trainX = [];
        trainY = [];
        for jj = trainInds
            trainX = [trainX;thisOff{jj}(randperm(size(thisOff{jj},1),...
                1200),:);thisAlc30{1,jj}(randperm(size(thisAlc30{1,jj},...
                1),1200),:)];
            trainY = [trainY;zeros(1200,1);ones(1200,1)];
        end
        % GLM
        mdlLOO_10 = fitglm(trainX,trainY,'distribution','binomial');
        pred = predict(mdlLOO_10,testX);
        [glmLOOX_10{k,ii},glmLOOY_10{k,ii},~,glmLOOA_10(k,ii)] = ...
            perfcurve(testY,pred,1);
        for jj = 1:216
           mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
           pred = predict(mdl,testX(:,jj));
           [~,~,~,sA10(k,ii,jj)] = perfcurve(testY,pred,1);
           sA10Sign(k,ii,jj) = sign(table2array(mdl.Coefficients(2,1))); 
        end
    end
    % 80/20
    allX = [];
    allY = [];
    for ii = 1:numel(thisOff)
        allX = [allX;thisOff{ii}(randperm(size(thisOff{ii},1),...
            1200),:);thisAlc30{1,ii}(randperm(size(thisAlc30{1,ii},1),...
            1200),:)];
        allY = [allY;zeros(1200,1);ones(1200,1)];
    end
    [trainX,trainY,testX,testY] = trainTest(allX,allY,0.2);
    % GLM
    mdl80_10 = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl80_10,testX);
    [glm80X_10{k},glm80Y_10{k},~,glm80A_10(k)] = perfcurve(testY,pred,1);
    disp(k)
end
save('c2DualOnOff_10_glm.mat','glmLOOX_10','glmLOOY_10','glmLOOA_10',...
    'glm80X_10','glm80Y_10','glm80A_10','sA10','sA10Sign')
%%
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'})';

load('c2DualOnOff_05_glm.mat','sA05','sA05Sign')
a05 = squeeze(mean(mean(sA05,1),2)).*sign(squeeze(mean(mean(sA05Sign))));
[~,sortA05Ind] = sort(abs(a05),'descend');
sortA05 = a05(sortA05Ind);
a05Feat = feat(sortA05Ind);

load('c2DualOnOff_10_glm.mat','sA10','sA10Sign')
a10 = squeeze(mean(sA10,[1,2])).*sign(squeeze(mean(sA10Sign,[1,2])));
[~,sortA10Ind] = sort(abs(a10),'descend');
sortA10 = a10(sortA10Ind);
a10Feat = feat(sortA10Ind);
%% Plot glm dual on vs. off models
% dual0; LOO
for ii = 1:size(glmLOOX_0,1)
    for jj = 1:size(glmLOOX_0,2)
        x = glmLOOX_0{ii,jj};
        y = glmLOOY_0{ii,jj};
        thisX_0(ii,jj,:) = interp1(linspace(0,1,numel(x)),x,...
            linspace(0,1,2401));
        thisY_0(ii,jj,:) = interp1(linspace(0,1,numel(y)),y,...
            linspace(0,1,2401));
    end
end
% dual0; 80/20
for ii = 1:size(glm80X_0,2)
    x = glm80X_0{ii};
    y = glm80Y_0{ii};
    thisX80_0(ii,:) = interp1(linspace(0,1,numel(x)),x,linspace(0,1,4321));
    thisY80_0(ii,:) = interp1(linspace(0,1,numel(y)),y,linspace(0,1,4321));
end
% dual05; LOO
for ii = 1:size(glmLOOX_05,1)
    for jj = 1:size(glmLOOX_05,2)
        x = glmLOOX_05{ii,jj};
        y = glmLOOY_05{ii,jj};
        thisX_05(ii,jj,:) = interp1(linspace(0,1,numel(x)),x,...
            linspace(0,1,2401));
        thisY_05(ii,jj,:) = interp1(linspace(0,1,numel(y)),y,...
            linspace(0,1,2401));
    end
end
% dual05; 80/20
for ii = 1:size(glm80X_05,2)
    x = glm80X_05{ii};
    y = glm80Y_05{ii};
    thisX80_05(ii,:) = interp1(linspace(0,1,numel(x)),x,...
        linspace(0,1,4321));
    thisY80_05(ii,:) = interp1(linspace(0,1,numel(y)),y,...
        linspace(0,1,4321));
end
% dual10; LOO
for ii = 1:size(glmLOOX_10,1)
    for jj = 1:size(glmLOOX_10,2)
        x = glmLOOX_10{ii,jj};
        y = glmLOOY_10{ii,jj};
        thisX_10(ii,jj,:) = interp1(linspace(0,1,numel(x)),x,...
            linspace(0,1,2401));
        thisY_10(ii,jj,:) = interp1(linspace(0,1,numel(y)),y,...
            linspace(0,1,2401));
    end
end
% dual10; 80/20
for ii = 1:size(glm80X_10,2)
    x = glm80X_10{ii};
    y = glm80Y_10{ii};
    thisX80_10(ii,:) = interp1(linspace(0,1,numel(x)),x,...
        linspace(0,1,4321));
    thisY80_10(ii,:) = interp1(linspace(0,1,numel(y)),y,...
        linspace(0,1,4321));
end
figure
hold on
plot(squeeze(mean(mean(thisX_0))),squeeze(mean(mean(thisY_0))))
plot(squeeze(mean(mean(thisX_05))),squeeze(mean(mean(thisY_05))))
plot(squeeze(mean(mean(thisX_10))),squeeze(mean(mean(thisY_10))))
title('dual on vs. off; GLM; LOO')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['0: ',num2str(round(mean(mean(glmLOOA_0)),2)),'\pm',...
    num2str(round(conf(reshape(glmLOOA_0,1,numel(glmLOOA_0)),0.95),2))],...
    ['5: ',num2str(round(mean(mean(glmLOOA_05)),2)),'\pm',...
    num2str(round(conf(reshape(glmLOOA_05,1,numel(glmLOOA_05)),0.95),2))],...
    ['10: ',num2str(round(mean(mean(glmLOOA_10)),2)),'\pm',...
    num2str(round(conf(reshape(glmLOOA_10,1,numel(glmLOOA_10)),0.95),2))]})
xlabel('FPR'); ylabel('TPR')
box off

figure
hold on
plot(squeeze(mean(thisX80_0)),squeeze(mean(thisY80_0)))
plot(squeeze(mean(thisX80_05)),squeeze(mean(thisY80_05)))
plot(squeeze(mean(thisX80_10)),squeeze(mean(thisY80_10)))
title('dual on vs. off; GLM; 80/20')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['0: ',num2str(round(mean(glm80A_0),2)),'\pm',...
    num2str(round(conf(glm80A_0,0.95),3))],...
    ['5: ',num2str(round(mean(glm80A_05),2)),'\pm',...
    num2str(round(conf(glm80A_05,0.95),3))],...
    ['10: ',num2str(round(mean(glm80A_10),2)),'\pm',...
    num2str(round(conf(glm80A_10,0.95),3))]})
xlabel('FPR'); ylabel('TPR')
box off
%% Apply cohort2 dual and control to cohort1 models
dualPreX = cat(1,c2PreData{1}{:},c2PreData{2}{:},c2PreData{3}{:});
conPreX = cat(1,c2PreData{4}{:},c2PreData{5}{:},c2PreData{6}{:});
dualPreY = ones(size(dualPreX,1),1);
conPreY = zeros(size(conPreX,1),1);

allX = [dualPreX;conPreX];
allY = [dualPreY;conPreY];
for ii = 1:20
    load(['F:/angelaAlcoholDual/c1DualCon/c1DualCon',num2str(ii),'.mat']...
        ,'acc')
    pred = cvglmnetPredict(acc{1}.mdl{1},zscore(allX));
    [x(ii,:),y(ii,:),~,a(ii)] = perfcurve(allY,pred,1);
end
%% Combine cohorts
% Reload cohort 1 data, with all animals
[c1Data,c1Samps,c1Files] = collateData(['F:\angelaAlcoholDual\',...
    'processedCohort1\'],{'Dual';'Con';'AE';'MIA'},{'pow','coh'},'trl','rel');
% Load cohort 2, pre data
load('F:\angelaAlcoholDual\cohort2.mat','c2PreData','c2PreFiles',...
    'c2PreSamps')
%% Go through each animal and pull out the same total number of samples
% Cohort 1
c1ID = cell(4,31);
for ii = 1:numel(c1Data)
    for jj = 1:numel(c1Data{ii})
        parts = strsplit(c1Files{ii}{jj},'_');
        c1ID{ii,jj} = parts{1};
        c1Group{ii,jj} = parts{2};
    end    
end
for ii = 1:size(c1ID,1)
   last = logicFind(0,cellfun(@(x) isempty(x),c1ID(1,:)),'==','last');
   c1uID{ii} = unique(c1ID(ii,1:last)); 
end
% Combine data from the same animal in cohort 1
c = 1;
for jj = 1:size(c1uID,2)
    for ii = 1:numel(c1uID{jj})
        c1DataNew{c,2} = cat(1,c1Data{jj}{contains(c1Files{jj}(:,1),...
            c1uID{jj}(ii))});
        c1DataNew{c,1} = c1uID{jj}(ii);
        c1DataNew{c,3} = c1Group{jj,ii};
        c = c+1;
    end
end
% Cohort 2
c2ID = cell(6,12);
thisGroup = cell(6,12);
c = 1;
c2parts = cell(1,4);
for ii = 1:numel(c2PreFiles)
    for jj = 1:numel(c2PreFiles{ii})
        parts = strsplit(c2PreFiles{ii}{jj},'_');
        c2ID{ii,jj} = parts{2};
        c2parts(c,:) = [parts(1:3),c2PreData{ii}{jj}];
        c = c+1;   
    end
end
% Grab first recording (no experience with etoh injections)
uID = unique(c2parts(:,2));
c2FirstData = cell(numel(uID),4);
for ii = 1:numel(uID)
    these = c2parts(contains(c2parts(:,2),uID{ii}),1:4);
    [~,ind] = sort(these(:,3));
    c2FirstData(ii,:) = these(ind(1),1:4);
end
% HARD CODE SEX
% Grab all male and female duals
maleDual = [c1DataNew(contains(c1DataNew(:,4),'M')&...
    contains(c1DataNew(:,3),'Dual'),1:2);...
    c2FirstData(contains(c2FirstData(:,1),'Dual')&...
    contains(c2FirstData(:,5),'M'),[2,4]),];
femaleDual = [c1DataNew(contains(c1DataNew(:,4),'F')&...
    contains(c1DataNew(:,3),'Dual'),1:2);...
    c2FirstData(contains(c2FirstData(:,1),'Dual')&...
    contains(c2FirstData(:,5),'F'),[2,4])];
% Grab all male and female cons
maleCon =  [c1DataNew(contains(c1DataNew(:,4),'M')&...
    contains(c1DataNew(:,3),'Con'),1:2);...
    c2FirstData(contains(c2FirstData(:,1),'Con')&...
    contains(c2FirstData(:,5),'M'),[2,4]),];
femaleCon = [c1DataNew(contains(c1DataNew(:,4),'F')&...
    contains(c1DataNew(:,3),'Con'),1:2);...
    c2FirstData(contains(c2FirstData(:,1),'Con')&...
    contains(c2FirstData(:,5),'F'),[2,4])];
% Grab all male and female AE
maleAE =  [c1DataNew(contains(c1DataNew(:,4),'M')&...
    contains(c1DataNew(:,3),'AE'),1:2);...
    c2FirstData(contains(c2FirstData(:,1),'AE')&...
    contains(c2FirstData(:,3),'M'),[2,4])];
femaleAE = [c1DataNew(contains(c1DataNew(:,4),'F')&...
    contains(c1DataNew(:,3),'AE'),1:2);...
    c2FirstData(contains(c2FirstData(:,1),'AE')&...
    contains(c2FirstData(:,3),'F'),[2,4])];
% Grab all male and female MIA
maleMIA =  [c1DataNew(contains(c1DataNew(:,4),'M')&...
    contains(c1DataNew(:,3),'MIA'),1:2);...
    c2FirstData(contains(c2FirstData(:,1),'MIA')&...
    contains(c2FirstData(:,3),'M'),[2,4]),];
femaleMIA = [c1DataNew(contains(c1DataNew(:,4),'F')&...
    contains(c1DataNew(:,3),'MIA'),1:2);...
    c2FirstData(contains(c2FirstData(:,1),'MIA')&...
    contains(c2FirstData(:,3),'F'),[2,4])];
% save('c1c2SexData.mat','c1DataNew','c2FirstData','maleDual',...
% 'femaleDual','maleCon','femaleCon','maleAE','femaleAE','maleMIA','femaleMIA')
c2uID{1} = unique([unique(c2ID(1,1:10)),unique(c2ID(2,1:11)),...
    unique(c2ID(3,1:12))]);
c2uID{2} = unique([unique(c2ID(4,1:9)),unique(c2ID(5,1:8)),...
    unique(c2ID(6,1:10))]);
% Grab data from duals
for jj = 1:numel(c1uID{1})
    inds = contains(c1ID(1,1:16),c1uID{1}(jj));
    dualCat{1,jj} = cat(1,c1Data{1}{inds});
end
last = [10,11,12];
for jj = 1:numel(c2uID{1})
    thisCat = [];
    for ii = 1:3
        inds = contains(c2ID(ii,1:last(ii)),c2uID{1}(jj));
        thisCat = [thisCat;cat(1,c2PreData{ii}{inds})];
    end
    dualCat{2,jj} = thisCat;
end
% Grab data from cons
for jj = 1:numel(c1uID{1})
    inds = contains(c1ID(2,1:31),c1uID{2}(jj));
    conCat{1,jj} = cat(1,c1Data{2}{inds});
end
last = [9,8,10];
for jj = 1:numel(c2uID{2})
    thisCat = [];
    for ii = 4:6
        inds = contains(c2ID(ii,1:last(ii-3)),c2uID{2}(jj));
        thisCat = [thisCat;cat(1,c2PreData{ii}{inds})];
    end
    conCat{2,jj} = thisCat;
end
% Combine cohorts
allDual = [dualCat(1,1:8),dualCat(2,:)];
allCon = [conCat(1,:),conCat(2,:)];

%% Take 1200 samples from each animal
% For LOO use one control and one dual, from either cohort
c = 1;
for ii = 1:16
    for jj = 1:17
        cmbs(c,:) = [ii,jj];
        c = c+1;
    end
end
for ii = 1:size(cmbs,1)
    conInd = cmbs(ii,1);
    dualInd = cmbs(ii,2);
    conTrainInd = ~ismember(1:16,conInd);
    dualTrainInd = ~ismember(1:17,dualInd);
    for jj = 1:16
        subCon{jj} = allCon{jj}(randperm(size(allCon{jj},1),1200),:);
    end
    for jj = 1:17
        subDual{jj} = allDual{jj}(randperm(size(allDual{jj},1),1200),:);
    end
    % Real
    mdlX = cat(1,subCon{conTrainInd},subDual{dualTrainInd});
    mdlTestX = cat(1,subCon{conInd},subDual{dualInd});
    mdl = fitglm(mdlX,[zeros(1200*15,1);ones(1200*16,1)],...
        'distribution','binomial');
    pred = predict(mdl,mdlTestX);
    [x{ii},y{ii},~,a(ii)] = perfcurve([zeros(1200,1);ones(1200,1)],...
        pred,1);
    % Single feature
    for k = 1:216
        mdl = fitglm(mdlX(:,k),[zeros(1200*15,1);ones(1200*16,1)],...
            'distribution','binomial');
        pred = predict(mdl,mdlTestX(:,k));
        [~,~,~,sA(ii,k)] = perfcurve([zeros(1200,1);ones(1200,1)],...
            pred,1);
        sASign(ii,k) = sign(table2array(mdl.Coefficients(2,1)));
    end
    % Permuted
    catX = cat(1,subCon{conTrainInd},subDual{dualTrainInd});
    mdl = fitglm(catX(randperm(size(catX,1),size(catX,1)),:),...
        [zeros(1200*15,1);ones(1200*16,1)],'distribution','binomial');
    pred = predict(mdl,cat(1,subCon{conInd},subDual{dualInd}));
    [thisX,thisY,~,aP(ii)] = perfcurve([zeros(1200,1);ones(1200,1)],...
        pred,1);
    xP(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
       linspace(0,1,2401));
    yP(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
       linspace(0,1,2401));
    disp(ii)
end
% save('F:\angelaAlcoholDual\c1c2_dualCon.mat','x','y','a','xP','yP',...
%     'aP','sA')
%%
for ii = 1:numel(x)
   newX(ii,:) = interp1(linspace(0,1,numel(x{ii})),x{ii},...
       linspace(0,1,2401));
   newY(ii,:) = interp1(linspace(0,1,numel(y{ii})),y{ii},...
       linspace(0,1,2401));
end
figure
hold on
plot(mean(newX),mean(newY),'-k')
plot(mean(xP),mean(yP),'--k')
legend({['LOO: ',num2str(round(mean(a),2)),'\pm',...
    num2str(round(conf(a,0.95),2))],['Permuted: ',...
    num2str(round(mean(aP),2)),'\pm',num2str(round(conf(aP,0.95),2))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
title('c1c2; dual vs. con; GLM; LOO')
box off
%% Build sex specific models with 1200 samples each rat, LOO and animal detector
% load data
load('F:\angelaAlcoholDual\c1c2SexData.mat')
% remove animals without histology
maleDual(7,:) = [];
femaleDual(4,:) = [];
femaleCon([2,9],:) = [];
% male dual vs. male con 11v8 (10v8 without missing histology)
c = 1;
clear maleCmb
for ii = 1:size(maleDual,1)
   for jj = 1:size(maleCon,1)
        maleCmb(c,:) = [ii,jj];
        c = c+1;
   end
end
for ii = 1:size(maleCmb,1)
    disp(ii)
    subSampCon = cellfun(@(x) x(randperm(size(x,1),1200),:),maleCon(:,2),'uniformoutput',0);
    subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),maleDual(:,2),'uniformoutput',0);
    thisDual = subSampDual{maleCmb(ii,1)};
    thisCon = subSampCon{maleCmb(ii,2)};
    testX = [thisDual;thisCon];
    testY = [ones(1200,1);zeros(1200,1)];
    otherDualInd = logicFind(1,~ismember(1:size(maleDual,1),maleCmb(ii,1)),'==');
    theseDual = cat(1,subSampDual{otherDualInd});
    otherConInd = logicFind(1,~ismember(1:size(maleCon,1),maleCmb(ii,2)),'==');
    theseCon = cat(1,subSampCon{otherConInd});
    trainX = [theseDual;theseCon];
    trainY = [ones(size(theseDual,1),1);zeros(size(theseCon,1),1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    mC(ii,:) = table2array(mdl.Coefficients(:,1));
    pred = predict(mdl,testX);
    [thisX,thisY,~,mA(ii)] = perfcurve(testY,pred,1);
    mX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    mY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
    for jj = 1:size(trainX,2)
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        mCS(ii,jj) = table2array(mdl.Coefficients(2,1));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,mCSA(ii,jj)] = perfcurve(testY,pred,1);
    end
end
% male dual vs. male con sub-sampled to 6 vs 6 (5v5 train)
for ii = 1:size(maleCmb,1)
    disp(ii)
    subSampCon = cellfun(@(x) x(randperm(size(x,1),1200),:),maleCon(:,2),'uniformoutput',0);
    subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),maleDual(:,2),'uniformoutput',0);
    thisDual = subSampDual{maleCmb(ii,1)};
    thisCon = subSampCon{maleCmb(ii,2)};
    testX = [thisDual;thisCon];
    testY = [ones(1200,1);zeros(1200,1)];
    otherDualInd = logicFind(1,~ismember(1:size(maleDual,1),maleCmb(ii,1)),'==');
    otherDualInd = otherDualInd(randperm(numel(otherDualInd),numel(otherDualInd)));
    theseDual = cat(1,subSampDual{otherDualInd(1:5)});
    otherConInd = logicFind(1,~ismember(1:size(maleCon,1),maleCmb(ii,2)),'==');
    otherConInd = otherConInd(randperm(numel(otherConInd),numel(otherConInd)));
    theseCon = cat(1,subSampCon{otherConInd(1:5)});
    trainX = [theseDual;theseCon];
    trainY = [ones(size(theseDual,1),1);zeros(size(theseCon,1),1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    mSubC(ii,:) = table2array(mdl.Coefficients(:,1));
    pred = predict(mdl,testX);
    [thisX,thisY,~,mSubA(ii)] = perfcurve(testY,pred,1);
    mSubX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    mSubY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
    for jj = 1:size(trainX,2)
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        mSubCS(ii,jj) = table2array(mdl.Coefficients(2,1));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,mSubCSA(ii,jj)] = perfcurve(testY,pred,1);
    end
end
% female dual vs. female con 6v9 (5v8 train); 5v7 without missing histology
c = 1;
clear femaleCmb
for ii = 1:size(femaleDual,1)
   for jj = 1:size(femaleCon,1)
        femaleCmb(c,:) = [ii,jj];
        c = c+1;
   end
end
for ii = 1:size(femaleCmb,1)
    disp(ii)
    subSampCon = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleCon(:,2),'uniformoutput',0);
    subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleDual(:,2),'uniformoutput',0);
    thisDual = subSampDual{femaleCmb(ii,1)};
    thisCon = subSampCon{femaleCmb(ii,2)};
    testX = [thisDual;thisCon];
    testY = [ones(1200,1);zeros(1200,1)];
    otherDualInd = logicFind(1,~ismember(1:size(femaleDual,1),femaleCmb(ii,1)),'==');
    theseDual = cat(1,subSampDual{otherDualInd});
    otherConInd = logicFind(1,~ismember(1:size(femaleCon,1),femaleCmb(ii,2)),'==');
    theseCon = cat(1,subSampCon{otherConInd});
    trainX = [theseDual;theseCon];
    trainY = [ones(size(theseDual,1),1);zeros(size(theseCon,1),1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    fC(ii,:) = table2array(mdl.Coefficients(:,1));
    pred = predict(mdl,testX);
    [thisX,thisY,~,fA(ii)] = perfcurve(testY,pred,1);
    fX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    fY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
    for jj = 1:size(trainX,2)
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        fCS(ii,jj) = table2array(mdl.Coefficients(2,1));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,fCSA(ii,jj)] = perfcurve(testY,pred,1);
    end
end
% female dual vs. female con sub-sampled to 6 vs 6 (5v5 train)
for ii = 1:size(femaleCmb,1)
    disp(ii)
    subSampCon = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleCon(:,2),'uniformoutput',0);
    subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleDual(:,2),'uniformoutput',0);
    thisDual = subSampDual{femaleCmb(ii,1)};
    thisCon = subSampCon{femaleCmb(ii,2)};
    testX = [thisDual;thisCon];
    testY = [ones(1200,1);zeros(1200,1)];
    otherDualInd = logicFind(1,~ismember(1:size(femaleDual,1),femaleCmb(ii,1)),'==');
    otherDualInd = otherDualInd(randperm(numel(otherDualInd),numel(otherDualInd)));
    theseDual = cat(1,subSampDual{otherDualInd(1:numel(otherDualInd))});
    otherConInd = logicFind(1,~ismember(1:size(femaleCon,1),femaleCmb(ii,2)),'==');
    otherConInd = otherConInd(randperm(numel(otherConInd),numel(otherConInd)));
    theseCon = cat(1,subSampCon{otherConInd(1:numel(otherConInd))});
    trainX = [theseDual;theseCon];
    trainY = [ones(size(theseDual,1),1);zeros(size(theseCon,1),1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    fSubC(ii,:) = table2array(mdl.Coefficients(:,1));
    pred = predict(mdl,testX);
    [thisX,thisY,~,fSubA(ii)] = perfcurve(testY,pred,1);
    fSubX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    fSubY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
    for jj = 1:size(trainX,2)
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        fSubCS(ii,jj) = table2array(mdl.Coefficients(2,1));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,fSubCSA(ii,jj)] = perfcurve(testY,pred,1);
    end
end
%%
% animal detector
allData = [maleDual;maleCon];
% 10 v 8 male
for ii = 1:100
    disp(ii)
    thisPerm = allData(randperm(size(allData,1),size(allData,1)),:);
    group1 = cat(1,thisPerm(1:10,2));
    group2 = cat(1,thisPerm(11:18,2));
    subG1 = cellfun(@(x) x(randperm(size(x,1),1200),:),group1,...
        'uniformoutput',0);
    subG2 = cellfun(@(x) x(randperm(size(x,1),1200),:),group2,...
        'uniformoutput',0);
    testX = [subG1{1};subG2{1}];
    testY = [ones(1200,1);zeros(1200,1)];
    trainX = [cat(1,subG1{2:end});cat(1,subG2{2:end})];
    trainY = [ones(size(cat(1,subG1{2:end}),1),1);...
        zeros(size(cat(1,subG2{2:end}),1),1)];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,aMR(ii)] = perfcurve(testY,pred,1);
    mrX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    mrY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
end
% 6 v 6 male sub-sampled
for ii = 1:100
    disp(ii)
    thisPerm = allData(randperm(size(allData,1),size(allData,1)),:);
    group1 = cat(1,thisPerm(1:6,2));
    group2 = cat(1,thisPerm(7:12,2));
    subG1 = cellfun(@(x) x(randperm(size(x,1),1200),:),group1,...
        'uniformoutput',0);
    subG2 = cellfun(@(x) x(randperm(size(x,1),1200),:),group2,...
        'uniformoutput',0);
    testX = [subG1{1};subG2{1}];
    testY = [ones(1200,1);zeros(1200,1)];
    trainX = [cat(1,subG1{2:end});cat(1,subG2{2:end})];
    trainY = [ones(size(cat(1,subG1{2:end}),1),1);...
        zeros(size(cat(1,subG2{2:end}),1),1)];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,aMSubR(ii)] = perfcurve(testY,pred,1);
    mrSubX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    mrSubY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
end
% 6 vs. 9 female (5 vs 7 without missing histology)
allData = [femaleDual;femaleCon];
for ii = 1:100
    disp(ii)
    thisPerm = allData(randperm(size(allData,1),size(allData,1)),:);
    group1 = cat(1,thisPerm(1:5,2));
    group2 = cat(1,thisPerm(6:12,2));
    subG1 = cellfun(@(x) x(randperm(size(x,1),1200),:),group1,...
        'uniformoutput',0);
    subG2 = cellfun(@(x) x(randperm(size(x,1),1200),:),group2,...
        'uniformoutput',0);
    testX = [subG1{1};subG2{1}];
    testY = [ones(1200,1);zeros(1200,1)];
    trainX = [cat(1,subG1{2:end});cat(1,subG2{2:end})];
    trainY = [ones(size(cat(1,subG1{2:end}),1),1);...
        zeros(size(cat(1,subG2{2:end}),1),1)];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,frA(ii)] = perfcurve(testY,pred,1);
    frX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    frY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
end
% 6 vs. 6 female sub-sampled
for ii = 1:100
    disp(ii)
    thisPerm = allData(randperm(size(allData,1),size(allData,1)),:);
    group1 = cat(1,thisPerm(1:6,2));
    group2 = cat(1,thisPerm(7:12,2));
    subG1 = cellfun(@(x) x(randperm(size(x,1),1200),:),group1,...
        'uniformoutput',0);
    subG2 = cellfun(@(x) x(randperm(size(x,1),1200),:),group2,...
        'uniformoutput',0);
    testX = [subG1{1};subG2{1}];
    testY = [ones(1200,1);zeros(1200,1)];
    trainX = [cat(1,subG1{2:end});cat(1,subG2{2:end})];
    trainY = [ones(size(cat(1,subG1{2:end}),1),1);...
        zeros(size(cat(1,subG2{2:end}),1),1)];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,frSubA(ii)] = perfcurve(testY,pred,1);
    frSubX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    frSubY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
end
%% Dual vs. con figure
figure
hold on
plot(mean(mX),mean(mY),'k')
plot(mean(mSubX),mean(mSubY),':k')
plot(mean(mrX),mean(mrY),'--k')
legend({['male: ',num2str(round(mean(mA),2)),'\pm',...
    num2str(round(conf(mA,0.95),2))],['male sub: ',...
    num2str(round(mean(mSubA),2)),'\pm',...
    num2str(round(conf(mSubA,0.95),2))],['animal detector: ',...
    num2str(round(mean(aMR),2)),'\pm',num2str(round(conf(aMR,0.95),2))]}...
    ,'location','se')
box off
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
title('male dual vs. con; 10 vs. 8; GLM; LOO')

figure
hold on
plot(mean(fX),mean(fY),'k')
plot(mean(fSubX),mean(fSubY),':k')
plot(mean(frX),mean(frY),'--k')
legend({['female: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],['female sub: ',...
    num2str(round(mean(fSubA),2)),'\pm',...
    num2str(round(conf(fSubA,0.95),2))],['animal detector: ',...
    num2str(round(mean(frA),2)),'\pm',num2str(round(conf(frA,0.95),2))]}...
    ,'location','se')
box off
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
title('female dual vs. con; 5 vs. 7; GLM; LOO')
%% Dual vs. MIA
% male dual vs. male mia
c = 1;
clear maleCmb
for ii = 1:size(maleDual,1)
   for jj = 1:size(maleMIA,1)
        maleCmb(c,:) = [ii,jj];
        c = c+1;
   end
end
% for ii = 1:size(maleCmb,1)
%     disp(ii)
%     subSampMIA = cellfun(@(x) x(randperm(size(x,1),1200),:),maleMIA(:,2),'uniformoutput',0);
%     subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),maleDual(:,2),'uniformoutput',0);
%     thisDual = subSampDual{maleCmb(ii,1)};
%     thisMIA = subSampMIA{maleCmb(ii,2)};
%     testX = [thisDual;thisMIA];
%     testY = [ones(1200,1);zeros(1200,1)];
%     otherDualInd = logicFind(1,~ismember(1:size(maleDual,1),maleCmb(ii,1)),'==');
%     theseDual = cat(1,subSampDual{otherDualInd});
%     otherMIAInd = logicFind(1,~ismember(1:size(maleMIA,1),maleCmb(ii,2)),'==');
%     theseMIA = cat(1,subSampMIA{otherMIAInd});
%     trainX = [theseDual;theseMIA];
%     trainY = [ones(size(theseDual,1),1);zeros(size(theseMIA,1),1)];
%     
%     mdl = fitglm(trainX,trainY,'distribution','binomial');
%     mCmia(ii,:) = table2array(mdl.Coefficients(:,2));
%     pred = predict(mdl,testX);
%     [thisX,thisY,~,mA(ii)] = perfcurve(testY,pred,1);
%     mXmia(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
%         linspace(0,1,2401));
%     mY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
%         linspace(0,1,2401));
%     for jj = 1:size(trainX,2)
%         mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
%         mCSmia(ii,jj) = table2array(mdl.Coefficients(2,2));
%         pred = predict(mdl,testX(:,jj));
%         [~,~,~,mCSAmia(ii,jj)] = perfcurve(testY,pred,1);
%     end
% end
% male dual vs. male mia sub-sampled to 3 vs 3 (2v2 train)
for ii = 1:size(maleCmb,1)
    disp(ii)
    subSampMIA = cellfun(@(x) x(randperm(size(x,1),1200),:),maleMIA(:,2),'uniformoutput',0);
    subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),maleDual(:,2),'uniformoutput',0);
    thisDual = subSampDual{maleCmb(ii,1)};
    thisMIA = subSampMIA{maleCmb(ii,2)};
    testX = [thisDual;thisMIA];
    testY = [ones(1200,1);zeros(1200,1)];
    otherDualInd = logicFind(1,~ismember(1:size(maleDual,1),maleCmb(ii,1)),'==');
    otherDualInd = otherDualInd(randperm(numel(otherDualInd),numel(otherDualInd)));
    theseDual = cat(1,subSampDual{otherDualInd(1:2)});
    otherMIAInd = logicFind(1,~ismember(1:size(maleMIA,1),maleCmb(ii,2)),'==');
    otherMIAInd = otherMIAInd(randperm(numel(otherMIAInd),numel(otherMIAInd)));
    theseMIA = cat(1,subSampMIA{otherMIAInd(1:2)});
    trainX = [theseDual;theseMIA];
    trainY = [ones(size(theseDual,1),1);zeros(size(theseMIA,1),1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    mSubMIAC(ii,:) = table2array(mdl.Coefficients(:,1));
    pred = predict(mdl,testX);
    [thisX,thisY,~,mSubAmia(ii)] = perfcurve(testY,pred,1);
    mSubXmia(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    mSubYmia(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
    for jj = 1:size(trainX,2)
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        mSubMIACS(ii,jj) = table2array(mdl.Coefficients(2,1));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,mSubMIACSA(ii,jj)] = perfcurve(testY,pred,1);
    end
end
% female dual vs. female mia
c = 1;
clear femaleCmb
for ii = 1:size(femaleDual,1)
   for jj = 1:size(femaleMIA,1)
        femaleCmb(c,:) = [ii,jj];
        c = c+1;
   end
end
% female dual vs. female mia 6 vs. 6 (5v5 train)
for ii = 1:size(femaleCmb,1)
    disp(ii)
    subSampMIA = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleMIA(:,2),'uniformoutput',0);
    subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleDual(:,2),'uniformoutput',0);
    thisDual = subSampDual{femaleCmb(ii,1)};
    thisMIA = subSampMIA{femaleCmb(ii,2)};
    testX = [thisDual;thisMIA];
    testY = [ones(1200,1);zeros(1200,1)];
    otherDualInd = logicFind(1,~ismember(1:size(femaleDual,1),femaleCmb(ii,1)),'==');
    otherDualInd = otherDualInd(randperm(numel(otherDualInd)));
    theseDual = cat(1,subSampDual{otherDualInd(1:numel(otherDualInd))});
    otherMIAInd = logicFind(1,~ismember(1:size(femaleMIA,1),femaleCmb(ii,2)),'==');
    otherMIAInd = otherMIAInd(randperm(numel(otherMIAInd)));
    theseMIA = cat(1,subSampMIA{otherMIAInd(1:numel(otherMIAInd))});
    trainX = [theseDual;theseMIA];
    trainY = [ones(size(theseDual,1),1);zeros(size(theseMIA,1),1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    fMIAC(ii,:) = table2array(mdl.Coefficients(:,1));
    pred = predict(mdl,testX);
    [thisX,thisY,~,fAmia(ii)] = perfcurve(testY,pred,1);
    fXmia(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    fYmia(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
    for jj = 1:size(trainX,2)
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        fMIACS(ii,jj) = table2array(mdl.Coefficients(2,1));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,fMIACSA(ii,jj)] = perfcurve(testY,pred,1);
    end
end
%% female dual vs. female mia sub-sampled to 3 vs 3 (2v2 train)
for ii = 1:size(femaleCmb,1)
    disp(ii)
    subSampMIA = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleMIA(:,2),'uniformoutput',0);
    subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleDual(:,2),'uniformoutput',0);
    thisDual = subSampDual{femaleCmb(ii,1)};
    thisMIA = subSampMIA{femaleCmb(ii,2)};
    testX = [thisDual;thisMIA];
    testY = [ones(1200,1);zeros(1200,1)];
    otherDualInd = logicFind(1,~ismember(1:size(femaleDual,1),femaleCmb(ii,1)),'==');
    otherDualInd = otherDualInd(randperm(numel(otherDualInd),numel(otherDualInd)));
    theseDual = cat(1,subSampDual{otherDualInd(1:2)});
    otherMIAInd = logicFind(1,~ismember(1:size(femaleMIA,1),femaleCmb(ii,2)),'==');
    otherMIAInd = otherMIAInd(randperm(numel(otherMIAInd),numel(otherMIAInd)));
    theseMIA = cat(1,subSampMIA{otherMIAInd(1:2)});
    trainX = [theseDual;theseMIA];
    trainY = [ones(size(theseDual,1),1);zeros(size(theseMIA,1),1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,fSubAmia(ii)] = perfcurve(testY,pred,1);
    fSubXmia(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    fSubYmia(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
end
% 3v3 male animal detector
allData = [maleDual;maleMIA];
for ii = 1:100
    disp(ii)
    thisPerm = allData(randperm(size(allData,1),size(allData,1)),:);
    group1 = cat(1,thisPerm(1:3,2));
    group2 = cat(1,thisPerm(4:6,2));
    subG1 = cellfun(@(x) x(randperm(size(x,1),1200),:),group1,...
        'uniformoutput',0);
    subG2 = cellfun(@(x) x(randperm(size(x,1),1200),:),group2,...
        'uniformoutput',0);
    testX = [subG1{1};subG2{1}];
    testY = [ones(1200,1);zeros(1200,1)];
    trainX = [cat(1,subG1{2:end});cat(1,subG2{2:end})];
    trainY = [ones(size(cat(1,subG1{2:end}),1),1);...
        zeros(size(cat(1,subG2{2:end}),1),1)];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,aMSubRmia(ii)] = perfcurve(testY,pred,1);
    mrSubXmia(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    mrSubYmia(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
end
% 3v3 female animal detector
allData = [femaleDual;femaleMIA];
for ii = 1:100
    disp(ii)
    thisPerm = allData(randperm(size(allData,1),size(allData,1)),:);
    group1 = cat(1,thisPerm(1:3,2));
    group2 = cat(1,thisPerm(4:6,2));
    subG1 = cellfun(@(x) x(randperm(size(x,1),1200),:),group1,...
        'uniformoutput',0);
    subG2 = cellfun(@(x) x(randperm(size(x,1),1200),:),group2,...
        'uniformoutput',0);
    testX = [subG1{1};subG2{1}];
    testY = [ones(1200,1);zeros(1200,1)];
    trainX = [cat(1,subG1{2:end});cat(1,subG2{2:end})];
    trainY = [ones(size(cat(1,subG1{2:end}),1),1);...
        zeros(size(cat(1,subG2{2:end}),1),1)];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,aFSubRmia(ii)] = perfcurve(testY,pred,1);
    frSubXmia(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    frSubYmia(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
end
%% Dual vs. mia figure
figure
hold on
plot(mean(mSubXmia),mean(mSubYmia),'k')
plot(mean(mrSubXmia),mean(mrSubYmia),'--k')
legend({['male: ',num2str(round(mean(mSubAmia),2)),'\pm',...
    num2str(round(conf(mSubAmia,0.95),2))],['animal detector: ',...
    num2str(round(mean(aMSubRmia),2)),'\pm',...
    num2str(round(conf(aMSubRmia,0.95),2))]},'location','se')
box off
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
title('male dual vs. mia; 3 vs. 3; GLM; LOO')

figure
hold on
plot(mean(fX),mean(fY),'k')
plot(mean(fSubXmia),mean(fSubYmia),':k')
plot(mean(frSubXmia),mean(frSubYmia),'--k')
legend({['female: ',num2str(round(mean(fA),2)),'\pm',...
    num2str(round(conf(fA,0.95),2))],['female sub-sample: ',...
    num2str(round(mean(fSubAmia),2)),'\pm',...
    num2str(round(conf(fSubAmia,0.95),2))],['animal detector: ',...
    num2str(round(mean(aFSubRmia),2)),'\pm',...
    num2str(round(conf(aFSubRmia,0.95),2))]},'location','se')
box off
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
title('female dual vs. mia; 3 vs. 3; GLM; LOO')
%% Dual vs. AE
% male dual vs. male ae
c = 1;
clear maleCmb
for ii = 1:size(maleDual,1)
   for jj = 1:size(maleAE,1)
        maleCmb(c,:) = [ii,jj];
        c = c+1;
   end
end
% for ii = 1:size(maleCmb,1)
%     disp(ii)
%     subSampAE = cellfun(@(x) x(randperm(size(x,1),1200),:),maleAE(:,2),'uniformoutput',0);
%     subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),maleDual(:,2),'uniformoutput',0);
%     thisDual = subSampDual{maleCmb(ii,1)};
%     thisAE = subSampAE{maleCmb(ii,2)};
%     testX = [thisDual;thisAE];
%     testY = [ones(1200,1);zeros(1200,1)];
%     otherDualInd = logicFind(1,~ismember(1:size(maleDual,1),maleCmb(ii,1)),'==');
%     theseDual = cat(1,subSampDual{otherDualInd});
%     otherAEInd = logicFind(1,~ismember(1:size(maleAE,1),maleCmb(ii,2)),'==');
%     theseAE = cat(1,subSampAE{otherAEInd});
%     trainX = [theseDual;theseAE];
%     trainY = [ones(size(theseDual,1),1);zeros(size(theseAE,1),1)];
%     
%     mdl = fitglm(trainX,trainY,'distribution','binomial');
%     mCAE(ii,:) = table2array(mdl.Coefficients(:,2));
%     pred = predict(mdl,testX);
%     [thisX,thisY,~,mA(ii)] = perfcurve(testY,pred,1);
%     mXAE(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
%         linspace(0,1,2401));
%     mY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
%         linspace(0,1,2401));
%     for jj = 1:size(trainX,2)
%         mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
%         mCSAE(ii,jj) = table2array(mdl.Coefficients(2,2));
%         pred = predict(mdl,testX(:,jj));
%         [~,~,~,mCSAAE(ii,jj)] = perfcurve(testY,pred,1);
%     end
% end
% male dual vs. male AE sub-sampled to 4 vs 4 (3v3 train)
for ii = 1:size(maleCmb,1)
    disp(ii)
    subSampAE = cellfun(@(x) x(randperm(size(x,1),1200),:),maleAE(:,2),'uniformoutput',0);
    subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),maleDual(:,2),'uniformoutput',0);
    thisDual = subSampDual{maleCmb(ii,1)};
    thisAE = subSampAE{maleCmb(ii,2)};
    testX = [thisDual;thisAE];
    testY = [ones(1200,1);zeros(1200,1)];
    otherDualInd = logicFind(1,~ismember(1:size(maleDual,1),maleCmb(ii,1)),'==');
    otherDualInd = otherDualInd(randperm(numel(otherDualInd),numel(otherDualInd)));
    theseDual = cat(1,subSampDual{otherDualInd(1:3)});
    otherAEInd = logicFind(1,~ismember(1:size(maleAE,1),maleCmb(ii,2)),'==');
    otherAEInd = otherAEInd(randperm(numel(otherAEInd),numel(otherAEInd)));
    theseAE = cat(1,subSampAE{otherAEInd(1:3)});
    trainX = [theseDual;theseAE];
    trainY = [ones(size(theseDual,1),1);zeros(size(theseAE,1),1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    mSubAEC(ii,:) = table2array(mdl.Coefficients(:,1));
    pred = predict(mdl,testX);
    [thisX,thisY,~,mSubAAE(ii)] = perfcurve(testY,pred,1);
    mSubXAE(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    mSubYAE(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
    for jj = 1:size(trainX,2)
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        mSubAECS(ii,jj) = table2array(mdl.Coefficients(2,1));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,mSubAECSA(ii,jj)] = perfcurve(testY,pred,1);
    end
end
% female dual vs. female AE
c = 1;
clear femaleCmb
for ii = 1:size(femaleDual,1)
   for jj = 1:size(femaleAE,1)
        femaleCmb(c,:) = [ii,jj];
        c = c+1;
   end
end
% for ii = 1:size(femaleCmb,1)
%     disp(ii)
%     subSampAE = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleAE(:,2),'uniformoutput',0);
%     subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleDual(:,2),'uniformoutput',0);
%     thisDual = subSampDual{femaleCmb(ii,1)};
%     thisAE = subSampAE{femaleCmb(ii,2)};
%     testX = [thisDual;thisAE];
%     testY = [ones(1200,1);zeros(1200,1)];
%     otherDualInd = logicFind(1,~ismember(1:size(femaleDual,1),femaleCmb(ii,1)),'==');
%     theseDual = cat(1,subSampDual{otherDualInd});
%     otherAEInd = logicFind(1,~ismember(1:size(femaleAE,1),femaleCmb(ii,2)),'==');
%     theseAE = cat(1,subSampAE{otherAEInd});
%     trainX = [theseDual;theseAE];
%     trainY = [ones(size(theseDual,1),1);zeros(size(theseAE,1),1)];
%     
%     mdl = fitglm(trainX,trainY,'distribution','binomial');
%     fC(ii,:) = table2array(mdl.Coefficients(:,2));
%     pred = predict(mdl,testX);
%     [thisX,thisY,~,fA(ii)] = perfcurve(testY,pred,1);
%     fX(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
%         linspace(0,1,2401));
%     fY(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
%         linspace(0,1,2401));
%     for jj = 1:size(trainX,2)
%         mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
%         fCS(ii,jj) = table2array(mdl.Coefficients(2,2));
%         pred = predict(mdl,testX(:,jj));
%         [~,~,~,fCSA(ii,jj)] = perfcurve(testY,pred,1);
%     end
% end
% female dual vs. female AE sub-sampled to 4 vs 4 (3v3 train)
for ii = 1:size(femaleCmb,1)
    disp(ii)
    subSampAE = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleAE(:,2),'uniformoutput',0);
    subSampDual = cellfun(@(x) x(randperm(size(x,1),1200),:),femaleDual(:,2),'uniformoutput',0);
    thisDual = subSampDual{femaleCmb(ii,1)};
    thisAE = subSampAE{femaleCmb(ii,2)};
    testX = [thisDual;thisAE];
    testY = [ones(1200,1);zeros(1200,1)];
    otherDualInd = logicFind(1,~ismember(1:size(femaleDual,1),femaleCmb(ii,1)),'==');
    otherDualInd = otherDualInd(randperm(numel(otherDualInd),numel(otherDualInd)));
    theseDual = cat(1,subSampDual{otherDualInd(1:3)});
    otherAEInd = logicFind(1,~ismember(1:size(femaleAE,1),femaleCmb(ii,2)),'==');
    otherAEInd = otherAEInd(randperm(numel(otherAEInd),numel(otherAEInd)));
    theseAE = cat(1,subSampAE{otherAEInd(1:3)});
    trainX = [theseDual;theseAE];
    trainY = [ones(size(theseDual,1),1);zeros(size(theseAE,1),1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    mSubAEC(ii,:) = table2array(mdl.Coefficients(:,1));
    pred = predict(mdl,testX);
    [thisX,thisY,~,fSubAAE(ii)] = perfcurve(testY,pred,1);
    fSubXAE(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    fSubYAE(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
    for jj = 1:size(trainX,2)
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        fSubAECS(ii,jj) = table2array(mdl.Coefficients(2,1));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,fSubAECSA(ii,jj)] = perfcurve(testY,pred,1);
    end
end
%% 4v4 male animal detector
allData = [maleDual;maleAE];
for ii = 1:100
    disp(ii)
    thisPerm = allData(randperm(size(allData,1),size(allData,1)),:);
    group1 = cat(1,thisPerm(1:4,2));
    group2 = cat(1,thisPerm(5:8,2));
    subG1 = cellfun(@(x) x(randperm(size(x,1),1200),:),group1,...
        'uniformoutput',0);
    subG2 = cellfun(@(x) x(randperm(size(x,1),1200),:),group2,...
        'uniformoutput',0);
    testX = [subG1{1};subG2{1}];
    testY = [ones(1200,1);zeros(1200,1)];
    trainX = [cat(1,subG1{2:end});cat(1,subG2{2:end})];
    trainY = [ones(size(cat(1,subG1{2:end}),1),1);...
        zeros(size(cat(1,subG2{2:end}),1),1)];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,aMSubRAE(ii)] = perfcurve(testY,pred,1);
    mrSubXAE(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    mrSubYAE(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
end
% 4v4 female animal detector
allData = [femaleDual;femaleAE];
for ii = 1:100
    disp(ii)
    thisPerm = allData(randperm(size(allData,1),size(allData,1)),:);
    group1 = cat(1,thisPerm(1:4,2));
    group2 = cat(1,thisPerm(5:8,2));
    subG1 = cellfun(@(x) x(randperm(size(x,1),1200),:),group1,...
        'uniformoutput',0);
    subG2 = cellfun(@(x) x(randperm(size(x,1),1200),:),group2,...
        'uniformoutput',0);
    testX = [subG1{1};subG2{1}];
    testY = [ones(1200,1);zeros(1200,1)];
    trainX = [cat(1,subG1{2:end});cat(1,subG2{2:end})];
    trainY = [ones(size(cat(1,subG1{2:end}),1),1);...
        zeros(size(cat(1,subG2{2:end}),1),1)];
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    pred = predict(mdl,testX);
    [thisX,thisY,~,aFSubRAE(ii)] = perfcurve(testY,pred,1);
    frSubXAE(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
        linspace(0,1,2401));
    frSubYAE(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
        linspace(0,1,2401));
end
%% Dual vs. ae
figure
hold on
plot(mean(mSubXAE),mean(mSubYAE),'k')
plot(mean(mrSubXAE),mean(mrSubYAE),'--k')
legend({['male: ',num2str(round(mean(mSubAAE),2)),'\pm',...
    num2str(round(conf(mSubAAE,0.95),2))],['animal detector: ',...
    num2str(round(mean(aMSubRAE),2)),'\pm',...
    num2str(round(conf(aMSubRAE,0.95),2))]},'location','se')
box off
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
title('male dual vs. ae; 4 vs. 4; GLM; LOO')

figure
hold on
plot(mean(fSubXAE),mean(fSubYAE),'k')
plot(mean(frSubXAE),mean(frSubYAE),'--k')
legend({['female: ',num2str(round(mean(fSubAAE),2)),'\pm',...
    num2str(round(conf(fSubAAE,0.95),2))],['animal detector: ',...
    num2str(round(mean(aFSubRAE),2)),'\pm',...
    num2str(round(conf(aFSubRAE,0.95),2))]},'location','se')
box off
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
title('female dual vs. ae; 4 vs. 4; GLM; LOO')
%% Male dual vs. mia
%% Compare 05, 10, and combined single feature models
% a05 = squeeze(mean(mean(sA05,1),2)).*sign(squeeze(mean(mean(sA05Sign))));
% a10 = squeeze(mean(mean(sA10,1),2)).*sign(squeeze(mean(mean(sA10Sign))));
% aComb = squeeze(mean(sA,1))'.*sign(squeeze(mean(sASign,1)))';
% feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
%     {'d','t','a','b','lg','hg'})';
% [~,sortA05Ind] = sort(abs(a05),'descend');
% sortA05 = a05(sortA05Ind);
% [~,sortA10Ind] = sort(abs(a10),'descend');
% sortA10 = a10(sortA10Ind);
% [~,sortACombInd] = sort(abs(aComb),'descend');
% sortAComb = aComb(sortACombInd);
% feat05 = feat(sortA05Ind);
% feat10 = feat(sortA10Ind);
% featComb = feat(sortACombInd);
% 
% figure
% hold on
% plot(a05(1:48),a10(1:48),'ks')
% plot(a05(49:end),a10(49:end),'ko')
% figure
% hold on
% plot(a10(1:48),aComb(1:48),'ks')
% plot(a10(49:end),aComb(49:end),'ko')
%% Build dual (off) vs. con (off) and dual (on) vs. con (off)
% Load output of runangelaMIAonoff.m
for ii = 1:100
    load(['F:/angelaAlcoholDual/offOff/dualOffConOff',num2str(ii),'.mat'])
    a05OffAll(ii,:) = a05off(ii,:);
    a05SingleOffAll(ii,:,:) = a05SingleOff(ii,:,:);
    a05SignOffAll(ii,:,:) = a05SignOff(ii,:,:);
    a05OnAll(ii,:) = a05on(ii,:);
    a05SingleOnAll(ii,:,:) = a05SingleOn(ii,:,:);
    a05SignOnAll(ii,:,:) = a05SignOn(ii,:,:);
    
    a10OffAll(ii,:) = a10off(ii,:);
    a10SingleOffAll(ii,:,:) = a10SingleOff(ii,:,:);
    a10SignOffAll(ii,:,:) = a10SignOff(ii,:,:);
    a10OnAll(ii,:) = a10on(ii,:);
    a10SingleOnAll(ii,:,:) = a10SingleOn(ii,:,:);
    a10SignOnAll(ii,:,:) = a10SignOn(ii,:,:);
    
    aBothOffAll(ii,:) = aBothOff(ii,:);
    aBothSingleOffAll(ii,:,:) = aBothSingleOff(ii,:,:);
    aBothSignOffAll(ii,:,:) = aBothSignOff(ii,:,:);
    aBothOnAll(ii,:) = aBothOn(ii,:);
    aBothSingleOnAll(ii,:,:) = aBothSingleOn(ii,:,:);
    aBothSignOnAll(ii,:,:) = aBothSignOn(ii,:,:);
end
%%
% Sort single features
nameVect = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'})';
a05OffMean = squeeze(mean(a05SingleOffAll.*a05SignOffAll,[1,2]));
[~,ind] = sort(abs(a05OffMean),'descend');
a05OffMeanSort = a05OffMean(ind);
a05OffFeat = nameVect(ind);

a05OnMean = squeeze(mean(a05SingleOnAll.*a05SignOnAll,[1,2]));
[~,ind] = sort(abs(a05OnMean),'descend');
a05OnMeanSort = a05OnMean(ind);
a05OnFeat = nameVect(ind);

a10OffMean = squeeze(mean(a10SingleOffAll.*a10SignOffAll,[1,2]));
[~,ind] = sort(abs(a10OffMean),'descend');
a10OffMeanSort = a10OffMean(ind);
a10OffFeat = nameVect(ind);

a10OnMean = squeeze(mean(a10SingleOnAll.*a10SignOnAll,[1,2]));
[~,ind] = sort(abs(a10OnMean),'descend');
a10OnMeanSort = a10OnMean(ind);
a10OnFeat = nameVect(ind);

aBothOffMean = squeeze(mean(aBothSingleOffAll.*aBothSignOffAll,[1,2]));
[~,ind] = sort(abs(aBothOffMean),'descend');
aBothOffMeanSort = aBothOffMean(ind);
aBothOffFeat = nameVect(ind);

aBothOnMean = squeeze(mean(aBothSingleOnAll.*aBothSignOnAll,[1,2]));
[~,ind] = sort(abs(aBothOnMean),'descend');
aBothOnMeanSort = aBothOnMean(ind);
aBothOnFeat = nameVect(ind);
%% 
figure
hold on
plot(squeeze(mean(x05,[1,2])),squeeze(mean(y05,[1,2])))
plot(squeeze(mean(x05on,[1,2])),squeeze(mean(y05on,[1,2])))
legend({['dual off vs. con off: ',num2str(round(mean(a05,[1,2]),2)),...
    '\pm',num2str(round(conf(reshape(a05,1,numel(a05)),0.95),3))],...
    ['dual on vs. con off: ',num2str(round(mean(a05On,[1,2]),2)),...
    '\pm',num2str(round(conf(reshape(a05On,1,numel(a05On)),0.95),3))]},...
    'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
title('5 mg/kg; built dual off vs. con off; LOO')
figure
hold on
plot(squeeze(mean(x10,[1,2])),squeeze(mean(y10,[1,2])))
plot(squeeze(mean(x10on,[1,2])),squeeze(mean(y10on,[1,2])))
legend({['dual off vs. con off: ',num2str(round(mean(a10,[1,2]),2)),...
    '\pm',num2str(round(conf(reshape(a10,1,numel(a10)),0.95),3))],...
    ['dual on vs. con off: ',num2str(round(mean(a10On,[1,2]),2)),...
    '\pm',num2str(round(conf(reshape(a10On,1,numel(a10On)),0.95),3))]},...
    'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
title('10 mg/kg; built dual off vs. con off; LOO')

figure
hold on
plot(squeeze(mean(xBoth,[1,2])),squeeze(mean(yBoth,[1,2])))
plot(squeeze(mean(xBothOn,[1,2])),squeeze(mean(yBothOn,[1,2])))
plot(squeeze(mean(xBothPerm,[1,2])),squeeze(mean(yBothPerm,[1,2])))
legend({['dual off vs. con off: ',num2str(round(mean(aBoth,[1,2]),2)),...
    '\pm',num2str(round(conf(reshape(aBoth,1,numel(aBoth)),0.95),3))],...
    ['dual on vs. con off: ',num2str(round(mean(aBothOn,[1,2]),2)),...
    '\pm',num2str(round(conf(reshape(aBothOn,1,numel(aBothOn)),0.95),3))...
    ],['permuted: ',num2str(round(mean(aBothPerm,[1,2]),2)),'\pm',...
    num2str(round(conf(reshape(aBothPerm,1,numel(aBothPerm)),0.95),3))]}...
    ,'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
title('5 & 10 mg/kg; built dual off vs. con off; LOO')
%% Build, and test, dual (on) vs. con (off)
for ii = 1:100
    for jj = 1:size(cmbsBoth,1)
        thisDual = cat(1,dualOn{cmbsBoth(jj,1),1:2});
        thisCon = conOff{cmbsBoth(jj,2)};
        testX = [thisDual(randperm(size(thisDual,1),1200),:);...
            thisCon(randperm(size(thisCon,1),1200),:)];
        testY = [ones(1200,1);zeros(1200,1)];
        dualTrainInds = logicFind(1,~ismember(1:9,cmbsBoth(jj,1)),'==');
        conTrainInds = logicFind(1,~ismember(1:7,cmbsBoth(jj,2)),'==');
        trainX = []; trainY = []; permTrainY = [];
        permDual = dualTrainInds(randperm(numel(dualTrainInds),...
            numel(dualTrainInds)/2));
        for k = dualTrainInds
            thisAnimal = cat(1,dualOn{k,1:2});
            trainX = [trainX;thisAnimal(randperm(size(...
                thisAnimal,1),1200),:)];
            trainY = [trainY;ones(1200,1)];
            if ismember(k,permDual)
               permTrainY =  [permTrainY;zeros(1200,1)];
            else
                permTrainY = [permTrainY;ones(1200,1)];
            end
        end
        permCon = conTrainInds(randperm(numel(conTrainInds),...
            numel(conTrainInds)/2));
        for k = conTrainInds
            trainX = [trainX;conOff{k}(randperm(size(...
                conOff{k},1),1200),:)];
            trainY = [trainY;zeros(1200,1)];
            if ismember(k,permCon)
               permTrainY =  [permTrainY;ones(1200,1)];
            else
                permTrainY = [permTrainY;zeros(1200,1)];
            end
        end
        % GLM
        mdl = fitglm(trainX,trainY,'distribution','binomial');
        pred = predict(mdl,testX);
        [thisX,thisY,~,a(ii,jj)] = perfcurve(testY,pred,1);
        x(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
            linspace(0,1,2401));
        y(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
            linspace(0,1,2401));
        % Single  feature
        for f = 1:216
            mdl = fitglm(trainX(:,f),trainY,'distribution','binomial');
            pred = predict(mdl,testX(:,f));
            [~,~,~,aSingle(ii,jj,f)] = perfcurve(testY,pred,1);
            aSign(ii,jj,f) = sign(table2array(mdl.Coefficients(2,1)));
        end
        % Permuted
%         mdl = fitglm(trainX,permTrainY,'distribution','binomial');
%         pred = predict(mdl,testX);
%         [thisX,thisY,~,aPerm(ii,jj)] = perfcurve(testY,pred,1);
%         xPerm(ii,jj,:) = interp1(linspace(0,1,numel(thisX)),thisX,...
%             linspace(0,1,2401));
%         yPerm(ii,jj,:) = interp1(linspace(0,1,numel(thisY)),thisY,...
%             linspace(0,1,2401));
    end
    disp(ii)
end
% save('F:/angelaAlcoholDual/dualOnConOff2.mat','x','y','a','aSingle',...
%     'aSign','aPerm','xPerm','yPerm')
%%
thisM = squeeze(mean(aSingle.*aSign,[1,2]));
[~,sInd] = sort(abs(thisM),'descend');
sA = thisM(sInd);
nameVect = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'})';
sortFeat = nameVect(sInd);
%%
figure
hold on
plot(squeeze(mean(x,[1,2])),squeeze(mean(y,[1,2])))
plot(squeeze(mean(xPerm,[1,2])),squeeze(mean(yPerm,[1,2])))
legend({['dual on vs. con off: ',num2str(round(mean(a,[1,2]),2)),'\pm',...
    num2str(round(conf(reshape(a,1,numel(a)),0.95),3))],...
    ['permuted: ',num2str(round(mean(aPerm,[1,2]),2)),'\pm',...
    num2str(round(conf(reshape(aPerm,1,numel(aPerm)),0.95),3))]},...
    'location','nw')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
title('5 & 10 mg/kg; built dual on vs. con off; LOO')
box off
%%
mInds = 6:16;
mSamps = cell2mat(onSamps(mInds,:));
otherSampes = cell2mat(onSamps([1:5,17:end],:));
%% subtract alcohol on from alcohol off (pre and post 45+ minutes)
offData = cell(3,39);
alcDiff = cell(3,39);
for ii = 1:3
    for jj = 1:numel(predata{ii})
        offData{ii,jj} = [predata{ii}{jj};post45{ii,jj}];
        alcDiff{ii,jj} = mean(offData{ii,jj})-alc30{ii,jj};
    end
end
% combine other groups and set aside mia (leave animals apart for easier
% train/test sets)
mia = alcDiff(:,6:16);
other = alcDiff(:,[1:5,17:end]);
% cheat hard code
last = [25,26,28];
% Leave one mia and two others out
for ii = 1:3
    for jj = 1:100
        miaInds = randperm(9);
        otherInds = randperm(last(ii));
        miaTrain = cat(1,mia{ii,miaInds(2:end)});
        otherTrain = cat(1,other{ii,otherInds(3:end)});
        thisTrainX = [miaTrain;otherTrain];
        thisTrainY = [ones(size(miaTrain,1),1);zeros(size(otherTrain,1),1)];
        miaTest = mia{ii,miaInds(1)};
        otherTest = cat(1,other{ii,otherInds(1:2)});
        thisTestX = [miaTest;otherTest];
        thisTestY = [ones(size(miaTest,1),1);zeros(size(otherTest,1),1)];
        
        mdl = fitglm(thisTrainX,thisTrainY,'distribution','binomial');
        prob = predict(mdl,thisTestX);
        [~,~,~,a(jj,ii)] = perfcurve(thisTestY,prob,1);
    end
end
%% load lasso data
for ii = 1:100
    load(['G:\GreenLab\data\angelaMiaAlcInj\onVoff\miaData',num2str(ii),...
        '.mat'],'a','aR','auc')
    allA(ii,:) = a;
    allAR(ii,:) = aR;
    singleA(ii,:,:) = auc;
end
%%
figure
hold on
histogram(allA(:,1),'binWidth',.01)
histogram(allA(:,2),'binWidth',.01)
histogram(allA(:,3),'binWidth',.01)

doubleHist(allA(:,3),allA(:,1))
title('alc 1 mg/kg vs. sal')
doubleHist(allA(:,2),allA(:,1))
title('alc 0.5 mg/kg vs. sal')
%% single features
nameVect = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'})';
mSingleFeat = squeeze(mean(singleA,1))';
%% angelaMiaAlcInj
for ii = 1:100
    load(['G:/GreenLab/data/angelaMiaAlcInj/model/',num2str(ii),'.mat'])
    aReal(ii) = a;
    aRand(ii) = aR;
end
doubleHist(aReal,aRand)
%%
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\angelaMIA\processed\'],{'mSaline';'mMIA';'fSaline';'fMIA'},...
    {'pow','coh'},'avg','rel');
mSal = cat(1,data{1,1}{:});
mP = cat(1,data{1,2}{:});
fSal = cat(1,data{1,3}{:});
fP = cat(1,data{1,4}{:});
save(['C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\analyzed\'...
    '2cohortAvgData.mat'],'data','files','samp','mSal','mP','fSal','fP')
%%
for ii = 1:100
    load(['C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\analyzed\'...
        'sexModels\',num2str(ii),'.mat'])
    fA(ii) = female.a;
    fAcc(ii,:) = 1-female.err;
    fARand(ii) = female.aRand;
    fAccRand(ii,:) = 1-female.errRand;
    fX(ii,:) = female.x;
    fY(ii,:) = female.y;
    fB(ii,:) = female.beta{1}.survBeta;
    fBS(ii,:) = female.beta{1}.signBeta;
    
    mA(ii) = male.a;
    mAcc(ii,:) = 1-male.err;
    mARand(ii) = male.aRand;
    mAccRand(ii,:) = 1-male.errRand;
    mX(ii,:) = male.x;
    mY(ii,:) = male.y;
    mB(ii,:) = male.beta{1}.survBeta;
    mBS(ii,:) = male.beta{1}.signBeta;
end
figure
hold on
plot(mean(mX,1),mean(mY,1))
plot(mean(fX,1),mean(fY,1))
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
box off
xlabel('FPR'); ylabel('TPR')
legend({['Male: ',num2str(round(mean(mA),2)),'\pm',...
    num2str(round(std(mA),2))];['Female: ',num2str(round(mean(fA),2)),...
    '\pm',num2str(round(std(fA),2))]},'location','se')
title('Saline vs. Poly I:C')

doubleHist(reshape(mAcc,1,10000),reshape(mAccRand,1,10000),...
    'xlab','Accuracy','main','Male');
doubleHist(reshape(fAcc,1,10000),reshape(fAccRand,1,10000),...
    'xlab','Accuracy','main','Female');
%% male female; logistics
for ii = 1:100
    disp(num2str(ii))
    sInds = randperm(size(mSal,1),3);
    pInds = randperm(size(mP,1),3);
    
    mTrainX = [mSal(~ismember(1:size(mSal,1),sInds),:);...
        mP(~ismember(1:size(mP,1),pInds),:)];
    mTrainY = [zeros(size(mSal,1)-3,1);ones(size(mP,1)-3,1)];
    mTestX = [mSal(sInds,:);mP(pInds,:)];
    mTestY = [0;0;0;1;1;1];
    
    sInds = randperm(size(mSal,1),3);
    pInds = randperm(size(mP,1),3);
    
    fTrainX = [fSal(~ismember(1:size(fSal,1),sInds),:);...
        fP(~ismember(1:size(fP,1),pInds),:)];
    fTrainY = [zeros(size(fSal,1)-3,1);ones(size(fP,1)-3,1)];
    fTestX = [fSal(sInds,:);fP(pInds,:)];
    fTestY = [0;0;0;1;1;1];
    for jj = 1:216
        mdl = fitglm(mTrainX(:,jj),mTrainY,'distribution','binomial');
        prob = predict(mdl,mTestX(:,jj));
        [mLogX{jj}(ii,:),mLogY{jj}(ii,:),~,mLogA(ii,jj)] = perfcurve(...
            mTestY,prob,1,'Tvals',0:.1:1,'UseNearest',0);
        
        mdl = fitglm(fTrainX(:,jj),fTrainY,'distribution','binomial');
        prob = predict(mdl,fTestX(:,jj));
        [fLogX{jj}(ii,:),fLogY{jj}(ii,:),~,fLogA(ii,jj)] = perfcurve(...
            fTestY,prob,1,'Tvals',0:.1:1,'UseNearest',0);
    end
end
fLogAM = mean(fLogA,1); 
mLogAM = mean(mLogA,1);
label = {'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'};
nameVect = names(label,{'d','t','a','b','lg','hg'});
[fFeatV,fInd] = sort(fLogAM','descend');
[mFeatV,mInd] = sort(mLogAM','descend');
fFeat = nameVect(fInd)';
mFeat = nameVect(mInd)';
mSurv = mean(mB,1);
fSurv = mean(fB,1);
mSurvS = mSurv(mInd)';
fSurvS = fSurv(fInd)';
% Get proxy sign of models by comparing means
for ii = 1:216
    mSignLog(ii,1) = mean(mP(:,ii))>mean(mSal(:,ii));
    fSignLog(ii,1) = mean(fP(:,ii))>mean(fSal(:,ii));
end
mSignLog = mSignLog(mInd);
fSignLog = fSignLog(fInd);
% Find tiers for log AUCs
mTier = tier(mLogA(:,mInd));
fTier = tier(fLogA(:,fInd));
%%
load('C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\analyzed\avgData.mat')
nameVect = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},{'d','t','a','b','lg','hg'});
%%
inds = 49:6:216;
s = data(1:15,:);
poly = data(16:36,:);
%%
for ii = 1:size(inds,2)
    [~,p(ii)] = ttest2(s(:,inds(ii)),poly(:,inds(ii)));
    if p(ii)*28 < 0.05
       figure
       plot(ones(1,size(poly,1)),poly(:,inds(ii)),'ob')
       hold on
       plot([0.75 1.25 ],repmat(mean(poly(:,inds(ii))),1,2),'-k')
       
       plot(zeros(1,size(s,1)),s(:,inds(ii)),'or')
       plot([-.25 .25],repmat(mean(s(:,inds(ii))),1,2),'-k')
       
       xlim([-0.5 1.5])
       title(nameVect{inds(ii)})
       set(gca,'XTick',[0 1],'XTickLabel',{'S','Poly'})
       box off
    end
end

