%% Collate raw data: LSD vs. Saline; acute effects from that days' baseline
[data,samp,files] = collateData(['F:\lsd\processed\imagAcute\'],{'lsd';...
    'sal'},{'pow','coh'},'trl','');
% For each animal, average that day's baseline and use to normalize post
allData = [data{1,1};data{1,2}];
mBase = cellfun(@(x) mean(x,1),allData(:,1),'UniformOutput',0);
for ii = 1:numel(mBase)
    dif{ii} = allData{ii,2}-repmat(mBase{ii},size(allData{ii,2},1),1);
end
allPostLSD = dif(1:numel(files{1,1}));
allPostSal = dif(numel(files{1,1})+1:end);
minSamp = min(cellfun(@(x) size(x,1),dif));
save('F:\lsd\normAcuteImag.mat','allPostLSD','allPostSal','minSamp','data','files')
% Modeling under runAcuteLSDvSaline.m
% cd D:\lsd\acuteLSDvSaline\
% x = []; y = [];
% xR = []; yR = [];
% for ii = 1:100
%    load(['acuteLSDvSal',num2str(ii),'.mat'],'acc','accR','hist','histR')
%    [x(ii,:),y(ii,:),~,allA(ii)] = perfcurve(hist.cfg.naive.testY,...
%        acc{1}.pred,1,'TVals',linspace(0,1,numel(acc{1}.pred)),...
%        'UseNearest',0);
%    [xR(ii,:),yR(ii,:),~,allAR(ii)] = perfcurve(histR.cfg.naive.testY,...
%        accR{1}.pred,1,'TVals',linspace(0,1,numel(accR{1}.pred)),...
%        'UseNearest',0);
% end
%% LOO simple logistic
x = []; y = []; a = [];
xP = []; yP = []; aP = [];
xS = []; yS = []; aS = [];
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
for n = 1:30
    c = 1;
    lsdN = size(allPostLSD,2);
    salN = size(allPostSal,2);
    for ii = 1:lsdN
        for jj = 1:salN
            inds(c,:) = [ii,jj];
            c = c+1;
        end
    end
    % LSD data
    trainLSD = logicFind(1,~ismember(1:lsdN,inds(n,1)),'==');
    testLSD = logicFind(0,~ismember(1:lsdN,inds(n,1)),'==');
    
    thisLSD24 = [];
    for ii = trainLSD
        rng(ii*n)
        thisLSD24 = [thisLSD24;allPostLSD{ii}(randperm(size(allPostLSD{ii},1),...
            minSamp),:)];
    end
    lsdTest = allPostLSD{testLSD}(randperm(size(allPostLSD{testLSD},1),...
            minSamp),:);
    % Saline data
    trainSal = logicFind(1,~ismember(1:salN,inds(n,2)),'==');
    testSal = logicFind(0,~ismember(1:salN,inds(n,2)),'==');
    thisSal = [];
    for ii = 1:numel(allPostSal)
        rng((ii+numel(allPostLSD))*n)
        thisSal = [thisSal;allPostSal{ii}(randperm(size(allPostSal{ii},1),...
            minSamp),:)];
    end
    salTest = allPostSal{testSal}(randperm(size(allPostSal{testSal},1),...
            minSamp),:);
    % Combine and build models
    trainX = [thisLSD24;thisSal];
    trainY = [ones(size(thisLSD24,1),1);zeros(size(thisSal,1),1)];
    testX = [lsdTest;salTest];
    testY = [ones(minSamp,1);zeros(minSamp,1)];
    
    mdl = fitglm(trainX(:,subInds),trainY,'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX(:,subInds));
    [x(n,:),y(n,:),~,a(n)] = perfcurve(testY,prob,1);
    beta(:,:,n) = table2array(mdl.Coefficients);
    % Permuted
    mdl = fitglm(trainX(:,subInds),trainY(randperm(numel(trainY),numel(trainY))),'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX(:,subInds));
    [xP(n,:),yP(n,:),~,aP(n)] = perfcurve(testY,prob,1);
    % Single feature models
    for jj = 1:60
        theseTrain = trainX(:,subInds);
        theseTest = testX(:,subInds);
        mdl = fitglm(theseTrain(:,jj),trainY,'distribution','binomial','binomialSize',numel(trainY));
        acuteBeta(jj,n) = table2array(mdl.Coefficients(2,1));
        prob = predict(mdl,theseTest(:,jj));
        [xS(n,jj,:),yS(n,jj,:),~,aS(n,jj)] = perfcurve(testY,prob,1);
        [xSP(n,jj,:),ySP(n,jj,:),~,aSP(n,jj)] = perfcurve(testY(randperm(numel(testY),numel(testY))),prob,1);
    end
end
%% Plot Log LOO
% load('F:\lsd\acuteLSDvSalineLogLOOImag.mat')
figure
hold on
plot(mean(x,1),mean(y,1),'-k')
plot(mean(xP,1),mean(yP,1),'--k')
xlabel('FPR'); ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
title('LSD v. Saline Acute: Log LOO')
legend({['Real: ',num2str(round(mean(a),2)),'\pm',...
    num2str(round(conf(a,0.95),2))],['Permuted: ',...
    num2str(round(mean(aP),2)),'\pm',num2str(round(conf(aP,0.95),2))]},...
    'location','se')
%% Single feature analysis - bar plot
mAS = mean(aS,1);
for ii = 1:60
    [~,p(ii)] = ttest(aS(:,ii)-0.5);
end
pAdj = p*216;
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
feat = feat(subInds);
sigFeat = feat(pAdj<=0.05);
pow = sum(pAdj(1:48)<=0.05);
coh = sum(pAdj(49:end)<=0.05);
for jj = 1:6
%     powFreq(jj) = sum(pAdj(1+(jj-1):6:48)<=0.05);
%     cohFreq(jj) = sum(pAdj(49+(jj-1):6:end)<=0.05);
    powFreq(jj) = sum(mAS(1+(jj-1):6:48)>=0.6);
    cohFreq(jj) = sum(pAdj(49+(jj-1):6:end)>=0.6);
end
figure
x = categorical({'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
x = reordercats(x,{'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
bar(x,[powFreq;cohFreq]','stacked')
box off
title('Acute features by frequency')
legend({'Power','Coherence'})
%% 24 hour effects: LSD vs. Saline
[data,samp,files] = collateData(['F:\lsd\processed\imag24Hr\'],{'lsd';...
    'sal'},{'pow','coh'},'trl','');
lsd24 = data{1,1};
sal24 = data{1,2};
minSamp = min([cellfun(@(x) size(x,1),lsd24);...
    cellfun(@(x) size(x,1),sal24)]);
save('F:\lsd\LSDvSal24HrImag.mat','lsd24','sal24','minSamp')
% Modeling under run24HrLSDvSaline.m
% cd D:\lsd\LSDvSal24Hr\
% x24 = []; y24 = [];
% x24R = []; y24R = [];
% for ii = 1:100
%    load(['LSDvSal24hr_',num2str(ii),'.mat'],'acc','accR','hist','histR')
%    [x24(ii,:),y24(ii,:),~,allA24(ii)] = perfcurve(hist.cfg.naive.testY,...
%        acc{1}.pred,1,'TVals',linspace(0,1,numel(acc{1}.pred)),...
%        'UseNearest',0);
%    [x24R(ii,:),y24R(ii,:),~,allA24R(ii)] = perfcurve(...
%        histR.cfg.naive.testY,accR{1}.pred,1,'TVals',...
%        linspace(0,1,numel(accR{1}.pred)),'UseNearest',0);
% end
%% 24 hour effects: drug vs. baseline days (LSD and saline separate)
% grab baseline files (using data with no behavior from water recordings)
load('F:\lsd\LSDvSal24HrImag.mat')
load('F:\lsd\normAcuteImag.mat','data')
% Load data in the same order as lsd24 and sal24 from above
[baseLSDData,baseLSDSamp,baseLSDFiles] = collateData(...
    'F:\lsd\processed\imagWater\',{'29','39','80','82','88','90'},...
    {'pow','coh'},'trl','');
ids = [29,39,80,82,88,90];
for ii = 1:numel(ids)
    inds = ~cellfun(@isempty,strfind(baseLSDFiles{1},['Water_',...
        num2str(ids(ii))]));
    allLSDBase{ii} = cat(1,baseLSDData{1}{inds},data{1}{ii,1});
    lsdDiff{ii} = lsd24{ii}-mean(allLSDBase{ii},1);
end
[baseSalData,baseSalSamp,baseSalFiles] = collateData(...
    'F:\lsd\processed\imagWater\',{'26';'30';'37';'38';'81'},...
    {'pow','coh'},'trl','');
for ii = 1:numel(baseSalData)
    allSalBase{ii} = cat(1,baseSalData{ii}{:},data{2}{ii,1});
    salDiff{ii} = sal24{ii}-mean(allSalBase{ii},1);
end
save('F:\lsd\lsdSal24Hr-BaseImag.mat','allLSDBase','allSalBase',...
    'baseLSDFiles','baseSalFiles','lsdDiff','salDiff','minSamp');
%% Base-LSD vs. Base-Sal
% load 'F:\lsd\lsdSal24Hr-BaseImag.mat'
c = 1; inds = [];
for ii = 1:6
    for jj = 1:5
        inds(c,:) = [ii,jj];
        c = c+1;
    end
end
salN = numel(salDiff);
lsdN = numel(lsdDiff);
% Angela headstage indices
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
for jj = 1:size(inds,1)
    disp(jj)
    trainSal = logicFind(1,~ismember(1:salN,inds(jj,2)),'==');
    testSal = logicFind(0,~ismember(1:salN,inds(jj,2)),'==');
    trainLSD = logicFind(1,~ismember(1:lsdN,inds(jj,1)),'==');
    testLSD = logicFind(0,~ismember(1:lsdN,inds(jj,1)),'==');
    
    thisSal = [];
    for ii = trainSal
        thisSal = [thisSal;salDiff{ii}(randperm(size(salDiff{ii},1),...
            minSamp),:)];
    end
    salTest = salDiff{testSal}(randperm(size(salDiff{testSal},1),...
        minSamp),:);
    thisLSD = [];
    for ii = trainLSD
        thisLSD = [thisLSD;lsdDiff{ii}(randperm(size(lsdDiff{ii},1),...
            minSamp),:)];
    end
    lsdTest = lsdDiff{testLSD}(randperm(size(lsdDiff{testLSD},1),...
        minSamp),:);
    % Combine
    trainX = [thisSal;thisLSD];
    trainY = [zeros(size(thisSal,1),1);ones(size(thisLSD,1),1)];
    testX = [salTest;lsdTest];
    testY = [zeros(minSamp,1);ones(minSamp,1)];
    mdl = fitglm(trainX(:,subInds),trainY,'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX(:,subInds));
    [x24(jj,:),y24(jj,:),~,a24(jj)] = perfcurve(testY,prob,1);
    % Permuted
    mdl = fitglm(trainX,trainY(randperm(numel(trainY),numel(trainY))),'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX);
    [xP24(jj,:),yP24(jj,:),~,aP24(jj)] = perfcurve(testY,prob,1);
end
% Plot
figure
hold on
plot(mean(x24,1),mean(y24,1),'-k')
plot(mean(xP24,1),mean(yP24,1),'--k')
xlabel('FPR'); ylabel('TPR')
title('Base-LSD vs. Base-Saline: Log LOO')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['Real: ',num2str(round(mean(a24),2)),'\pm',...
    num2str(round(conf(a24,0.95),3))],['Permuted: ',...
    num2str(round(mean(aP24),2)),'\pm',...
    num2str(round(conf(aP24,0.95),2))]},'location','se')
%% Plot in 3D feature space
% Get ordered list of single features
load('F:\lsd\acuteLSDvSalineLogLOOImag.mat','aS')
mAS = mean(aS,1);
[sorted,indSort] = sort(mAS,'descend');
% Subset features
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,...
    175:180];
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
feat = feat(subInds);
% Get acute data ready
load('F:\lsd\normAcuteImag.mat')

allSal = cat(1,allPostSal{:});
allSal = allSal(:,subInds);
allLSD = cat(1,allPostLSD{:});
allLSD = allLSD(:,subInds);

% Get 24Hr post data ready
load('lsdSal24Hr-BaseImag.mat')
postSal = cat(1,salDiff{:});
postLSD = cat(1,lsdDiff{:});

% Plot - all samples
figure
hold on
scatter3(allSal(:,indSort(1)),allSal(:,indSort(2)),allSal(:,indSort(3)),...
    '.k')
scatter3(allLSD(:,indSort(1)),allLSD(:,indSort(2)),allLSD(:,indSort(3)),...
    '.r')
scatter3(postSal(:,indSort(1)),postSal(:,indSort(2)),...
    postSal(:,indSort(3)),'.b')
scatter3(postLSD(:,indSort(1)),postLSD(:,indSort(2)),...
    postLSD(:,indSort(3)),'.g')
xlabel(feat{indSort(1)})
ylabel(feat{indSort(2)})
zlabel(feat{indSort(3)})
legend({'acuteSal','acuteLSD','24Sal','24LSD'})
% Plot - mean within animal
figure
hold on
scatter3(cellfun(@(x) mean(x(:,indSort(1))),allPostSal),...
    cellfun(@(x) mean(x(:,indSort(2))),allPostSal),...
    cellfun(@(x) mean(x(:,indSort(3))),allPostSal),'.k')
scatter3(cellfun(@(x) mean(x(:,indSort(1))),allPostLSD),...
    cellfun(@(x) mean(x(:,indSort(2))),allPostLSD),...
    cellfun(@(x) mean(x(:,indSort(3))),allPostLSD),'.r')
scatter3(cellfun(@(x) mean(x(:,indSort(1))),salDiff),...
    cellfun(@(x) mean(x(:,indSort(2))),salDiff),...
    cellfun(@(x) mean(x(:,indSort(3))),salDiff),'.b')
scatter3(cellfun(@(x) mean(x(:,indSort(1))),lsdDiff),...
    cellfun(@(x) mean(x(:,indSort(2))),lsdDiff),...
    cellfun(@(x) mean(x(:,indSort(3))),lsdDiff),'.g')
xlabel(feat{indSort(1)})
ylabel(feat{indSort(2)})
zlabel(feat{indSort(3)})
legend({'acuteSal','acuteLSD','24Sal','24LSD'})
view([163 26])
%% LSD Stim effects
% eoi = base1,base2(wash),stim1
% Raw
[data,samp1,files,time1] = collateData(['D:\dualSite\processed\'...
    'toUseSingleImag\'],{'IL','in'},{'pow','coh'},'trl','');

% eoi = base,wash,stim
% Raw
[data2,samp2,files2,time2] = collateData(['F:\lsdStim\processed\'...
    'baseImag\'],{'mPFC','in'},{'pow','coh'},'trl','');

% Combine
for ii = 1%:3
    allData{ii} = [data{1,ii};data2{1,ii}];
    allFiles{ii} = [files{ii,1};files2{ii,1}];
    relTime{ii} = [time1.rel{1,ii};time2.rel{1,ii}];
    absTime{ii} = [time1.abs{1,ii};time2.abs{1,ii}];   
    % check for any with fewer than 216 features
    feats = cellfun(@(x) size(x,2),allData{ii});
    allData{ii}(feats(:,1)<216,:) = [];
    allFiles{ii}(feats(:,1)<216,:) = [];
    relTime{ii}(feats(:,1)<216,:) = [];
    absTime{ii}(feats(:,1)<216,:) = [];
    samps{ii} = cellfun(@(x) size(x,1),allData{ii});
    % minSamps from base and wash
    minSamps{ii} = min(samps{ii}(:,[1,2]),[],2);
    % minSamps from base and stim
%     minSamps{ii} = min(samps{ii}(:,[1,3]),[],2);
end

% Only grab first 30 minutes of stim data
% minutes = 30;
% data30min = allData;
% for ii = 1:size(allData{1},1)
%     toGrab = absTime{1,1}{ii,3}(:,1)<=absTime{1}{ii,3}(1)+60*minutes;
%     data30min{1}{ii,3} = allData{1}{ii,3}(toGrab,:);    
% end

% Grab lsd data; eoi = base,wash,stim
% Imag
[lsdData,lsdSamp,lsdFiles,lsdTime] = collateData(['F:\lsdStim'...
    '\processed\postLSDImag\'],{'.mat'},{'pow','coh'},'trl','');
% % Grab only first 30 minutes of stim data
% minutes = 30;
% lsdData30min = lsdData;
% for ii = 1:size(lsdData{1},1)
%     toGrab = lsdTime.abs{1}{ii,3}(:,1)<=lsdTime.abs{1}{ii,3}(1)+60*minutes;
%     lsdData30min{1}{ii,3} = lsdData{1}{ii,3}(toGrab,:);
% end
% % Replace allData and lsdData with 30 minute versions
% allData = data30min; lsdData = lsdData30min;
% Calculate stim-baseline control vs. stim-baseline drug
ids = {'IRDM14';'IRDM15';'IRDM16';'IRDM21';'IRDM2_';'IRDM5';'IRDM6'};
thisBaseDiff = cell(1,numel(ids));
for ii = 1:numel(ids)
    theseFiles = logicFind(1,~cellfun(@isempty,strfind(allFiles{1},ids{ii})),'==');
    for jj = theseFiles
        thisBaseDiff{ii} = [thisBaseDiff{ii};allData{1}{jj,3}-mean(allData{1}{jj,1},1)];
    end
end
thisLSDDiff = cell(1,numel(ids));
for ii = 1:numel(ids)
    thisLSDDiff{ii} = lsdData{1}{ii,3}-mean(lsdData{1}{ii,1},1);
end
% Remove IRDM2 (animal 5) due to low samples and IRDM14 (animal 1) for too
% few features
thisBaseDiff = thisBaseDiff([2:4,6:7]);
thisLSDDiff = thisLSDDiff([2:4,6:7]);
minSamp = min([cellfun(@(x) size(x,1),thisBaseDiff),cellfun(@(x) size(x,1),thisLSDDiff)]);
%%
x = []; y = []; t = []; a = []; 
xS = []; yS = []; tS = []; aS = []; 
xP = []; yP = []; tP = []; aP = [];
xSP = []; ySP = []; tSP = []; aSP = [];
betaStim = [];
% Subset indices to match other electrode array
% subInds = [1:12,37:54,79:90,115:126,211:216];
subInds = 1:216;
for n = 1:100
    thisBase = []; thisLSD = [];
    for ii = 1:size(thisBaseDiff,2)
        thisBase = [thisBase;thisBaseDiff{ii}(randperm(size(thisBaseDiff{ii},1),minSamp),:)];
        thisLSD = [thisLSD;thisLSDDiff{ii}(randperm(size(thisLSDDiff{ii},1),minSamp),:)];
    end
    [trainX,trainY,testX,testY] = trainTest([thisBase;thisLSD],...
        [zeros(size(thisBase,1),1);ones(size(thisLSD,1),1)],0.20);
    mdl = fitglm(trainX(:,subInds),trainY,'distribution','binomial','binomialSize',numel(trainY));
%     betas(:,:,n) = table2array(mdl.Coefficients);
    prob = predict(mdl,testX(:,subInds));
    [x{n},y{n},t{n},a(n)] = perfcurve(testY,prob,1);
    % Permuted
    mdl = fitglm(trainX(:,subInds),trainY(randperm(numel(trainY),numel(trainY)),:),'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX(:,subInds));
    [xP{n},yP{n},tP{n},aP(n)] = perfcurve(testY,prob,1);
    c = 1;
    for ii = subInds
        mdl = fitglm(trainX(:,ii),trainY,'distribution','binomial','binomialSize',numel(trainY));
        betaStim(c,n) = table2array(mdl.Coefficients(2,1));
        prob = predict(mdl,testX(:,ii));
        [xS{n,c},yS{n,c},tS{n,c},aS(n,c)] = perfcurve(testY,prob,1);
        [xSP{n,c},ySP{n,c},tSP{n,c},aSP(n,c)]= perfcurve(testY(randperm(numel(testY),numel(testY))),prob,1);
        c  = c+1;
    end
end
% Single feature analysis
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'})';
% feat = names({'lmPFC','rmPFc','lNAcC','rNAcC'},...
%     {'d','t','a','b','lg','hg'})';
[sortA,sortInd] = sort(mean(aS,1)','descend');
sortFeat = feat(sortInd);
sortBeta = mean(betaStim(sortInd,:),2);
% Plot
% Interpolate x and y vectors to be the same size
for ii = 1:100
    interptX(ii,:) = interp1(linspace(0,1,numel(x{ii})),x{ii},...
        linspace(0,1,180));
    interptY(ii,:) = interp1(linspace(0,1,numel(y{ii})),y{ii},...
        linspace(0,1,180));
    interptXP(ii,:) = interp1(linspace(0,1,numel(xP{ii})),xP{ii},...
        linspace(0,1,180));
    interptYP(ii,:) = interp1(linspace(0,1,numel(yP{ii})),yP{ii},...
        linspace(0,1,180));
end
figure
hold on
plot(mean(interptX,1),mean(interptY,1),'-k')
plot(mean(interptXP,1),mean(interptYP,1),'--k')
xlabel('FPR'); ylabel('TPR')
title('Base-LSD vs. Base-Saline Stim: Log')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['Real: ',num2str(round(mean(a),2)),'\pm',...
    num2str(round(conf(a,0.95),3))],['Permuted: ',...
    num2str(round(mean(aP),2)),'\pm',...
    num2str(round(conf(aP,0.95),2))]},'location','se')
%% Actual leave one out (only five models) - shitty
x = []; y = []; t = []; a = []; 
xS = []; yS = []; tS = []; aS = []; 
xP = []; yP = []; tP = []; aP = [];
xSP = []; ySP = []; tSP = []; aSP = [];
betaStim = [];
% Subset indices to match other electrode array
% subInds = [1:12,37:54,79:90,115:126,211:216];
subInds = 1:216;

baseN = numel(thisBaseDiff);
lsdN = numel(thisLSDDiff);
cmbs = nchoosek(1:5,4);

for n = 1:5
    trainBaseX = []; trainLSDX = [];
    for ii = cmbs(n,:)
        trainBaseX = [trainBaseX;thisBaseDiff{ii}(randperm(size(thisBaseDiff{ii},1),minSamp),:)];
        trainLSDX = [trainLSDX;thisLSDDiff{ii}(randperm(size(thisLSDDiff{ii},1),minSamp),:)];
    end
    trainBaseY = zeros(size(trainBaseX,1),1);
    trainLSDY = ones(size(trainLSDX,1),1);
    
    lo = 6-n;
    testBaseX = thisBaseDiff{lo}(randperm(size(thisBaseDiff{lo},1),minSamp),:);
    testLSDX = thisLSDDiff{lo}(randperm(size(thisLSDDiff{lo},1),minSamp),:);
    testBaseY = zeros(size(testBaseX,1),1);
    testLSDY = ones(size(testLSDX,1),1);
    
    trainX = [trainBaseX;trainLSDX];
    trainY = [trainBaseY;trainLSDY];
    testX = [testBaseX;testLSDX];
    testY = [testBaseY;testLSDY];
    
    mdl = fitglm(trainX(:,subInds),trainY,'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX(:,subInds));
    [x{n},y{n},t{n},a(n)] = perfcurve(testY,prob,1);
    % Permuted
    mdl = fitglm(trainX(:,subInds),trainY(randperm(numel(trainY),numel(trainY)),:),'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX(:,subInds));
    [xP{n},yP{n},tP{n},aP(n)] = perfcurve(testY,prob,1);
    c = 1;
    for ii = subInds
        mdl = fitglm(trainX(:,ii),trainY,'distribution','binomial','binomialSize',numel(trainY));
        betaStim(c,n) = table2array(mdl.Coefficients(2,1));
        prob = predict(mdl,testX(:,ii));
        [xS{n,c},yS{n,c},tS{n,c},aS(n,c)] = perfcurve(testY,prob,1);
        [xSP{n,c},ySP{n,c},tSP{n,c},aSP(n,c)]= perfcurve(testY(randperm(numel(testY),numel(testY))),prob,1);
        c  = c+1;
    end
end
% Single feature analysis
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'})';
% feat = names({'lmPFC','rmPFc','lNAcC','rNAcC'},...
%     {'d','t','a','b','lg','hg'})';
[sortA,sortInd] = sort(mean(aS,1)','descend');
sortFeat = feat(sortInd);
sortBeta = mean(betaStim(sortInd,:),2);
% Plot

figure
hold on
plot(mean(cat(2,x{:}),2),mean(cat(2,y{:}),2),'-k')
plot(mean(cat(2,xP{:}),2),mean(cat(2,yP{:}),2),'--k')
xlabel('FPR'); ylabel('TPR')
title('Base-LSD vs. Base-Saline Stim: Log LOO')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['Real: ',num2str(round(mean(a),2)),'\pm',...
    num2str(round(conf(a,0.95),3))],['Permuted: ',...
    num2str(round(mean(aP),2)),'\pm',...
    num2str(round(conf(aP,0.95),2))]},'location','se')
%% Keep data sets separate; saline base v stim & lsd base v stim
% Rel
[data,samp1,files,time1] = collateData(['D:\dualSite\processed\'...
    'toUseSingleImag\'],{'IL','in'},{'pow','coh'},'trl','');
% Rel
[data2,samp2,files2,time2] = collateData(['F:\lsdStim\processed\'...
    'baseImag\'],{'mPFC','in'},{'pow','coh'},'trl','');
% Combine
for ii = 1%:3
    allData{ii} = [data{1,ii};data2{1,ii}];
    allFiles{ii} = [files{ii,1};files2{ii,1}];
    relTime{ii} = [time1.rel{1,ii};time2.rel{1,ii}];
    absTime{ii} = [time1.abs{1,ii};time2.abs{1,ii}];   
    % check for any with fewer than 216 features
    feats = cellfun(@(x) size(x,2),allData{ii});
    allData{ii}(feats(:,1)<216,:) = [];
    allFiles{ii}(feats(:,1)<216,:) = [];
    relTime{ii}(feats(:,1)<216,:) = [];
    absTime{ii}(feats(:,1)<216,:) = [];
    samps{ii} = cellfun(@(x) size(x,1),allData{ii});
    % minSamps from base and stim
    minSamps{ii} = min(samps{ii}(:,[1,3]),[],2);
end
% Raw
[lsdData,lsdSamp,lsdFiles,lsdTime] = collateData(['F:\lsdStim'...
    '\processed\postLSDImag\'],{'.mat'},{'pow','coh'},'trl','');
lsdData = lsdData{1}(:,[1,3]);
% manually set minSamp
minSamp = 80;
% Saline: base v. stim model
ids = {'IRDM14';'IRDM15';'IRDM16';'IRDM21';'IRDM2_';'IRDM5';'IRDM6'};
thisSaline = cell(numel(ids),2);
for ii = 1:numel(ids)
    theseFiles = logicFind(1,~cellfun(@isempty,strfind(allFiles{1},...
        ids{ii})),'==');
    for jj = theseFiles
        thisSaline{ii,1} = [thisSaline{ii,1};allData{1}{jj,1}];
        thisSaline{ii,2} = [thisSaline{ii,2};allData{1}{jj,3}];
    end
end
% Remove animals 1 and 5 (1 has missing features in LSD condition, and 
% 5 has too few samples in both) 
saline = thisSaline([2:4,6:7],:);

%%
[thisTrainX,thisTrainY,thisTestX,thisTestY] = deal(cell(1,6));
[xSaline,ySaline,tSaline] = deal(cell(1,100));
[aSaline] = deal(zeros(1,100));
[betaSaline,aSalineSingle] = deal(zeros(216,100));
[xSalineSingle,ySalineSingle,tSalineSingle] = deal(cell(216,100));
for ii = 1:100
    for jj = 1:5
        thisBase = saline{jj,1}(randperm(size(saline{jj,1},1),minSamp),:);
        thisStim = saline{jj,2}(randperm(size(saline{jj,2},1),minSamp),:);
        [thisTrainX{jj},thisTrainY{jj},thisTestX{jj},thisTestY{jj}] = ...
            trainTest([thisBase;thisStim],[zeros(minSamp,1);...
            ones(minSamp,1)],0.2);
    end
    % Normalize x values
    trainX = zscore(cat(1,thisTrainX{:}));
    trainY = cat(1,thisTrainY{:});
    testX = zscore(cat(1,thisTestX{:}));
    testY = cat(1,thisTestY{:});
%     [trainX,trainY,testX,testY] = trainTest([cat(1,saline{:,1});cat(1,saline{:,2})],[zeros(sum(cellfun(@(x) size(x,1),saline(:,1))),1);ones(sum(cellfun(@(x) size(x,1),saline(:,2))),1)],0.2);

    
    mdl = fitglm(trainX,trainY,'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX);
    [xSaline{ii},ySaline{ii},tSaline{ii},aSaline(ii)] = perfcurve(testY,prob,1);
    for n = 1:216
        mdl = fitglm(trainX(:,n),trainY,'distribution','binomial','binomialSize',numel(trainY));
        betaSaline(n,ii) = table2array(mdl.Coefficients(2,1));
        prob = predict(mdl,testX(:,n));
        [xSalineSingle{n,ii},ySalineSingle{n,ii},tSalineSingle{n,ii},aSalineSingle(n,ii)] = perfcurve(testY,prob,1);
    end
end
% LSD: base v. stim model
lsd = lsdData([2:4,6:7],:);
[thisTrainX,thisTrainY,thisTestX,thisTestY] = deal(cell(1,6));
[xLSD,yLSD,tLSD,xLSDP,yLSDP,tLSDP] = deal(cell(1,100));
[aLSD,aLSDP] = deal(zeros(1,100));
[betaLSD,aLSDsingle,aLSDsingleP] = deal(zeros(216,100));
[xLSDsingle,yLSDsingle,tLSDsingle,xLSDsingleP,yLSDsingleP,tLSDsingleP] = deal(cell(216,100));
for ii = 1:100
    for jj = 1:5
        thisBase = lsd{jj,1}(randperm(size(lsd{jj,1},1),minSamp),:);
        thisStim = lsd{jj,2}(randperm(size(lsd{jj,2},1),minSamp),:);
        [thisTrainX{jj},thisTrainY{jj},thisTestX{jj},thisTestY{jj}] = ...
            trainTest([thisBase;thisStim],[zeros(minSamp,1);...
            ones(minSamp,1)],0.2);
    end
    % Normalize x values
    trainX = zscore(cat(1,thisTrainX{:}));
    trainY = cat(1,thisTrainY{:});
    testX = zscore(cat(1,thisTestX{:}));
    testY = cat(1,thisTestY{:});
%     [trainX,trainY,testX,testY] = trainTest([cat(1,lsd{:,1});cat(1,lsd{:,2})],[zeros(sum(cellfun(@(x) size(x,1),lsd(:,1))),1);ones(sum(cellfun(@(x) size(x,1),lsd(:,2))),1)],0.2);
    
    mdl = fitglm(trainX,trainY,'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX);
    [xLSD{ii},yLSD{ii},tLSD{ii},aLSD(ii)] = perfcurve(testY,prob,1);
    [xLSDP{ii},yLSDP{ii},tLSDP{ii},aLSDP(ii)] = perfcurve(testY(randperm(numel(testY),numel(testY))),prob,1);
    for n = 1:216
        mdl = fitglm(trainX(:,n),trainY,'distribution','binomial','binomialSize',numel(trainY));
        betaLSD(n,ii) = table2array(mdl.Coefficients(2,1));
        prob = predict(mdl,testX(:,n));
        [xLSDsingle{n,ii},yLSDsingle{n,ii},tLSDsingle{n,ii},aLSDsingle(n,ii)] = perfcurve(testY,prob,1);
        
        prob = predict(mdl,testX(randperm(size(testX,1)),n));
        [xLSDsingleP{n,ii},yLSDsingleP{n,ii},tLSDsingleP{n,ii},aLSDsingleP(n,ii)] = perfcurve(testY,prob,1);
    end
end

% Single feature analysis
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
aLSDsingleM = mean(aLSDsingle,2);
[aLSDsingleMsort,lsdSortInd] = sort(aLSDsingleM,'descend');
lsdFeat = feat(lsdSortInd)';

aSalineSingleM = mean(aSalineSingle,2);
[aSalineSingleMsort,salineSortInd] = sort(aSalineSingleM,'descend');
salFeat = feat(salineSortInd)';
%% Plot
% Interpolate x and y vectors to be the same size
for ii = 1:100
    interptXsaline(ii,:) = interp1(linspace(0,1,numel(xSaline{ii})),xSaline{ii},...
        linspace(0,1,160));
    interptYsaline(ii,:) = interp1(linspace(0,1,numel(ySaline{ii})),ySaline{ii},...
        linspace(0,1,160));
    
    interptXlsd(ii,:) = interp1(linspace(0,1,numel(xLSD{ii})),xLSD{ii},...
        linspace(0,1,160));
    interptYlsd(ii,:) = interp1(linspace(0,1,numel(yLSD{ii})),yLSD{ii},...
        linspace(0,1,160));
end
figure
hold on
plot(mean(interptXsaline,1),mean(interptYsaline,1),'-k')
plot(mean(interptXlsd,1),mean(interptYlsd,1),'--k')
plot(mean(cat(2,xLSDP{:}),2),mean(cat(2,yLSDP{:}),2),':k')
xlabel('FPR'); ylabel('TPR')
title('Saline and LSD: Base vs Stim Log LOO')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['Saline: ',num2str(round(mean(aSaline),2)),'\pm',...
    num2str(round(conf(aSaline,0.95),3))],['LSD: ',...
    num2str(round(mean(aLSD),2)),'\pm',...
    num2str(round(conf(aLSD,0.95),3))],['Permuted: ',...
    num2str(round(mean(aLSDP),2)),'\pm',...
    num2str(round(conf(aLSDP,0.95),2))]},'location','se')
%% Single feature analysis
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
lsdMeanA = mean(aLSDsingle,2);
[lsdSortA,lsdSortInd] = sort(lsdMeanA,'descend');
lsdFeat = feat(lsdSortInd)';
betaLSDsort = mean(betaLSD(lsdSortInd,:),2);
[~,lsdP] = ttest2(aLSDsingle',aLSDsingleP'); 
lsdPadj = lsdP.*216;
lsdPadjSort = lsdPadj(lsdSortInd)';

salineMeanA = mean(aSalineSingle,2);
[salineSortA,salineSortInd] = sort(salineMeanA,'descend');
salineFeat = feat(salineSortInd)';
betaSalineSort = betaSaline(salineSortInd);
[~,salP] = ttest(aSalineSingle'-0.5);
salPadj = salP.*216;
salPadjSort = salPadj(salineSortInd)';

cutoffA = 0.55;
cutoffP = 0.05;
for jj = 1:6
    powFreqLSD(jj) = sum(lsdMeanA(1+(jj-1):6:48)>=cutoffA);
    cohFreqLSD(jj) = sum(lsdMeanA(49+(jj-1):6:end)>=cutoffA);
    
    powFreqSaline(jj) = sum(salineMeanA(1+(jj-1):6:48)>=cutoffA);
    cohFreqSaline(jj) = sum(salineMeanA(49+(jj-1):6:end)>=cutoffA);
    
    powFreqLSDP(jj) = sum(lsdPadj(1+(jj-1):6:48)<=cutoffP);
    cohFreqLSDP(jj) = sum(lsdP(49+(jj-1):6:end)<=cutoffP);
    
    powFreqSalineP(jj) = sum(salPadj(1+(jj-1):6:48)<=cutoffP);
    cohFreqSalineP(jj) = sum(salPadj(49+(jj-1):6:end)<=cutoffP);
end

figure
x = categorical({'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
x = reordercats(x,{'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
bar(x,[powFreqLSDP;cohFreqLSDP]','stacked')
box off
title('LSD stim features by frequency')
legend({'Power','Coherence'})

figure
x = categorical({'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
x = reordercats(x,{'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
bar(x,[powFreqSalineP;cohFreqSalineP]','stacked')
box off
title('Saline stim features by frequency')
legend({'Power','Coherence'})

%% Frequency networks
cmbs = nchoosek(1:8,2);
freqs = {'d','t','a','b','lg','hg'};
thisAS = mean(aS,1);
[~,p] = ttest2(aS,repmat(aP',1,60));
pAdj = p.*60;
for jj = 1:6
    adjMat = [];
    count = 48+jj;
    for ii = 1:size(cmbs,1)
        % Make sure beta coherence is significant
        if pAdj(count) <= 0.05
            adjMat(cmbs(ii,1),cmbs(ii,2)) =  thisAS(count);
        else
            adjMat(cmbs(ii,1),cmbs(ii,2)) = 0;
        end
        count = count+6;
    end
    adjMat = [adjMat;zeros(1,8)];
    % Find scaling factor based on minimum
    s = 1-min(min(abs(adjMat(adjMat~=0))));
    % Compute graph with scale factor
    g = graph(adjMat,{'lPFC','rPFC','lOFC','rOFC','lNAs','rNAs','lNAc','rNAc'},'upper');
    % Set up node coordinates based on regular octagon
    x = [1,sqrt(2)/2,0,-sqrt(2)/2,-1,-sqrt(2)/2,0,sqrt(2)/2];
    y = [0,sqrt(2)/2,1,sqrt(2)/2,0,-sqrt(2)/2,-1,-sqrt(2)/2];
    % Set edge color based on edge sign (+ = red; - = blue)
    edgeSign = g.Edges.Weight>=0;
    edgeColor = zeros(numel(edgeSign),3);
    edgeColor(logicFind(1,edgeSign,'=='),:) = repmat([1 0 0],sum(edgeSign),1);
    edgeColor(logicFind(0,edgeSign,'=='),:) = repmat([0 0 1],sum(~edgeSign),1);
    % Set node color based on node sign (+ = red; - = blue)
    nodeSign = thisAS(jj:6:48)>=0;
    nodeColor = zeros(numel(nodeSign),3);
    nodeColor(logicFind(1,nodeSign,'=='),:) = repmat([1 0 0],sum(nodeSign),1);
    nodeColor(logicFind(0,nodeSign,'=='),:) = repmat([0 0 1],sum(~nodeSign),1);
    % Check node significance
    sigNode = pAdj(jj:6:48) <= 0.05;
    % If not significant, set to white and size 10
    nodeColor(~sigNode,:) = repmat([1 1 1],sum(~sigNode),1);
    theseNodes = abs(thisAS(jj:6:48));
    theseNodes(~sigNode) = 10;
    % Plot
    figure('position',[895,400,667,571])
    h = plot(g,'XData',x,'YData',y,'markersize',theseNodes,'linewidth',abs(g.Edges.Weight)+s,'edgecolor',edgeColor,'nodecolor',nodeColor);
    title(freqs{jj})
    axis off
%     print(['LSD',freqs{jj},'.eps'],'-dwinc','-painters')
end
%% single feature analysis
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
feat = feat(subInds);

aSM = mean(aS,1);
[aSMsort,inds] = sort(aSM','descend');
featSort = feat(inds)';
betaM = mean(betaStim,2);
betaMS = betaM(inds);

for ii = 1:60
    [~,p(ii)] = ttest(aS(:,ii)-0.5);
end
pAdj = p*60;

sigFeat = feat(pAdj<=0.05);
pow = sum(pAdj(1:24)<=0.05);
coh = sum(pAdj(25:end)<=0.05);
for jj = 1:6
%     powFreq(jj) = sum(pAdj(1+(jj-1):6:48)<=0.05);
%     cohFreq(jj) = sum(pAdj(49+(jj-1):6:end)<=0.05);
    powFreq(jj) = sum(aSM(1+(jj-1):6:48)>=0.6);
    cohFreq(jj) = sum(aSM(49+(jj-1):6:end)>=0.6);
end
figure
x = categorical({'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
x = reordercats(x,{'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
bar(x,[powFreq;cohFreq]','stacked')
box off
title('Stim features by frequency')
legend({'Power','Coherence'})
%% LSD stim + IRDM
% Skip IRDM36 since no corresponding baseline stim IRDM recordings
[lsdIRDM,lsdIRDMsamps,lsdIRDMfiles] = collateData('F:\lsdIRDM\processed\',{'27';'32';'34';'37';'41'},{'pow','coh'},'trl','');
[baseIRDM,baseIRDMsamps,baseIRDMfiles] = collateData('F:\lsdIRDM\processed\control\',{'27';'32';'34';'37';'41'},{'pow','coh'},'trl','');
% Will return an error since IRDM41 only has 2 LSD recordings
for ii = 1:5
    for jj = 1:3
        baseDiff{ii,jj} = mean(baseIRDM{ii}{jj,1})-baseIRDM{ii}{jj,3};
        lsdDiff{ii,jj} = mean(lsdIRDM{ii}{jj,1})-lsdIRDM{ii}{jj,3};
    end
end
minSamp = min([sum(cellfun(@(x) size(x,1),baseDiff),2);sum(cellfun(@(x) size(x,1),lsdDiff),2)]);
% Concatenate across days
for ii = 1:5
    allLSDDiff{ii} = cat(1,lsdDiff{ii,:});
    allBaseDiff{ii} = cat(1,baseDiff{ii,:});
end
%% Combine and build models
for ii = 1:5
    testInd = ii;
    trainInd = 1:5;
    trainInd = trainInd(~ismember(trainInd,testInd));
    thisLSD = []; thisBase = [];
    for jj = trainInd
        thisLSD = [thisLSD;allLSDDiff{jj}(randperm(size(allLSDDiff{jj},1),minSamp),:)];
        thisBase = [thisBase;allBaseDiff{jj}(randperm(size(allBaseDiff{jj},1),minSamp),:)];
    end
    trainX = [thisLSD;thisBase];
    trainY = [ones(size(thisLSD,1),1);zeros(size(thisBase,1),1)];
    testX = [allLSDDiff{testInd}(randperm(size(allLSDDiff{testInd},1),minSamp),:);allBaseDiff{testInd}(randperm(size(allBaseDiff{testInd},1),minSamp),:)];
    testY = [ones(minSamp,1);zeros(minSamp,1)];
    
    mdl = fitglm(trainX,trainY,'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX);
    [x(ii,:),y(ii,:),t(ii,:),a(ii)] = perfcurve(testY,prob,1);
    % Permuted
    mdl = fitglm(trainX,trainY(randperm(numel(trainY),numel(trainY))),'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX);
    [xP(ii,:),yP(ii,:),tP(ii,:),aP(ii)] = perfcurve(testY,prob,1);
    % Single feature
%     for jj = 1:216
%         mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial','binomialSize',numel(trainY));
%         prob = predict(mdl,testX(:,jj));
%         [xS(ii,jj,:),yS(ii,jj,:),tS(ii,jj,:),aS(ii,jj)] = perfcurve(testY,prob,1);
%     end
end
%% PCA
% [coeff,score,latent,~,explained] = pca(allSal);
% allLSDpca = bsxfun(@minus,allLSD,mean(allLSD))/coeff';
% 
% figure
% hold on
% scatter3(score(:,1),score(:,2),score(:,3),'.')
% scatter3(allLSDpca(:,1),allLSDpca(:,2),allLSDpca(:,3),'.')
% scatter3(mean(score(:,1)),mean(score(:,2)),mean(score(:,3)),'ob')
% scatter3(mean(allLSDpca(:,1)),mean(allLSDpca(:,2)),mean(allLSDpca(:,3)),'or')
%% LOO simple logisitic
% load('F:\lsd\LSDvSal24HrRelImag.mat')
% x24 = []; y24 = []; a24 = [];
% subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
% for n = 1:30
%     c = 1;
%     lsdN = size(lsd24,1);
%     salN = size(sal24,1);
%     for ii = 1:lsdN
%         for jj = 1:salN
%             inds(c,:) = [ii,jj];
%             c = c+1;
%         end
%     end
%     % LSD data
%     trainLSD = logicFind(1,~ismember(1:lsdN,inds(n,1)),'==');
%     testLSD = logicFind(0,~ismember(1:lsdN,inds(n,1)),'==');
%     
%     thisLSD24 = [];
%     for ii = trainLSD
%         rng(ii*n)
%         thisLSD24 = [thisLSD24;lsd24{ii}(randperm(size(lsd24{ii},1),...
%             minSamp),:)];
%     end
%     lsdTest = lsd24{testLSD}(randperm(size(lsd24{testLSD},1),...
%             minSamp),:);
%     % Saline data
%     trainSal = logicFind(1,~ismember(1:salN,inds(n,2)),'==');
%     testSal = logicFind(0,~ismember(1:salN,inds(n,2)),'==');
%     thisSal = [];
%     for ii = 1:numel(sal24)
%         rng((ii+numel(sal24))*n)
%         thisSal = [thisSal;sal24{ii}(randperm(size(sal24{ii},1),...
%             minSamp),:)];
%     end
%     salTest = sal24{testSal}(randperm(size(sal24{testSal},1),...
%             minSamp),:);
%     % Combine and build models
%     trainX = [thisLSD24;thisSal];
%     trainY = [ones(size(thisLSD24,1),1);zeros(size(thisSal,1),1)];
%     testX = [lsdTest;salTest];
%     testY = [ones(minSamp,1);zeros(minSamp,1)];
%     
%     mdl = fitglm(trainX(:,subInds),trainY,'distribution','binomial',...
%         'binomialSize',numel(trainY));
%     prob = predict(mdl,testX(:,subInds));
%     [x24(n,:),y24(n,:),~,a24(n)] = perfcurve(testY,prob,1);
% %     beta24(:,:,n) = table2array(mdl.Coefficients);
%     % Permuted
%     mdl = fitglm(trainX(:,subInds),trainY(randperm(numel(trainY),numel(trainY))),...
%         'distribution','binomial','binomialSize',numel(trainY));
%     prob = predict(mdl,testX(:,subInds));
%     [x24P(n,:),y24P(n,:),~,a24P(n)] = perfcurve(testY,prob,1);
% end
% %%
% % Plot
% load('D:\lsd\LSDvSaline24HrLogLOO.mat')
% figure
% hold on
% plot(mean(x24,1),mean(y24,1),'-k')
% plot(mean(x24P,1),mean(y24P,1),'--k')
% xlabel('FPR'); ylabel('TPR')
% title('LSD v. Saline 24 Hr: Log LOO')
% set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
% legend({['Real: ',num2str(round(mean(a24),2)),'\pm',...
%     num2str(round(conf(a24,0.95),2))],['Permuted: ',...
%     num2str(round(mean(a24P),2)),'\pm',...
%     num2str(round(conf(a24P,0.95),2))]},'location','se')
%% Modeling under run24HrDrugvBaseline.m
% %% Simple log LOO
% load('F:\lsd\lsdSalBaseline.mat')
% load('F:\lsd\LSDvSal24HrRel.mat')       
% % LSD data
% lsd24N = numel(lsd24);
% lsdBaseN = numel(allLSDBase);
% c = 1; inds = [];
% for ii = 1:lsd24N
%     for jj = 1:lsdBaseN
%         inds(c,:) = [ii,jj];
%         c = c+1;
%     end
% end
% lsdX = []; lsdY = []; lsdA = []; betaLSD = [];
% lsdXP = []; lsdYP = []; lsdAP = [];
% for n = 1:6%size(inds,1)
%     trainLSD24 = logicFind(1,~ismember(1:lsd24N,inds(n,1)),'==');
%     testLSD24 = logicFind(0,~ismember(1:lsd24N,inds(n,1)),'==');
%     
%     trainLSDBase = logicFind(1,~ismember(1:lsdBaseN,inds(n,1)),'==');
%     testLSDBase = logicFind(0,~ismember(1:lsdBaseN,inds(n,1)),'==');
%     % LSD 24
%     thisLSD24 = [];
%     for ii = trainLSD24
%         thisLSD24 = [thisLSD24;lsd24{ii}(randperm(size(lsd24{ii},1),...
%             minSamp),:)];
%     end
%     lsd24Test = lsd24{testLSD24}(randperm(size(lsd24{testLSD24},1),...
%         minSamp),:);
%     % LSD Base
%     thisLSDBase = [];
%     for ii = trainLSDBase
%         thisLSDBase = [thisLSDBase;allLSDBase{ii}(randperm(size(...
%             allLSDBase{ii},1),minSamp),:)];
%     end
%     lsdBaseTest = allLSDBase{testLSDBase}(randperm(size(...
%         allLSDBase{testLSDBase},1),minSamp),:);
%     % Combine and build models
%     trainLSDX{n} = [thisLSD24;thisLSDBase];
%     trainLSDY{n} = [ones(size(thisLSD24,1),1);zeros(size(thisLSDBase,1),1)];
%     testLSDX{n} = [lsd24Test;lsdBaseTest];
%     testLSDY{n} = [ones(minSamp,1);zeros(minSamp,1)];
%     
%     mdlLSD{n} = fitglm(trainLSDX{n},trainLSDY{n},'distribution','binomial',...
%         'binomialSize',numel(trainLSDY{n}));
%     prob = predict(mdlLSD{n},testLSDX{n});
%     [lsdX(n,:),lsdY(n,:),~,lsdA(n)] = perfcurve(testLSDY{n},prob,1);
%     betaLSD(:,:,n) = table2array(mdlLSD{n}.Coefficients);
%     % Permuted
%     mdl = fitglm(trainLSDX{n},trainLSDY{n}(randperm(numel(trainLSDY{n}),numel(trainLSDY{n}))),...
%         'distribution','binomial','binomialSize',numel(trainLSDY{n}));
%     prob = predict(mdl,testLSDX{n});
%     [lsdXP(n,:),lsdYP(n,:),~,lsdAP(n)] = perfcurve(testLSDY{n},prob,1);
% end
% % Saline data
% sal24N = numel(sal24);
% salBaseN = numel(allSalBase);
% c = 1; inds = [];
% for ii = 1:sal24N
%     for jj = 1:salBaseN
%         inds(c,:) = [ii,jj];
%         c = c+1;
%     end
% end
% salX = []; salY = []; salA = []; betaSal = [];
% salXP = []; salYP = []; salAP = [];
% for n = 1:5%size(inds,1)
%     trainSal24 = logicFind(1,~ismember(1:sal24N,inds(n,1)),'==');
%     testSal24 = logicFind(0,~ismember(1:sal24N,inds(n,1)),'==');
%     
%     trainSalBase = logicFind(1,~ismember(1:salBaseN,inds(n,1)),'==');
%     testSalBase = logicFind(0,~ismember(1:salBaseN,inds(n,1)),'==');
%     % Sal 24
%     thisSal24 = [];
%     for ii = trainSal24
%         thisSal24 = [thisSal24;sal24{ii}(randperm(size(sal24{ii},1),...
%             minSamp),:)];
%     end
%     sal24Test = sal24{testSal24}(randperm(size(sal24{testSal24},1),...
%         minSamp),:);
%     % Sal base
%     thisSalBase = [];
%     for ii = trainSalBase
%         thisSalBase = [thisSalBase;allSalBase{ii}(randperm(size(...
%             allSalBase{ii},1),minSamp),:)];
%     end
%     salBaseTest = allSalBase{testSalBase}(randperm(size(...
%         allSalBase{testSalBase},1),minSamp),:);
%     % Combine and build models
%     trainSalX{n} = [thisSal24;thisSalBase];
%     trainSalY{n} = [ones(size(thisSal24,1),1);zeros(size(thisSalBase,1),1)];
%     testSalX{n} = [sal24Test;salBaseTest];
%     testSalY{n} = [ones(minSamp,1);zeros(minSamp,1)];
%     
%     mdlSal{n} = fitglm(trainSalX{n},trainSalY{n},'distribution','binomial',...
%         'binomialSize',numel(trainSalY{n}));
%     prob = predict(mdlSal{n},testSalX{n});
%     [salX(n,:),salY(n,:),~,salA(n)] = perfcurve(testSalY{n},prob,1);
%     betaSal(:,:,n) = table2array(mdlSal{n}.Coefficients);
%     % Permuted
%     mdl = fitglm(trainSalX{n},trainSalY{n}(randperm(numel(trainSalY{n}),numel(trainSalY{n}))),...
%         'distribution','binomial','binomialSize',numel(trainSalY{n}));
%     prob = predict(mdl,testSalX{n});
%     [salXP(n,:),salYP(n,:),~,salAP(n)] = perfcurve(testSalY{n},prob,1);
% end
% % Apply lsd mdl to sal and vice versa
% c = 1;
% for ii = 1:numel(mdlSal)
%     for jj = 1:numel(mdlLSD)
%         % LSD model to saline data
%         prob = predict(mdlLSD{jj},testSalX{ii});
%         [lsdSalX(c,:),lsdSalY(c,:),~,lsdSalA(c)] = perfcurve(testSalY{ii},prob,1);
%         
%         % saline model to LSD data
%         prob = predict(mdlSal{ii},testLSDX{jj});
%         [salLSDX(c,:),salLSDY(c,:),~,salLSDA(c)] = perfcurve(testLSDY{jj},prob,1);
%         c = c+1;
%     end
% end
% % Plot
% % LSD
% figure
% hold on
% plot(mean(lsdX,1),mean(lsdY,1),'-k')
% plot(mean(lsdXP,1),mean(lsdYP,1),'--k')
% plot(mean(salLSDX,1),mean(salLSDY,1),'-.k')
% xlabel('FPR'); ylabel('TPR')
% title('LSD 24 Hr vs. Baseline: Log LOO')
% set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
% legend({['Real: ',num2str(round(mean(lsdA),2)),'\pm',...
%     num2str(round(conf(lsdA,0.95),2))],['Permuted: ',...
%     num2str(round(mean(lsdAP),2)),'\pm',...
%     num2str(round(conf(lsdAP,0.95),2))],['Saline Model: ',...
%     num2str(round(mean(salLSDA),2)),'\pm',...
%     num2str(round(conf(salLSDA,0.95),2))]},'location','se')
% % Saline
% figure
% hold on
% plot(mean(salX,1),mean(salY,1),'-k')
% plot(mean(salXP,1),mean(salYP,1),'--k')
% plot(mean(lsdSalX,1),mean(lsdSalY,1),'-.k')
% xlabel('FPR'); ylabel('TPR')
% title('Sal 24 Hr vs. Baseline: Log LOO')
% set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
% legend({['Real: ',num2str(round(mean(salA),2)),'\pm',...
%     num2str(round(conf(salA,0.95),2))],['Permuted: ',...
%     num2str(round(mean(salAP),2)),'\pm',...
%     num2str(round(conf(salAP,0.95),2))],['LSD Model: ',...
%     num2str(round(mean(lsdSalA),2)),'\pm',...
%     num2str(round(conf(lsdSalA,0.95),2))]},'location','se')

