% eoi = base1,base2(wash),stim1
[data,samp1,files,time1] = collateData('D:\dualSite\processed\toUseSingle\',...
    {'IL','in';'NAcS','in';'OFC','in'},{'pow','coh'},'trl','');
% eoi = base,wash,stim
[data2,samp2,files2,time2] = collateData(['F:\lsdStim\processed\'...
    'base\'],{'mPFC','in';'NAcS','in';'OFC','in'},{'pow','coh'},'trl',...
    '');
% combine
for ii = 1:3
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
% grab lsd data; eoi = base,wash,stim
[lsdData,lsdSamp,lsdFiles,lsdTime] = collateData(['F:\lsdStim'...
    '\processed\postLSD\'],{'.mat'},{'pow','coh'},'trl','');
allLSDBase = cat(1,lsdData{1}{2:7,1});
allLSDStim = cat(1,lsdData{1}{2:7,3});
allLSDWash = cat(1,lsdData{1}{2:7,2});
lsdAbsTime = lsdTime.abs{1}(2:7,:);
lsdData = lsdData{1}(2:7,:);
lsdMinSamps = min(cellfun(@(x) size(x,1),lsdSamp{1}(2:7,:)),[],2);
lsdSamp = lsdSamp{1}(2:7,:);
%% get washout start times
for ii = 1:size(allFiles{1},1)
    split = strsplit(allFiles{1}{ii},'_');
    id(ii) = split(1);
    if ii <= 18
        sdir = 'D:\dualSite\processed\toUseSingle\';
    else
        sdir = 'G:\GreenLab\data\lsdStim\processed\base\';
    end
    load([sdir,allFiles{1}{ii}],'hist')
    if ii <= 18
        eI = eventInd(hist.eventTs,{'base2',[0 5]});
    else
        eI = eventInd(hist.eventTs,{'washout',[0 5]});
    end
    washStart(ii) = hist.eventTs.t{eI(1)};
end
lsdFiles = lsdFiles{1}(2:7,:);
for ii = 1:size(lsdFiles,1)
   load(['G:\GreenLab\data\lsdStim\processed\postLSD\',lsdFiles{ii}],'hist')
   lsdWashStart(ii) = hist.eventTs.t{5};
end
%% Washout vs. Baseline models (Con and LSD)
% LSD
[lsdBaseTrainX,lsdBaseTestX,lsdWashTrainX,lsdWashTestX,lsdTrainX,...
    lsdTrainY,lsdTestX,lsdTestY] = deal([]);
for ii = 1:6
    trainInds = randperm(lsdMinSamps(ii),floor(lsdMinSamps(ii)*.8));
    testInds = 1:lsdMinSamps(ii);
    testInds = testInds(~ismember(testInds,trainInds));
    
    lsdBaseTrainX = [lsdBaseTrainX;lsdData{ii,1}(trainInds,:)];
    lsdBaseTestX = [lsdBaseTestX;lsdData{ii,1}(testInds,:)];
    lsdWashTrainX = [lsdWashTrainX;lsdData{ii,2}(trainInds,:)];
    lsdWashTestX = [lsdWashTestX;lsdData{ii,2}(testInds,:)];
    
    lsdTrainX = [lsdTrainX;lsdBaseTrainX;lsdWashTrainX];
    lsdTestX = [lsdTestX;lsdBaseTestX;lsdWashTestX];
    lsdTrainY = [lsdTrainY;zeros(numel(trainInds),1);ones(numel(trainInds),1)];
    lsdTestY = [lsdTestY;zeros(numel(testInds),1);ones(numel(testInds),1)];
end
% CON
cmbs = nchoosek(1:52,6);
for ii = 1:nchoosek(52,6)
    s(ii) = sum(minSamps{1}(cmbs(ii,:)));
end
% get combinations that are within +/- 50 samples of lsd data (552)
theseCmbs = s>502 & s<602;
figure
hold on
histogram(s(theseCmbs),'edgecolor','r','binwidth',2)
histogram(s(~theseCmbs),'binwidth',2)
plot([552 552],get(gca,'ylim'),'--b')
cmbsInds = randperm();
for ii = 1:100
    
end
%% Build data sets of stim not (control and LSD)
for jj = 1:100
    thisTrainX = []; thisLSDTrainX = [];
    thisTrainY = []; thisLSDTrainY = [];
    thisTestX = []; thisLSDTestX = [];
    thisTestY =[]; thisLSDTestY = [];
    [thisWashTestX,thisLSDWashTestX,thisWashTestY,thisLSDWashTestY,...
        thisWashTime,thisLSDWashTime,thisLSDWashTestSynX,...
        thisLSDWashTestSynY] = deal(cell(1,12));
    for ii = 1:size(minSamps{1},1)
        % calculate test and train sizes
        testN = ceil(0.2*minSamps{1}(ii));
        trainN = floor(0.8*minSamps{1}(ii));
        % generate indices
        % base
        testIndBase = randperm(samps{1}(ii,1),testN);
        indsLeft = 1:samps{1}(ii,1);
        indsLeft = indsLeft(~ismember(indsLeft,testIndBase));
        trainIndBase = randperm(numel(indsLeft),trainN);
        trainIndBase = indsLeft(trainIndBase);
        % stim
        testIndStim = randperm(samps{1}(ii,3),testN);
        indsLeft = 1:samps{1}(ii,3);
        indsLeft = indsLeft(~ismember(indsLeft,testIndStim));
        trainIndStim = randperm(numel(indsLeft),trainN);
        trainIndStim = indsLeft(trainIndStim);
        % grab data
        thisTrainX = [thisTrainX;allData{1}{ii,1}(trainIndBase,:);allData{1}{ii,3}(trainIndStim,:)];
        thisTrainY = [thisTrainY;zeros(trainN,1);ones(trainN,1)];
        thisTestX = [thisTestX;allData{1}{ii,1}(testIndBase,:);allData{1}{ii,3}(testIndStim,:)];
        thisTestY = [thisTestY;zeros(testN,1);ones(testN,1)];
        % get times for washout data
        allWashTime = mean(absTime{1}{ii,2},2)-washStart(ii);
        % subset for every 5 seconds up to 1 minute
        c = 1;
        for k = 5:5:60
            theseWashTimeInds{c} = allWashTime<k & allWashTime>k-5;
            theseWashTime{c} = allWashTime(theseWashTimeInds{c});
            thisWashTime{c} = [thisWashTime{c};theseWashTime{c}];
            thisWashTestX{c} = [thisWashTestX{c};allData{1}{ii,1}(testIndBase,:);allData{1}{ii,2}(theseWashTimeInds{c},:)];
            thisWashTestY{c} = [thisWashTestY{c};zeros(testN,1);ones(sum(theseWashTimeInds{c}),1)];
            c = c+1;
        end
    end
    % store in array
    trainX{jj} = thisTrainX;
    trainY{jj} = thisTrainY;
    testX{jj} = thisTestX;
    testY{jj} = thisTestY;
    washTestX{jj} = thisWashTestX;
    washTestY{jj} = thisWashTestY;
    washTime{jj} = thisWashTime;
    % get LSD data
    for ii = 1:size(lsdMinSamps,1)
        testN = ceil(0.2*lsdMinSamps(ii));
        trainN = floor(0.8*lsdMinSamps(ii));
        % get base indices
        testIndLSDBase = randperm(size(lsdSamp{ii,1},1),testN);
        indsLeft = 1:size(lsdSamp{ii,1},1);
        indsLeft = indsLeft(~ismember(indsLeft,testIndLSDBase));
        trainIndLSDBase = randperm(numel(indsLeft),trainN);
        trainIndLSDBase = indsLeft(trainIndLSDBase);
        % get stim indices
        testIndLSDStim = randperm(size(lsdSamp{ii,3},1),testN);
        indsLeft = 1:size(lsdSamp{ii,3},1);
        indsLeft = indsLeft(~ismember(indsLeft,testIndLSDStim));
        trainIndLSDStim = randperm(numel(indsLeft),trainN);
        trainIndLSDStim = indsLeft(trainIndLSDStim);
         % grab data
        thisLSDTrainX = [thisLSDTrainX;lsdData{ii,1}(trainIndLSDBase,:);lsdData{ii,3}(trainIndLSDStim,:)];
        thisLSDTrainY = [thisLSDTrainY;zeros(trainN,1);ones(trainN,1)];
        thisLSDTestX = [thisLSDTestX;lsdData{ii,1}(testIndLSDBase,:);lsdData{ii,3}(testIndLSDStim,:)];
        thisLSDTestY = [thisLSDTestY;zeros(testN,1);ones(testN,1)];
        % get times for washout data
        allLSDWashTime = mean(lsdAbsTime{ii,2},2)-lsdWashStart(ii);
        % subset for every 5 seconds uo to 1 minute
        c = 1;
        for k = 5:5:60
            theseLSDWashTimeInds{c} = allLSDWashTime<k & allLSDWashTime>k-5;
            theseLSDWashTime{c} = allLSDWashTime(theseLSDWashTimeInds{c});
            thisLSDWashTime{c} = [thisLSDWashTime{c};theseLSDWashTime{c}];
            thisLSDWashTestX{c} = [thisLSDWashTestX{c};lsdData{ii,1}(testIndLSDBase,:);lsdData{ii,2}(theseLSDWashTimeInds{c},:)];
            thisLSDWashTestY{c} = [thisLSDWashTestY{c};zeros(testN,1);ones(sum(theseLSDWashTimeInds{c}),1)];
            c = c+1;
        end
    end
%     [out_featuresSyn,out_labelsSyn] = ADASYN(thisLSDWashTestX{c},thisLSDWashTestY{c},1);
%     synInd = randperm(size(out_featuresSyn,1),testN);
%     thisLSDWashTestSynX{c} = [thisLSDWashTestSynX{c};out_featuresSyn(synInd,:);lsdData{ii,1}(testIndLSDBase,:)];
%     thisLSDWashTestSynY{c} = [thisLSDWashTestSynY{c};out_labelsSyn(synInd,:);zeros(testN,1)];
    % store
    lsdTrainX{jj} = thisLSDTrainX;
    lsdTrainY{jj} = thisLSDTrainY;
    lsdTestX{jj} = thisLSDTestX;
    lsdTestY{jj} = thisLSDTestY;
    lsdWashTestX{jj} = thisLSDWashTestX;
    lsdWashTestY{jj} = thisLSDWashTestY;
    lsdWashTime{jj} = thisLSDWashTime;
end
% save('lsdStimNot.mat','trainX','trainY','testX','testY','washTestX','washTestY','washTime','lsdTrainX','lsdTrainY','lsdTestX','lsdTestY','lsdWashTestX','lsdWashTestY','lsdWashTime','-v7.3')
%% Open models
for ii = 1:100
    load(['G:\GreenLab\data\lsdStim\lsdStimNot\',num2str(ii),'.mat'],'cwA','lwA','clA','acc','lsdAcc','lsdHist','hist')
    allCLA(ii) = clA; 
    allCWA(ii,:) = cwA;
    allLWA(ii,:) = lwA;
    allCA(ii) = acc{1}.acc;
    allLA(ii) = lsdAcc{1}.acc;
    % apply lsd model to control data
    [predY] = cvglmnetPredict(lsdAcc{1}.mdl,testX{ii},['lambda_',lsdHist.cfg.minTerm],'response');
    [~,~,~,allLCA(ii)] = perfcurve(testY{ii},predY,1);
    % apply models to random
    [predY] = cvglmnetPredict(acc{1}.mdl,testX{ii},['lambda_',hist.cfg.minTerm],'response');
    [~,~,~,allCRA(ii)] = perfcurve(testY{ii}(randperm(numel(testY{ii}),numel(testY{ii}))),predY,1);
    % apply models to random
    [predY] = cvglmnetPredict(lsdAcc{1}.mdl,lsdTestX{ii},['lambda_',lsdHist.cfg.minTerm],'response');
    [~,~,~,allLRA(ii)] = perfcurve(lsdTestY{ii}(randperm(numel(lsdTestY{ii}),numel(lsdTestY{ii}))),predY,1);
end
% figure
% hold on
% shadedErrorBar(5:5:60,mean(allCWA,1),std(allCWA,[],1))
% shadedErrorBar(5:5:60,mean(allLWA,1),std(allLWA,[],1))
%% Single features
load('F:\lsdStim\lsdStimNot.mat')
for ii = 1:100
    disp(ii)
    for k = 1:216
        mdl = fitglm(trainX{ii}(:,k),trainY{ii},'distribution','binomial','binomialsize',size(trainX{ii},1));
        s(k,ii) = table2array(mdl.Coefficients(2,1));
        prob = predict(mdl,testX{ii}(:,k));
        [~,~,~,a(ii,k)] = perfcurve(testY{ii},prob,1);
        [~,~,~,aR(ii,k)] = perfcurve(testY{ii}(randperm(numel(testY{ii}),numel(testY{ii}))),prob,1);
        
        lsdMdl = fitglm(lsdTrainX{ii}(:,k),lsdTrainY{ii},'distribution','binomial','binomialsize',size(lsdTrainX{ii},1));
        lsdS(k,ii) = table2array(lsdMdl.Coefficients(2,1));
        prob = predict(lsdMdl,lsdTestX{ii}(:,k));
        [~,~,~,lA(ii,k)] = perfcurve(lsdTestY{ii},prob,1);
        [~,~,~,lAR(ii,k)] = perfcurve(testY{ii}(randperm(numel(lsdTestY{ii}),numel(lsdTestY{ii}))),prob,1);
%         for jj = 1:12
%             washProb{ii,jj}(k,:) = predict(mdl,washTestX{ii}{jj}(:,k));
%             [~,~,~,aW(ii,jj,k)] = perfcurve(washTestY{ii}{jj},washProb{ii,jj},1);
%             
%             lsdWashProb{ii,jj}(k,:) = predict(lsdMdl,lsdWashTestX{ii}{jj}(:,k));
%             [~,~,~,aLSDW(ii,jj,k)] = perfcurve(lsdWashTestY{ii}{jj},lsdWashProb{ii,jj},1);
%         end
    end
end
%%
[sA,sAI] = sort(mean(a,1),'descend');
[sLA,sLAI] = sort(mean(lA,1),'descend');
figure
hold on
plot(1:216,sA,'.k')
plot(1:216,sLA(sAI),'.r')
figure
hold on
plot(1:216,sA(sLAI),'.k')
plot(1:216,sLA,'.r')
figure
hold on
plot(1:216,sA,'.k')
plot(1:216,sLA,'.r')
%%
for ii = 1:216
    [~,p(ii)] = ttest2(a(:,ii),aR(:,ii));
    [~,lP(ii)] = ttest2(lA(:,ii),lAR(:,ii));
end
pAdj = p*432;
lPAdj = lP*432;
%%
bins = [0.55:0.05:0.75;0.6:0.05:0.8];
[conBins,lsdBins,conPowBins,conCohBins,lsdPowBins,lsdCohBins] = ...
    deal(zeros(1,size(bins,2)));
for ii = 1:size(bins,2)
    conBins(ii) = sum(mean(a,1)>bins(1,ii) & mean(a,1)<=bins(2,ii));
    lsdBins(ii) = sum(mean(lA,1)>bins(1,ii) & mean(lA,1)<=bins(2,ii));
    conPowBins(ii) =  sum(mean(a(:,1:48),1)>bins(1,ii) & mean(a(:,1:48),1)<=bins(2,ii));
    conCohBins(ii) = sum(mean(a(:,49:end),1)>bins(1,ii) & mean(a(:,49:end),1)<=bins(2,ii));
    lsdPowBins(ii) = sum(mean(lA(:,1:48),1)>bins(1,ii) & mean(lA(:,1:48),1)<=bins(2,ii));
    lsdCohBins(ii) = sum(mean(lA(:,49:end),1)>bins(1,ii) & mean(lA(:,49:end),1)<=bins(2,ii));
    for jj = 1:6
        conPowFreqBins(ii,jj) = sum(mean(a(:,1+(jj-1):6:48),1)>bins(1,ii) & mean(a(:,1+(jj-1):6:48),1)<=bins(2,ii));
        conCohFreqBins(ii,jj) = sum(mean(a(:,49+(jj-1):6:end),1)>bins(1,ii) & mean(a(:,49+(jj-1):6:end),1)<=bins(2,ii));
        lsdPowFreqBins(ii,jj) = sum(mean(lA(:,1+(jj-1):6:48),1)>bins(1,ii) & mean(lA(:,1+(jj-1):6:48),1)<=bins(2,ii));
        lsdCohFreqBins(ii,jj) = sum(mean(lA(:,49+(jj-1):6:end),1)>bins(1,ii) & mean(lA(:,49+(jj-1):6:end),1)<=bins(2,ii));
    end
end
conFreqBins = conPowFreqBins+conCohFreqBins;
lsdFreqBins = lsdPowFreqBins+lsdCohFreqBins;
figure
bar(categorical({'CON','LSD'}),[conBins;lsdBins],'stacked')
title('All features')
box off

figure
x = categorical({'CON Pow','LSD Pow','CON Coh','LSD Coh'});
x = reordercats(x,{'CON Pow','LSD Pow','CON Coh','LSD Coh'});
bar(x,[conPowBins;lsdPowBins;conCohBins;lsdCohBins],'stacked')
title('Feature type')
box off

figure
x = categorical({'CON \delta','LSD \delta','CON \theta','LSD \theta','CON \alpha','LSD \alpha','CON \beta','LSD \beta','CON l\gamma','LSD l\gamma','CON h\gamma','LSD h\gamma'});
x = reordercats(x,{'CON \delta','LSD \delta','CON \theta','LSD \theta','CON \alpha','LSD \alpha','CON \beta','LSD \beta','CON l\gamma','LSD l\gamma','CON h\gamma','LSD h\gamma'});
bar(x,reshape([conFreqBins;lsdFreqBins],5,12)','stacked')
title('All features by frequency')
box off
%%
conSPos = mean(s>0,2);
lsdSPos = mean(lsdS>0,2);

figure
hold on
plot(1:216,conSPos(sAI),'.k')
plot(1:216,lsdSPos(sAI),'.r')
%%
sA>sLA;
%% 
for ii = 1:100
   % grab washProbs of washout (9testY == 1)
   theseProbs(:,1,ii) = washProb{ii}(washTestY{ii} == 1);
   theseProbs(:,2,ii) = washTime{ii};
end
%% plotting stim effects through time
% average brain activity within animals
uID = unique(id);
xi = linspace(0,1,1000);
ti = 0:5000;
this = cell(numel(uID),3);
for ii = 1:numel(uID)
    these = logicFind(uID{ii},id,'==');
    for k = 1:3
        c = 1;
        for jj = these
            if numel(relTime{1}{jj,k})>2
                this{ii,k}(c,:) = interp1(mean(relTime{1}{jj,k},2),...
                    allData{1}{jj,k}(:,1),xi,'linear');
            end
            if k == 2
               absT{ii,c} = mean(absTime{1}{jj,k}-washStart(jj),2);
               absTi{ii}(c,:) = interp1(mean(absTime{1}{jj,k}-...
                   washStart(jj),2),allData{1}{jj,k}(:,1),ti,'linear');
            end
            c = c+1;
        end
    end
    figure
    plot(0:5000,mean(absTi{ii},1))
    hold on
    xlimit = get(gca,'xlim');
    plot(xlimit,[mean(mean(this{ii,1},'omitnan'),'omitnan') mean(mean(this{ii,1},'omitnan'),'omitnan')],'k')
    plot(xlimit,[mean(mean(this{ii,1},'omitnan'),'omitnan')+std(mean(this{ii,1},'omitnan')) mean(mean(this{ii,1},'omitnan'),'omitnan')+std(mean(this{ii,1},'omitnan'))],'--k')
    plot(xlimit,[mean(mean(this{ii,1},'omitnan'),'omitnan')-std(mean(this{ii,1},'omitnan')) mean(mean(this{ii,1},'omitnan'),'omitnan')-std(mean(this{ii,1},'omitnan'))],'--k')
    plot(xlimit,[mean(mean(this{ii,3},'omitnan'),'omitnan') mean(mean(this{ii,3},'omitnan'),'omitnan')],':r')
    xlim([0 500])
    title(uID{ii})
end
%% plot lsd stim effects
feature = 1;
for ii = 1:7
    load(['G:\GreenLab\data\lsdStim\processed\postLSD\',lsdFiles{1}{ii}],...
        'hist')
   figure
   hold on
   plot(mean(lsdTime.abs{1}{ii,2}-hist.eventTs.t{5},2),lsdData{1}{ii,2}(:,feature))
   xlim([0 500])
   plot(get(gca,'xlim'),[mean(lsdData{1}{ii,1}(:,feature)) mean(lsdData{1}{ii,1}(:,feature))],'k')
   plot(get(gca,'xlim'),[mean(lsdData{1}{ii,1}(:,feature))+std(lsdData{1}{ii,1}(:,feature)) mean(lsdData{1}{ii,1}(:,feature))+std(lsdData{1}{ii,1}(:,feature))],'--k')
   plot(get(gca,'xlim'),[mean(lsdData{1}{ii,1}(:,feature))-std(lsdData{1}{ii,1}(:,feature)) mean(lsdData{1}{ii,1}(:,feature))-std(lsdData{1}{ii,1}(:,feature))],'--k')
   plot(get(gca,'xlim'),[mean(lsdData{1}{ii,3}(:,feature)) mean(lsdData{1}{ii,3}(:,feature))],':r')
   title(lsdFiles{1}{ii})
end
%% model building
for ii = 1:3
    for jj = 1:100
        disp(jj)
        allBase = []; allStim = [];
        for k = 1:size(allData{ii},1)
            theseBaseInds = randperm(samps{ii}(k,1),minSamps{ii}(k));
            theseStimInds = randperm(samps{ii}(k,3),minSamps{ii}(k));
            allBase = [allBase;allData{ii}{k,1}(theseBaseInds,:)];
            allStim = [allStim;allData{ii}{k,3}(theseStimInds,:)];
        end
        trainN = ceil(0.8*size(allBase,1));
        testN = floor(0.2*size(allBase,1));
        
        inds = randperm(size(allBase,1));
        
        trainX = [allBase(inds(1:trainN),:);allStim(inds(1:trainN),:)];
        trainY = [zeros(trainN,1);ones(trainN,1)];
        
        testX = [allBase(inds(trainN+1:end),:);allStim(inds(trainN+1:end),:)];
        testY = [zeros(testN,1);ones(testN,1)];
        % single features
        for k = 1:216
            mdl = fitglm(trainX(:,k),trainY,'distribution','binomial',...
                'binomialsize',numel(trainY));
            prob = predict(mdl,testX(:,k));
            [~,~,~,a(ii,jj,k)] = perfcurve(testY,prob,1);
        end
        % lasso
%         cfg = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se',[]);
%         [alpha,lambda{ii,jj},beta,fits,acc,hist] = lassoNet(trainX,trainY,'binomial',...
%             'deviance',1,10,1,cfg);
    end
end
%% features
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
%%
mA = mean(a,2);
[smA,ind] = sort(mA,'descend');
sFeat = feat(ind);
base = cat(1,allData{:,1});
stim = cat(1,allData{:,2});
figure
scatter3(base(:,ind(1)),base(:,ind(2)),base(:,ind(3)),'.')
hold on
scatter3(stim(:,ind(1)),stim(:,ind(2)),stim(:,ind(3)),'.')
%%
xlabel(feat(ind(1))); ylabel(feat(ind(2))); zlabel(feat(ind(3)))
for ii = 1:8
    figure
    hold on
    histogram(base(:,ind(ii)))
    histogram(stim(:,ind(ii)))
    title(feat(ind(ii)))
end
%% Load population model
for ii = 1:100
    load(['G:\GreenLab\data\stimNot\pop2\',num2str(ii),'.mat'])
    allA(:,:,ii) = a;
    for jj = 1:3
        allAcc(jj,ii) = acc{jj}{1}.acc;
    end
end
mAllA = mean(allA,3);
for ii = 1:3
    [smAllA(ii,:),inds(ii,:)] = sort(mAllA(ii,:),'descend');
end
figure
pcolor(mAllA)
% imagesc(mAllA)
colormap viridis
ylim([0.5 1.5])
xlim([0 50])
%%
for ii= 1:100
    load(['G:\GreenLab\data\stimNot\pop2\',num2str(ii),'.mat'])
    [x(ii,:),y(ii,:),~,thisA(ii)] = perfcurve(hist{1}.cfg.naive.testY,acc{1}{1}.pred,1);
end
figure 
plot(mean(x,1),mean(y,1))
xlabel('FPR'); ylabel('TPR');
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
box off
%% Load individual models
allA = []; allS = [];
for ii = 1:100
    load(['G:\GreenLab\data\stimNot\ind\',num2str(ii),'.mat'])
    allA(:,:,:,ii) = a;
    allS(:,:,:,ii) = s;
    for jj = 1:3
        for k = logicFind(1,double(~cellfun(@isempty,acc(jj,:))),'==')
            allAcc(jj,k,ii) = acc{jj,k}{1}.acc;
        end
    end
end
allA(allA==0) = NaN;
mAllA = mean(allA,4,'omitnan');
for ii = 1:3
   figure
   colormap viridis
   imagesc(squeeze(mAllA(ii,:,:))')
   set(gca,'xticklabel',feat(1:6:216),'xtick',1:6:216)
   xtickangle(45)
end
%% look at all psds
sites = {'IL','NAcS','OFC'};
basePSD = cell(1,3); stimPSD = cell(1,3); washPSD = cell(1,3);
for ii = 1:3
    files = fileSearch('D:\dualSite\processed\toUseSingle\',sites{ii},...
        'in');
    for jj = 1:size(files,2)
        load(files{jj})
        if size(psdTrls{1,1}.Pow,1) == 8
            basePSD{ii} = cat(3,basePSD{ii},psdTrls{1,1}.Pow);
            washPSD{ii} = cat(3,washPSD{ii},psdTrls{1,2}.Pow);
            stimPSD{ii} = cat(3,stimPSD{ii},psdTrls{1,3}.Pow);
        end
    end
end
%%
basem = mean(basePSD{1},3);
stimm = mean(stimPSD{1},3);
washm = mean(washPSD{1},3);
figure
hold on
plot(basem(1,:))
plot(stimm(1,:))
plot(washm(1,:))