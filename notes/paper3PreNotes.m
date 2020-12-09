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
fNames = fileSearch(['F:\paper2\preBinge_240sec_redo\'],'.mat');
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
    notCoh = coh{end-1};
    notPow = psdTrls{end-1};
    notTrls = trls{end-1};
    
    preCoh = coh(1:end-2);
    prePow = psdTrls(1:end-2);
    preTrls = trls(1:end-2);
    % Grab times of binges
    load(['D:\paper2\toProcess\',names{ii,1},'_',names{ii,2},'.mat'],...
        'eventTs')
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
        if ~isempty(preTrls{k})
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
    end
    % Find notTrls that do not overlap with preTrls
    samps = [];
    for jj = 1:size(trls,2)
        if ~isempty(preTrls{jj})
        samps = [samps;trls{1,jj}.sampleinfo]; %#ok<AGROW>
        end
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
    
    save(['F:\paper3\preBingeCombinedNew\',names{ii,1},'_',names{ii,2},...
        '.mat'],'trls','psdTrls','coh')
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
    load(['D:\paper3\analyzed\preBinge240\',num2str(ii),'.mat'])
    lassoA(ii,:) = cellfun(@(x) x.auc,preBingeLasso);
    logA(ii,:) = cellfun(@(x) x.auc,preBingeLog);
    load(['D:\paper3\analyzed\preBinge240\',num2str(ii),'_Rand.mat'])
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