%% fix event labels for stim
files = fileSearch('G:/GreenLab/data/irdmNew/split','stim');
for ii = 2:size(files,2)
    disp(ii)
    load(files{ii})
    oldEventTs = eventTs;
    eventTs.t = oldEventTs.t([3:10,1,2]);
    eventTs.label = oldEventTs.label([3:10,1,2]);
    save(files{ii},'adfreq','eventTs','LFPTs','pl2')
end
%% get behavioral data
% files = fileSearch('E:/processed/ddt/','.mat');
% files = fileSearch('G:\GreenLab\data\irdmNew\processedContinuous\','.mat');
files = fileSearch('F:\irdmRound2\processedContinuous3\checkTrial\','.mat');
for ii = 385:size(files,2)
    disp(ii)
    disp(files{ii})
    cd('F:\irdmRound2\processedContinuous3\checkTrial\')
    load(files{ii})
    % only process files that haven't already been processed
    if ~exist('trials','var')
%         if isempty(eventTs.t{1})
%            keyboard
%            eventTs.t = eventTs.t(3:10);
%            eventTs.label = eventTs.label(3:10);
%         end
        cd('E:\processedMat\')
        load([files{ii}(1:end-8),'.mat'],'trials')
        trialsOld = trials;
        % Calculate new trials
        trials = ddtTrials(hist.eventTs,LFPTs,0);
        % check equal trials
        if size(trials,1) ~= size(trialsOld,1)
            trials = ddtTrials(hist.eventTs,LFPTs,1);
            warning('Inequal number of trials')
            keyboard
        end
        % check equal auc
        if trials(1).auc ~= trialsOld(1).auc
            trials = ddtTrials(hist.eventTs,LFPTs,1);
            warning('Inequal aucs')
            keyboard
        end
%         for jj = 1:size(trialsOld,1)
%             trials(jj).trapzAUC = trialsOld(jj).trapzAUC;
%             trials(jj).delay = trialsOld(jj).delay;
%         end
        save(['F:\irdmRound2\processedContinuous3\checkTrial\',files{ii}],'-append','trials','hist')
    end
    clear adfreq eventTs LFPTs pl2 trials
end
%% calculate auc using delays and trapz
% files = fileSearch('E:/processed/ddt/0-1-3-6-12_delay/','.mat')';
delays = [0,1,3,6,12];
% files = fileSearch('E:/processed/ddt/0-2.5-5-7.5-10_delay/','.mat')';
% delays = [0,2.5,5,7.5,10];
% files = fileSearch('E:/processed/ddt/0-8-16-32-60_delay/','.mat')';
% delays = [0;8;16;32;60];
% files = fileSearch('E:/processed/ddt/0-15-30-45-60_delay/','.mat')';
% delays = [0;15;30;45;60];
for fi = 1:size(files,1)
    load(files{fi})
    % check if file with no event (trial) data
    if numel(trials) == 1
        disp(files{fi})
        % Create new trial structure
        blockPercent = input('Block percent delay = ');
        auc = repmat(trials.auc,5,1);
        trapzAUC = repmat(trapz(delays,blockPercent),5,1);
        
        this = [blockPercent,auc,trapzAUC,delays];
        newTrials = cell2struct(num2cell(this)',{'blockPercentDelay','auc','trapzAUC','delay'});
        trials = newTrials;
    else
        thisTrial = struct2cell(trials);
        
        for ii = 1:size(trials,1)
            trials(ii).trapzAUC = trapz(delays,cell2mat(thisTrial(12,cell2mat(thisTrial(5,:)) == 1)));
        end
        for ii = 1:5
            blockInds = cell2mat(thisTrial(4,:)) == ii;
            for jj = logicFind(1,blockInds,'==')
                trials(jj).delay = delays(ii);
            end
        end
    end
    save(files{fi},'-append','trials')
end
%%
% [data,samps,files] = collateData('E:\processed\ddt\0-1-3-6-12_delay\',{'.mat'},{'pow','coh'},'avg','rel');
% catData = cat(1,data{1}{:});
% catSamps = samps{1};
% catFiles = files{1};
% [data,samps,files] = collateData('E:\processed\ddt\0-2.5-5-7.5-10_delay\',{'.mat'},{'pow','coh'},'avg','rel');
% catData = [catData;cat(1,data{1}{:})];
% catSamps = [catSamps;samps{1}];
% catFiles = [catFiles;files{1}];
% [data,samps,files] = collateData('E:\processed\ddt\0-4-8-16-32_delay\',{'.mat'},{'pow','coh'},'avg','rel');
% inds = cellfun(@numel,data{1})==216;
% catData = [catData;cat(1,data{1}{inds})];
% catSamps = [catSamps;samps{1}(inds)];
% catFiles = [catFiles;files{1}(inds)];
% [data,samps,files] = collateData('E:\processed\ddt\0-8-16-32-60_delay\',{'.mat'},{'pow','coh'},'avg','rel');
% catData = [catData;cat(1,data{1}{:})];
% catSamps = [catSamps;samps{1}];
% catFiles = [catFiles;files{1}];
% [data,samps,files] = collateData('E:\processed\ddt\0-15-30-45-60_delay\',{'.mat'},{'pow','coh'},'avg','rel');
% inds = cellfun(@numel,data{1})==216;
% catData = [catData;cat(1,data{1}{inds})];
% catSamps = [catSamps;samps{1}(inds)];
% catFiles = [catFiles;files{1}(inds)];
[data,samps,files] = collateData('F:/irdmRound2/processedContinuous3/',...
    {'IRDM11';'IRDM14';'IRDM15';'IRDM16';'IRDM18';'IRDM21';'IRDM22'},...
    {'pow','coh'},'avg','rel');
%% state - trapzAUC - all
load('E:\processed\ddt\catData.mat')
% get AUCs
for ii = 1:size(catFiles,1)
   load(catFiles{ii},'trials') 
   auc(ii) = trials(1).trapzAUC;
end
% model data
trainInds = randperm(size(auc,2),ceil(size(auc,2)*.8));
testInds = ismember(1:size(auc,2),trainInds);
trainX = catData(trainInds,:);
trainY = auc(trainInds);
testX = catData(testInds,:);
testY = auc(testInds);
cfg = lassoNetCfg({testX,testY'},[],'n','y','n',100,'1se',[]);
[~,lam,beta,fits,acc,hist] = lassoNet(trainX,trainY','poisson','deviance',1,10,1,cfg);

cfg = lassoNetCfg({testX,testY'},[],'y','y','n',100,'1se',[]);
[~,lamRand,betaRand,fitsRand,accRand,histRand] = lassoNet(trainX,trainY','poisson','deviance',1,10,1,cfg);

testY = testY(randperm(numel(testY)));
cfg = lassoNetCfg({testX,testY'},[],'n','y','n',100,'1se',[]);
[~,lamRandMan,betaRandMan,fitsRandMan,accRandMan,histRandMan] = lassoNet(trainX,trainY','poisson','deviance',1,10,1,cfg);
%% trait analysis
% y data from excel sheet Average AUC (0-32)
y = [1;1;1;1;0;0;1;1;0;0;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0;1;1;1;1;0;0;0;0;0;0];
% fedBase
[fedData,fedSamps,fedFiles] = collateData('E:/processed/fedBase/',{'.mat'},{'pow','coh'},'avg','rel');
fedData = cat(1,fedData{1}{:});
% depBase
[depData,depSamps,depFiles] = collateData('E:/processed/depBase/',{'.mat'},{'pow','coh'},'avg','rel');
depData = cat(1,depData{1}{:});
save('G:\GreenLab\data\irdmNew\fedDepBaseModel\depFedBaseModelData.mat','fedData','depData','y')
% run basic models: binomial, median-split y, internal cross-validated
% error
cfg = lassoNetCfg([],[],'n','y','n',100,'1se',[]);
% food deprived baseline model
[~,depLam,depBeta,depFits,depAcc,depHist] = lassoNet(depData,y,'binomial','class',1,5,1,cfg);
depA = 1-depLam{1}.allErr;
% fed baseline model
[~,fedLam,fedBeta,fedFits,fedAcc,fedHist] = lassoNet(fedData,y,'binomial','class',1,5,1,cfg);
fedA = 1-fedLam{1}.allErr;
% two animals each group left out
oneInd = logicFind(1,y,'==');
zeroInd = logicFind(0,y,'==');
oneCmbs = nchoosek(1:2:numel(oneInd),2);
zeroCmbs = nchoosek(1:2:numel(zeroInd),2);
c = 1;
for ii = 1:size(oneCmbs,1)
    for jj = 1:size(zeroCmbs,1)
        cmbs(c,:) = [oneInd(oneCmbs(ii,:)),oneInd(oneCmbs(ii,:)+1),zeroInd(zeroCmbs(jj,:)),zeroInd(zeroCmbs(jj,:)+1)];
        c = c+1;
    end
end
% single feature trait models
for ii = 1:100
    rng(ii)
    thisCmb = randi(size(cmbs,1),1);
    trainXdep = depData(~ismember(1:46,cmbs(thisCmb,:)),:);
    trainXfed = fedData(~ismember(1:46,cmbs(thisCmb,:)),:);
    trainY = y(~ismember(1:46,cmbs(thisCmb,:)));
    testXdep = depData(cmbs(thisCmb,:),:);
    testXfed = fedData(cmbs(thisCmb,:),:);
    testY = y(cmbs(thisCmb,:));
    for jj = 1:216
        fedMdl = fitglm(trainXfed(:,jj),trainY,'distribution','binomial','binomialsize',size(trainY,1));
        fedProb = predict(fedMdl,testXfed(:,jj));
        [~,~,~,fedA(ii,jj)] = perfcurve(testY,fedProb,1);
        fedSign(ii,jj) = table2array(fedMdl.Coefficients(2,1))/abs(table2array(fedMdl.Coefficients(2,1)));
        
        depMdl = fitglm(trainXdep(:,jj),trainY,'distribution','binomial','binomialsize',size(trainY,1));
        depProb = predict(depMdl,testXdep(:,jj));
        [~,~,~,depA(ii,jj)] = perfcurve(testY,depProb,1);
        depSign(ii,jj) = table2array(depMdl.Coefficients(2,1))/abs(table2array(depMdl.Coefficients(2,1)));        
    end
end
fedAM = mean(fedA,1).*(mean(fedSign,1)./abs(mean(fedSign,1)));
depAM = mean(depA,1).*(mean(depSign,1)./abs(mean(depSign,1)));
%%
for ii = 1:40
    load(['G:\GreenLab\data\irdmNew\fedDepBaseModel\all_leave2out\traitBinary_depFedBase_randAnimalDetector_',num2str(ii),'.mat'])
    fed(ii) = fedAcc{1}.acc;
    fedRand(ii) = fedAccRand{1}.acc;
    fedAD(ii) = fedAccAD{1}.acc;
    
    dep(ii) = depAcc{1}.acc;
    depRand(ii) = depAccRand{1}.acc;
    depAD(ii) = depAccAD{1}.acc;
end
%% plot trait models
load('G:\GreenLab\data\irdmNew\fedDepBaseModel\traitModels.mat')
for ii = 1:100
   load(['G:\GreenLab\data\irdmNew\fedDepBaseModel\animalDetectorNew\traitBinary_depFedBase_randAnimalDetector_',num2str(ii),'.mat'])
   depRand(ii,:) = 1-depLamRand{1}.allErr;
   fedRand(ii,:) = 1-fedLamRand{1}.allErr;
   depAD(ii,:) = 1-depLamAD{1}.allErr;
   fedAD(ii,:) = 1-fedLamAD{1}.allErr;
end
figure
this = histfit(mean(fedAD,2));
thisX = this(2).XData;
thisY = this(2).YData;
close(gcf)
doubleHist(fedA,mean(fedRand,2),'main','Fed Baseline vs. Animal Detector','xlab','accuracy')
hold on
plot(thisX,thisY./100,'--','color',[0.5 0.5 0.5],'lineWidth',2)

figure
this = histfit(mean(depAD,2));
thisX = this(2).XData;
thisY = this(2).YData;
close(gcf)
doubleHist(depA,mean(depRand,2),'main','Dep Baseline vs. Animal Detector','xlab','accuracy')
hold on
plot(thisX,thisY./100,'--','color',[0.5 0.5 0.5],'lineWidth',2)
%% state analysis
load('E:\processed\ddt\trapzAUCmodel_32delay.mat')
% get AUCs
for ii = 1:size(catFiles,1)
   load(catFiles{ii},'trials') 
   auc(ii) = trials(1).trapzAUC;
end
% all data
for ii = 1:100
    load(['G:\GreenLab\data\irdmNew\trapzAUCmodel\full\trapzAUC_',num2str(ii),'.mat'],'acc','accRand','accRandMan')
    a(ii,:) = acc{1}.acc;
    aR(ii,:) = accRand{1}.acc;
    aRM(ii,:) = accRandMan{1}.acc;
end
doubleHist(mean(a,2),mean(aRM,2))
% delay 32 - continuous
for ii = 1:100
    load(['G:\GreenLab\data\irdmNew\trapzAUCmodel\32DelayNew\trapzAUC32delay_',num2str(ii),'.mat'],'acc','accRand')
    a32(ii,:) = acc{1}.acc;
    pred = cvglmnetPredict(accRand{1}.mdl,histRand.cfg.naive.testX);
    aR32(ii,:) = accRand{1}.acc;
end
doubleHist(mean(a32,2),mean(aR32,2),'main','Default delay: Continuous trapz AUC','xlab','error')
% delay 32 - median split
for ii = 1:100
    load(['G:\GreenLab\data\irdmNew\trapzAUCmodel\32DelayBinary\trapzAUC32delay_medianSplit_',num2str(ii),'.mat'],'acc','accRand','accRandMan')
    [x(ii,:),y(ii,:),~,auc(ii)] = perfcurve(hist.cfg.naive.testY,acc{1}.pred,1,'TVals',0:0.001:1,'UseNearest',0);
    [xR(ii,:),yR(ii,:),~,aucR(ii)] = perfcurve(histRand.cfg.naive.testY,accRand{1}.pred,1,'TVals',0:0.001:1,'UseNearest',0);
%     auc32(ii,:) = acc{1}.acc;
%     aucR32(ii,:) = accRand{1}.acc;
%     aucRM32(ii,:) = accRandMan{1}.acc;
end
figure
hold on
plot(mean(x,1),mean(y,1),'k')
plot(mean(xR,1),mean(yR,1),'--k')
legend({['Real: \mu = ',num2str(round(mean(auc),2)),'\pm',...
    num2str(round(conf(auc,0.95),2))],...
    ['Real: \mu = ',num2str(round(mean(aucR),2)),'\pm',...
    num2str(round(conf(aucR,0.95),2))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('Trapz AUC median split: default delays')
% doubleHist(mean(auc32,2),mean(aucR32,2),'main','Default delay: Median split trapz AUC','xlab','auc')
%% single feature delay 32 median split
load('G:\GreenLab\data\irdmNew\trapzAUCmodel\trapzAUCmodel_32delay.mat')
aucBin = auc>=median(auc);
for ii = 1:100
    disp(ii)
    rng(ii)
    trainInds = randperm(size(aucBin,2),ceil(size(aucBin,2)*.8));
    testInds = ismember(1:size(aucBin,2),trainInds);
    trainX = catData(trainInds,:);
    trainY = aucBin(trainInds);
    testX = catData(testInds,:);
    testY = aucBin(testInds);
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial',...
            'binomialsize',numel(trainY));
        pred = predict(mdl,testX(:,jj));
        [~,~,~,a(ii,jj)] = perfcurve(testY,pred,1);
        sign(ii,jj) = table2array(mdl.Coefficients(2,1))/abs(table2array(mdl.Coefficients(2,1)));
    end
end
aM = mean(a,1).*(mean(sign,1)./abs(mean(sign,1)));
%% compare single features
load('G:\GreenLab\data\irdmNew\singleFeatures.mat')
c = distinguishable_colors(6);
figure
hold on
% scatter(fedAM(1:48),depAM(1:48),[],repmat(c,8,1),'o')
% scatter(fedAM(49:end),depAM(49:end),[],repmat(c,28,1),'s')
% scatter(fedAM(1:48),aM(1:48),[],repmat(c,8,1),'o')
% scatter(fedAM(49:end),aM(49:end),[],repmat(c,28,1),'s')
scatter(delayImmediateGenAM(1:48),aM(1:48),[],repmat(c,8,1),'o')
scatter(delayImmediateGenAM(49:end),aM(49:end),[],repmat(c,28,1),'s')
xlim([-1 1]); ylim([-1 1])
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
plot([-1 -0.5],[-0.5 -0.5],':k')
plot([1 0.5],[-0.5 -0.5],':k')
plot([-1 -0.5],[0.5 0.5],':k')
plot([1 0.5],[0.5 0.5],':k')

plot([-0.5 -0.5],[-1 -0.5],':k')
plot([-0.5 -0.5],[1 0.5],':k')
plot([0.5 0.5],[-1 -0.5],':k')
plot([0.5 0.5],[1 0.5],':k')
% xlabel('AUC Fed Baseline')
xlabel('Delay (1) v. Immediate (0)')
% ylabel('AUC Dep Baseline')
ylabel('AUC State 32 Delay')
%% find features that don't move quadrant
pos = aM > 0.5 & fedAM > 0.5 & delayImmediateGenAM > 0.5;
neg = aM < -0.5 & fedAM < -0.5 & delayImmediateGenAM < -0.5;
either = pos | neg;
powInd = logicFind(1,either(1:48),'==');
cohInd = logicFind(1,either(49:end),'==')+48;
figure
hold on
c = repmat(distinguishable_colors(6),36,1);
scatter(delayImmediateGenAM(powInd),aM(powInd),[],c(powInd,:),'o')
scatter(delayImmediateGenAM(cohInd),aM(cohInd),[],c(cohInd,:),'s')
feats = [fedAM(either);aM(either);delayImmediateGenAM(either)]';
%%
inds = logicFind(1,pos,'==');
n = 138;
delayImmediateGenAM(n)
fedAM(n)
aM(n)
%%
for ii = 1:100
    load(['G:\GreenLab\data\irdmNew\delayImmediate\genModel\delayImmediateGen',num2str(ii),'.mat'],'sign','a')
    delayImmediateGenA (ii,:) = a(:,end);
    delayImmediateGenSign(ii,:) = sign(:,end);
end
delayImmediateGenSignM = mean(delayImmediateGenSign,1);
delayImmediateGenAM = mean(delayImmediateGenA,1).*(delayImmediateGenSignM./abs(delayImmediateGenSignM));
%% immediate vs delay levers
for ii = 1:100
    load(['G:\GreenLab\data\irdmNew\delayImmediate\genModel\delayImmediateGen',num2str(ii),'.mat'],'accArray','accArrayP','hist','histP','sign','a')
    [idX(ii,:),idY(ii,:),~,idAuc(ii)] = perfcurve(hist.cfg.naive.testY,accArray{1}.pred,1);
    idSign(ii,:) = sign;
    [idXP(ii,:),idYP(ii,:),~,idAucP(ii)] = perfcurve(histP.cfg.naive.testY,accArrayP{1}.pred,1);
end
figure
hold on
plot(mean(idX,1),mean(idY),'-k')
plot(mean(idXP,1),mean(idYP,1),'--k')
legend({['Real: \mu = ',num2str(round(mean(idAuc),2)),'\pm',...
    num2str(round(conf(idAuc,0.95),3))],...
    ['Permuted: \mu = ',num2str(round(mean(idAucP),2)),'\pm',...
    num2str(round(conf(idAucP,0.95),3))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('Immediate vs. Delay Lever Choice')
%% individual animal immediate vs delay
load('delayImmediateModelData.mat')
u = unqiue(allData(1,:));
for ii = 1:numel(u)
    nD(ii) = sum(cell2mat(allData(7,logicFind(u{ii},allData(1,:),'==')))==1);
    nI(ii) = sum(cell2mat(allData(7,logicFind(u{ii},allData(1,:),'==')))==0);
end
figure
bar([nD;nI]','stacked')
%%
animals = {'IRDM11','IRDM14','IRDM15','IRDM16','IRDM18','IRDM21','IRDM22','IRDM27','IRDM28','IRDM31','IRDM32','IRDM39'};
% collect data
for ii = 1:100
    for jj = 1:numel(animals)
        load(['G:\GreenLab\data\irdmNew\delayImmediate\indModel\delayImmediate',animals{jj},'_',num2str(ii),'.mat'])
        a(ii,jj) = accArray{1}.acc;
        aP(ii,jj) = accArrayP{1}.acc;
        mdls(ii,jj) = accArray{1}.mdl;
        testX{ii,jj} = hist.cfg.naive.testX;
        testY{ii,jj} = hist.cfg.naive.testY;
        [indX{ii,jj},indY{ii,jj},~,indA(ii,jj)] = perfcurve(testY{ii,jj},accArray{1}.pred,1);
    end
end
% ii = train
for ii = 1:numel(animals)
    % jj = test
    for jj = 1:numel(animals)
        % k = model iteration
        for k = 1:100
            pred = cvglmnetPredict(mdls(k,ii),testX{k,jj},'lambda_1se','response');
            [waffleX(ii,jj,k,:),waffleY(ii,jj,k,:),~,waffleA(ii,jj,k)] = perfcurve(testY{k,jj},pred,1);
        end
    end
end
%%
for ii = 1:12
    
end
%%
c = distinguishable_colors(numel(animals));
figure
hold on
for ii = 1:numel(animals)
    plot(squeeze(mean(indX(:,ii,:),1)),squeeze(mean(indY(:,ii,:),1)),'color',c(ii,:))
end
legend(animals)

figure
pcolor(padarray(mean(waffleA,3)',[1 1],'post'))
set(gca,'xtick',1.5:12.5,'xticklabel',animals,'ytick',1.5:12.5,'yticklabel',animals)
xtickangle(45)
colormap viridis
xlabel('Train')
ylabel('Test')
%% compare continuous to discrete
load('G:\GreenLab\data\irdmNew\processedContinuous\IRDM11_DDT_2019-03-15.mat')
contPow = psdTrls;
contCoh = coh;
load('E:\processed\IRDM11_DDT_2019-03-15_all.mat')
%%
figure
subplot(2,2,1)
set(gca,'ylim',[-70 -20])
% plot band lines
for ii = 1:6
    rectangle('Position',[hist.bands{ii,2}(1),-70,diff(hist.bands{ii,2}),diff(get(gca,'ylim'))],'FaceColor',[0.9 0.9 0.9])
end
hold on
plot(psdTrls{1}.f,mean(psdTrls{1}.Pow(2,:,:),3))
plot(contPow{1}.f,squeeze(mean(contPow{1}.Pow(2,:,1,:),4,'omitnan')),'-r')
legend({'Discrete','Continuous'})


subplot(2,2,2)
imagesc(psdTrls{1}.avgRelPow)
colormap viridis
title('Discrete Normalized Power')
set(gca,'yticklabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})

subplot(2,2,3)
imagesc(squeeze(mean(contPow{1}.avgRelPow,4,'omitnan')))
colormap viridis
title('Continuous Normalized Power')
set(gca,'yticklabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})

subplot(2,2,4)
imagesc(squeeze(mean(contPow{1}.avgRelPow,4,'omitnan'))-psdTrls{1}.avgRelPow)
colormap viridis
title('Continuous - Discrete: Normalized Power')
set(gca,'yticklabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})

figure
subplot(2,2,1)
hold on
set(gca,'ylim',[0 1])
% plot band lines
for ii = 1:6
    rectangle('Position',[hist.bands{ii,2}(1),0,diff(hist.bands{ii,2}),diff(get(gca,'ylim'))],'FaceColor',[0.9 0.9 0.9])
end
plot(coh{1}.f,squeeze(mean(coh{1}.Cxy(1,:,:),3)))
plot(contCoh{1}.f,squeeze(mean(contCoh{1}.Cxy(1,:,1,:),4,'omitnan')))
legend({'Discrete','Continuous'})

subplot(2,2,2)
imagesc(squeeze(mean(coh{1}.normBandCoh,3,'omitnan')))
colormap viridis
title('Discrete Normalized Coherence')
set(gca,'xtick',1:6,'xticklabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})

subplot(2,2,3)
imagesc(squeeze(squeeze(mean(contCoh{1}.normBandCoh,4,'omitnan'))))
colormap viridis
title('Continuous Normalized Coherence')
set(gca,'xtick',1:6,'xticklabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})

subplot(2,2,4)
imagesc(squeeze(squeeze(mean(contCoh{1}.normBandCoh,4,'omitnan')))-squeeze(mean(coh{1}.normBandCoh,3,'omitnan')))
colormap viridis
title('Continuous - Discrete: Normalized Coherence')
set(gca,'xtick',1:6,'xticklabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'})