%% Collate raw data: LSD vs. Saline; acute effects from that days' baseline
[data,samp,files] = collateData(['G:\Greenlab\data\lsd\processed\',...
    'imagAcute\'],{'lsd';'sal'},{'pow','coh'},'trl','');
% For each animal, average that day's baseline and use to normalize post
allData = [data{1,1};data{1,2}];
mBase = cellfun(@(x) mean(x,1),allData(:,1),'UniformOutput',0);
for ii = 1:numel(mBase)
    dif{ii} = allData{ii,2}-repmat(mBase{ii},size(allData{ii,2},1),1);
end
allPostLSD = dif(1:numel(files{1,1}));
allPostSal = dif(numel(files{1,1})+1:end);
minSamp = min(cellfun(@(x) size(x,1),dif));
% save('G:\GreenLab\data\lsd\normAcuteImag.mat','allPostLSD','allPostSal','minSamp',...
%     'data','files')
% Only use data from the last 30 minutes
c = 1;
for ii = 1:2
    for jj = 1:size(data{1,ii},1)
        allData30(c,:) = [data{1,ii}(jj,1),...
            data{1,ii}{jj,2}(nearest_idx3(samp{1,ii}{jj,2}(end,1)...
            -30*60*400,samp{1,ii}{jj,2}(:,1)):end,:)];
        c = c+1;
    end
end
% mBase doesn't change
for ii = 1:numel(mBase)
    dif30{ii} = allData30{ii,2}-repmat(mBase{ii},size(allData30{ii,2},1),1);
end
allPostLSD30 = dif30(1:numel(files{1,1}));
allPostSal30 = dif30(numel(files{1,1})+1:end);
minSamp30 = min(cellfun(@(x) size(x,1),dif30));
% save('G:\GreenLab\data\lsd\normAcuteImag30.mat','allPostLSD30',...
%     'allPostSal30','minSamp30','data','files','allData30')
%% power changes in LSD relative to SAL
salFiles = fileSearch('G:\GreenLab\data\lsd\processed\imagAcute\','Sal');
[salPow,salCoh] = deal(cell(numel(salFiles),2));
for ii = 1:numel(salFiles)
    load(salFiles{ii},'psdTrls','coh')
    salPow{ii,1} = psdTrls{1}.Pow;
    salPow{ii,2} = psdTrls{2}.Pow;
    
    [b,c,t] = size(psdTrls{1}.bandPow);
    salBand{ii,1} = reshape(psdTrls{1}.bandPow,b*c,t)';  
    [b,c,t] = size(psdTrls{2}.bandPow);
    salBand{ii,2} = reshape(psdTrls{2}.bandPow,b*c,t)';
    
    salCoh{ii,1} = coh{1}.Cxy;
    salCoh{ii,2} = coh{2}.Cxy;
end
lsdFiles = fileSearch('G:\GreenLab\data\lsd\processed\imagAcute\','LSD');
[lsdPow,lsdCoh] = deal(cell(numel(lsdFiles),2));
for ii = 1:numel(lsdFiles)
    load(lsdFiles{ii},'psdTrls','coh')
    lsdPow{ii,1} = psdTrls{1}.Pow;
    lsdPow{ii,2} = psdTrls{2}.Pow;
    
    [b,c,t] = size(psdTrls{1}.bandPow);
    lsdBand{ii,1} = reshape(psdTrls{1}.bandPow,b*c,t)';  
    [b,c,t] = size(psdTrls{2}.bandPow);
    lsdBand{ii,2} = reshape(psdTrls{2}.bandPow,b*c,t)';
    
    lsdCoh{ii,1} = coh{1}.Cxy;
    lsdCoh{ii,2} = coh{2}.Cxy;
end
%% double check band power
bInd = [1 4;5 10;11 14;15 30;45 65;70 90];
for ii = 1:size(salPow,1)
    for jj = 1:size(salPow,2)
        these = [];
        for k = 1:size(bInd,1)
            data = squeeze(trapz(salPow{ii,jj}(:,bInd(k,1):bInd(k,2),:),2));
            these(k,:,:) = data;
        end
        [b,c,t] = size(these);
        salBandChk{ii,jj} = reshape(these,b*c,t)';
        if any(salBandChk{ii,jj}~=salBand{ii,jj})
           disp('fuck') 
        end
    end
end
for ii = 1:size(lsdPow,1)
    for jj = 1:size(lsdPow,2)
        these = [];
        for k = 1:size(bInd,1)
            data = squeeze(trapz(lsdPow{ii,jj}(:,bInd(k,1):bInd(k,2),:),2));
            these(k,:,:) = data;
        end
        [b,c,t] = size(these);
        lsdBandChk{ii,jj} = reshape(these,b*c,t)';
        if any(lsdBandChk{ii,jj}~=lsdBand{ii,jj})
           disp('fuck') 
        end
    end
end
%% plot pre post lsd and sal - power
allLSDpre = cat(3,lsdPow{:,1});
allLSDpost = cat(3,lsdPow{:,2});
mLSDpre = mean(allLSDpre,3);
mLSDpost = mean(allLSDpost,3);
lsdDiff = mLSDpost-mLSDpre;
lsdDiffS = (abs(std(allLSDpost,[],3)+std(allLSDpre,[],3)-2*...
    std(cat(3,allLSDpre,allLSDpost),[],3))).^0.5;

allSALpre = cat(3,salPow{:,1});
allSALpost = cat(3,salPow{:,2});
mSALpre = mean(allSALpre,3);
mSALpost = mean(allSALpost,3);
salDiff = mSALpost-mSALpre;
salDiffS = (abs(std(allSALpost,[],3)+std(allSALpre,[],3)-2*...
    std(cat(3,allSALpre,allSALpost),[],3))).^0.5;
x = (abs(std(cat(3,allSALpost,allLSDpost),[],3)+std(cat(3,allSALpost,allLSDpost),[],3)-2*...
    std(cat(3,allLSDpre,allLSDpost,allSALpre,allSALpost),[],3))).^0.5;
lsdSalDiff = lsdDiff-salDiff;
lsdSalDiffS = (abs(lsdDiffS+salDiffS-2*(x))).^0.5;
sites = {'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'};
figure
for ii = 1:8
    subplot(4,4,ii)
    hold on
    shadedErrorBar(1:100,lsdDiff(ii,:),lsdDiffS(ii,:),'r',1)
    shadedErrorBar(1:100,salDiff(ii,:),salDiffS(ii,:),'b',1)
%     shadedErrorBar(1:100,lsdSalDiff(ii,:),lsdSalDiffS(ii,:),'k',1)
    plot([1 1],[-6.5 2],':k')
    plot([4 4],[-6.5 2],':k')
    plot([5 5],[-6.5 2],':k')
    plot([10 10],[-6.5 2],':k')
    plot([11 11],[-6.5 2],':k')
    plot([14 14],[-6.5 2],':k')
    plot([15 15],[-6.5 2],':k')
    plot([30 30],[-6.5 2],':k')
    plot([45 45],[-6.5 2],':k')
    plot([65 65],[-6.5 2],':k')
    plot([70 70],[-6.5 2],':k')
    plot([90 90],[-6.5 2],':k')
    plot([0 100],[0 0],'--k')
    set(gca,'xtick',0:10:100)
    xlabel('frequency (Hz)')
    ylabel('a.u.')
    ylim([-6.6 1])
    title(sites{ii})
end
for ii = 1:8
   subplot(4,4,8+ii)
   hold on
    shadedErrorBar(1:100,lsdSalDiff(ii,:),lsdSalDiffS(ii,:),'k',1)
     set(gca,'xtick',0:10:100)
     ylim([-5 1])
     plot([0 100],[0 0],'--k')
end
%% coh - coherogram differences
allLSDpre = cat(3,lsdCoh{:,1});
allLSDpost = cat(3,lsdCoh{:,2});

% allLSDpre = cellfun(@(x) mean(x,3),lsdCoh(:,1),'UniformOutput',0);
% allLSDpre = cat(3,allLSDpre{:});
% 
% allLSDpost = cellfun(@(x) mean(x,3),lsdCoh(:,2),'UniformOutput',0);
% allLSDpost = cat(3,allLSDpost{:});
mLSDpre = mean(allLSDpre,3);
mLSDpost = mean(allLSDpost,3);
lsdDiff = mLSDpost-mLSDpre;
% Linear propogation of error
lsdDiffS = (abs(std(allLSDpost,[],3)+std(allLSDpre,[],3)-2*...
    std(cat(3,allLSDpre,allLSDpost),[],3))).^0.5;

allSALpre = cat(3,salCoh{:,1});
allSALpost = cat(3,salCoh{:,2});
% allSALpre = cellfun(@(x) mean(x,3),salCoh(:,1),'UniformOutput',0);
% allSALpre = cat(3,allSALpre{:});
% 
% allSALpost = cellfun(@(x) mean(x,3),salCoh(:,2),'UniformOutput',0);
% allSALpost = cat(3,allSALpost{:});
mSALpre = mean(allSALpre,3);
mSALpost = mean(allSALpost,3);
salDiff = mSALpost-mSALpre;
% Linear propogation of error
salDiffS = (abs(std(allSALpost,[],3)+std(allSALpre,[],3)-2*...
    std(cat(3,allSALpre,allSALpost),[],3))).^0.5;

x = (abs(std(cat(3,allSALpost,allLSDpost),[],3)+std(cat(3,allSALpost,allLSDpost),[],3)-2*...
    std(cat(3,allLSDpre,allLSDpost,allSALpre,allSALpost),[],3))).^0.5;
lsdSalDiff = lsdDiff-salDiff;
lsdSalDiffS = (abs(lsdDiffS+salDiffS-2*(x))).^0.5;
sites = {'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'};
cmbs = nchoosek(1:8,2);
cmbs = cmbs([3,6,7,21,22,28],:);
figure
for ii = 1:size(cmbs,1)
    subplot(3,2,ii)
    hold on
    shadedErrorBar(1:100,decimate(lsdDiff(ii,:),5),...
        decimate(lsdDiffS(ii,:),5),'r',1)
    shadedErrorBar(1:100,decimate(salDiff(ii,:),5),...
        decimate(salDiffS(ii,:),5),'b',1)
    shadedErrorBar(1:100,decimate(lsdSalDiff(ii,:),5),...
        decimate(lsdSalDiffS(ii,:),5),'k',1)
    plot([1 1],[-1 1],':k')
    plot([4 4],[-1 1],':k')
    plot([5 5],[-1 1],':k')
    plot([10 10],[-1 1],':k')
    plot([11 11],[-1 1],':k')
    plot([14 14],[-1 1],':k')
    plot([15 15],[-1 1],':k')
    plot([30 30],[-1 1],':k')
    plot([45 45],[-1 1],':k')
    plot([65 65],[-1 1],':k')
    plot([70 70],[-1 1],':k')
    plot([90 90],[-1 1],':k')
    plot([0 100],[0 0],'--k')
    set(gca,'xtick',0:10:100)
    ylim([-0.3 0.3])
    xlabel('frequency (Hz)')
    ylabel('a.u.')
    title([sites{cmbs(ii,1)},'-',sites{cmbs(ii,2)}])
end
%% coh - ecdfs
[data,samp,files] = collateData(['G:\Greenlab\data\lsd\processed\',...
    'imagAcute\'],{'lsd';'sal'},{'pow','coh'},'trl','');

subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
feat = feat(subInds);
for ii = 1:size(data{1},1)
    lsdDiff{ii} = mean(data{1}{ii,1},1)-data{1}{ii,2};
end
allLSDdiff = cat(1,lsdDiff{:});
allLSDdiff = allLSDdiff(:,subInds);
for ii = 1:size(data{2},1)
    salDiff{ii} = mean(data{2}{ii,1},1)-data{2}{ii,2};
end
allSALdiff = cat(1,salDiff{:});
allSALdiff = allSALdiff(:,subInds);

figure
for ii = 1:36
    subplot(6,6,ii)
%     plotSpread({allSALdiff(:,24+ii),allLSDdiff(:,24+ii)},...
%         'distributionColor',{'b','r'})
    ecdf(allSALdiff(:,24+ii))
    hold on
    ecdf(allLSDdiff(:,24+ii))
    xlim([-0.5 0.5])
    title(feat{24+ii})
end
%% LSD vs. Base and SAL vs. Base - Log LOO
load('G:\GreenLab\data\lsd\normAcuteImag30.mat')
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
minSamp30 = min(min(cellfun(@(x) size(x,1),allData30)));
% LSD vs. base - LOO
for ii = 1:6
    testX = cat(1,allData30{ii,1}(randperm(size(allData30{ii,1},1),...
        minSamp30),subInds),...
        allData30{ii,2}(randperm(size(allData30{ii,2},1),minSamp30),...
        subInds));
    testY = cat(1,zeros(minSamp30,1),ones(minSamp30,1));
    otherInd = logicFind(1,~ismember(1:6,ii),'==');
    trainX = [];
    trainY = [];
    for jj = otherInd
        trainX = [trainX;cat(1,...
            allData30{jj,1}(randperm(size(allData30{jj,1},1),...
            minSamp30),subInds),...
            allData30{jj,2}(randperm(size(allData30{jj,2},1),...
            minSamp30),subInds))];
        trainY = [trainY;cat(1,zeros(minSamp30,1),ones(minSamp30,1))];
    end
    mdl = fitglm(zscore(trainX),trainY,'Distribution','binomial');
    prob = predict(mdl,zscore(testX));
    [~,~,~,aLSD(ii)] = perfcurve(testY,prob,1);
end
% LSD vs. base - 80:20
for ii = 1:100
    for jj = 1:6
        thisData{jj,1} = allData30{jj,1}(randperm(size(allData30{jj,1},1),minSamp30),subInds);
        thisData{jj,2} = allData30{jj,2}(randperm(size(allData30{jj,2},1),minSamp30),subInds);
    end
    thisX = cat(1,thisData{:,1},thisData{:,2});
    thisY = cat(1,zeros(minSamp30*6,1),ones(minSamp30*6,1));
    [trainX,trainY,testX,testY] = trainTest(thisX,thisY,0.2);
    mdl = fitglm(zscore(trainX),trainY,'Distribution','binomial');
    prob = predict(mdl,zscore(testX));
    [~,~,~,aLSD80(ii)] = perfcurve(testY,prob,1);
end
% permuted LSD vs. base
for ii = 1:6
    testX = cat(1,allData30{ii,1}(randperm(size(allData30{ii,1},1),...
        minSamp30),subInds),...
        allData30{ii,2}(randperm(size(allData30{ii,2},1),minSamp30),...
        subInds));
    testY = cat(1,zeros(minSamp30,1),ones(minSamp30,1));
    otherInd = logicFind(1,~ismember(1:6,ii),'==');
    otherInd1 = nchoosek(otherInd,3);
    for jj = 1:size(otherInd1,1)
        otherInd2(jj,:) = otherInd(~ismember(otherInd,otherInd1(jj,:)));
    end
    for jj = 1:size(otherInd1,1)
        trainX = [];
        trainY = [];
        for k = 1:3
            trainX = [trainX;allData30{otherInd1(jj,k),2}(randperm(size(allData30{otherInd1(jj,k),2},1),minSamp30),subInds);allData30{otherInd1(jj,k),1}(randperm(size(allData30{otherInd1(jj,k),1},1),minSamp30),subInds)];
            trainY = [trainY;ones(minSamp30,1);zeros(minSamp30,1)];
        end
        for k = 1:2
            trainX = [trainX;allData30{otherInd1(jj,k),2}(randperm(size(allData30{otherInd1(jj,k),2},1),minSamp30),subInds);allData30{otherInd1(jj,k),1}(randperm(size(allData30{otherInd1(jj,k),1},1),minSamp30),subInds)];
            trainY = [trainY;zeros(minSamp30,1);ones(minSamp30,1)];
        end
        mdl = fitglm(zscore(trainX),trainY,'Distribution','binomial');
        prob = predict(mdl,zscore(testX));
        [~,~,~,aLSDp(ii,jj,1)] = perfcurve(testY,prob,1);
        trainX = [];
        trainY = [];
        for k = 1:3
            trainX = [trainX;allData30{otherInd1(jj,k),2}(randperm(size(allData30{otherInd1(jj,k),2},1),minSamp30),subInds);allData30{otherInd1(jj,k),1}(randperm(size(allData30{otherInd1(jj,k),1},1),minSamp30),subInds)];
            trainY = [trainY;zeros(minSamp30,1);ones(minSamp30,1)];
        end
        for k = 1:2
            trainX = [trainX;allData30{otherInd1(jj,k),2}(randperm(size(allData30{otherInd1(jj,k),2},1),minSamp30),subInds);allData30{otherInd1(jj,k),1}(randperm(size(allData30{otherInd1(jj,k),1},1),minSamp30),subInds)];
            trainY = [trainY;ones(minSamp30,1);zeros(minSamp30,1)];
        end
        mdl = fitglm(zscore(trainX),trainY,'Distribution','binomial');
        prob = predict(mdl,zscore(testX));
        [~,~,~,aLSDp(ii,jj,2)] = perfcurve(testY,prob,1);
    end
end
% sal vs. base
c = 1;
for ii = 7:11
    testX = cat(1,allData30{ii,1}(randperm(size(allData30{ii,1},1),...
        minSamp30),subInds),...
        allData30{ii,2}(randperm(size(allData30{ii,2},1),minSamp30),subInds));
    testY = cat(1,zeros(minSamp30,1),ones(minSamp30,1));
    otherInd = logicFind(1,~ismember(7:11,ii),'==')+6;
    trainX = [];
    trainY = [];
    for jj = otherInd
        trainX = [trainX;cat(1,...
            allData30{jj,1}(randperm(size(allData30{jj,1},1),...
            minSamp30),subInds),...
            allData30{jj,2}(randperm(size(allData30{jj,2},1),...
            minSamp30),subInds))];
        trainY = [trainY;cat(1,zeros(minSamp30,1),ones(minSamp30,1))];
    end
    mdl = fitglm(zscore(trainX),trainY,'Distribution','binomial');
    prob = predict(mdl,zscore(testX));
    [~,~,~,aSAL(c)] = perfcurve(testY,prob,1);
    c = c+1;
end
% permuted sal vs. base
for ii = 7:11
    testX = cat(1,allData30{ii,1}(randperm(size(allData30{ii,1},1),...
        minSamp30),subInds),...
        allData30{ii,2}(randperm(size(allData30{ii,2},1),minSamp30),...
        subInds));
    testY = cat(1,zeros(minSamp30,1),ones(minSamp30,1));
    otherInd = logicFind(1,~ismember(7:11,ii),'==')+6;
    otherInd1 = nchoosek(otherInd,2);
    for jj = 1:size(otherInd1,1)
        otherInd2(jj,:) = otherInd(~ismember(otherInd,otherInd1(jj,:)));
    end
    for jj = 1:size(otherInd1,1)
        trainX = [];
        trainY = [];
        for k = 1:2
            trainX = [trainX;allData30{otherInd1(jj,k),2}(randperm(size(allData30{otherInd1(jj,k),2},1),minSamp30),subInds);allData30{otherInd1(jj,k),1}(randperm(size(allData30{otherInd1(jj,k),1},1),minSamp30),subInds)];
            trainY = [trainY;zeros(minSamp30,1);ones(minSamp30,1)];
        end
        for k = 1:2
            trainX = [trainX;allData30{otherInd1(jj,k),2}(randperm(size(allData30{otherInd1(jj,k),2},1),minSamp30),subInds);allData30{otherInd1(jj,k),1}(randperm(size(allData30{otherInd1(jj,k),1},1),minSamp30),subInds)];
            trainY = [trainY;ones(minSamp30,1);zeros(minSamp30,1)];
        end
        mdl = fitglm(zscore(trainX),trainY,'Distribution','binomial');
        prob = predict(mdl,zscore(testX));
        [~,~,~,aSALp(ii,jj,1)] = perfcurve(testY,prob,1);
        trainX = [];
        trainY = [];
        for k = 1:2
            trainX = [trainX;allData30{otherInd1(jj,k),2}(randperm(size(allData30{otherInd1(jj,k),2},1),minSamp30),subInds);allData30{otherInd1(jj,k),1}(randperm(size(allData30{otherInd1(jj,k),1},1),minSamp30),subInds)];
            trainY = [trainY;ones(minSamp30,1);zeros(minSamp30,1)];
        end
        for k = 1:2
            trainX = [trainX;allData30{otherInd1(jj,k),2}(randperm(size(allData30{otherInd1(jj,k),2},1),minSamp30),subInds);allData30{otherInd1(jj,k),1}(randperm(size(allData30{otherInd1(jj,k),1},1),minSamp30),subInds)];
            trainY = [trainY;zeros(minSamp30,1);ones(minSamp30,1)];
        end
        mdl = fitglm(zscore(trainX),trainY,'Distribution','binomial');
        prob = predict(mdl,zscore(testX));
        [~,~,~,aSALp(ii,jj,2)] = perfcurve(testY,prob,1);
    end
end
(sum(aLSDp>mean(aLSD),[1,2,3])+1)/(numel(aLSDp)+1)
% save('G:\GreenLab\data\lsd\LSDvBase_SALvBase_acute.mat','aLSD','aSAL','aLSDp','aSALp')
%% Figures
figure 
subplot(3,1,1)
hold on
[f,xi,bw] = ksdensity(aLSD); 
fill(xi,f*bw,'w')
plotSpread(aLSD','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aLSD) mean(aLSD)],[0 0.25],'-')
yticks(0:0.05:0.25)
xlim([0 1])
ylim([0 0.25])
subplot(3,1,2)
hold on
[f,xi,bw] = ksdensity(aSAL); 
fill(xi,f*bw,'w')
plotSpread(aSAL','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aSAL) mean(aSAL)],[0 0.35],'-')
yticks(0:0.05:0.35)
xlim([0 1])
ylim([0 0.35])
subplot(3,1,3)
hold on
[f,xi,bw] = ksdensity(reshape(aLSDp,1,120)); 
fill(xi,f*bw,'w')
plotSpread(reshape(aLSDp,1,120)','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(reshape(aLSDp,1,120)) mean(reshape(aLSDp,1,120))],[0 0.25],'-')
yticks(0:0.05:0.25)
xlim([0 1])
ylim([0 0.25])
%% LSD v SAL - Log LOO - using last 30 minutes of SAL and LSD data
load('G:\GreenLab\data\lsd\normAcuteImag30.mat')
x = []; y = []; a = [];
xP = []; yP = []; aP = [];
xS = []; yS = []; aS = [];
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
count = 1;
c = 1;
lsdN = size(allPostLSD30,2);
salN = size(allPostSal30,2);
for ii = 1:lsdN
    for jj = 1:salN
        inds(c,:) = [ii,jj];
        c = c+1;
    end
end
for n = 1:30
    disp(n)
        % LSD data
        trainLSD = logicFind(1,~ismember(1:lsdN,inds(n,1)),'==');
        testLSD = logicFind(0,~ismember(1:lsdN,inds(n,1)),'==');

        thisLSD24 = [];
        for ii = trainLSD
            rng(ii*count)
            thisLSD24 = [thisLSD24;allPostLSD30{ii}(randperm(size(allPostLSD30{ii},1),...
                minSamp30),:)];
        end
        lsdTest = allPostLSD30{testLSD}(randperm(size(allPostLSD30{testLSD},1),...
            minSamp30),:);
        % Saline data
        trainSal = logicFind(1,~ismember(1:salN,inds(n,2)),'==');
        testSal = logicFind(0,~ismember(1:salN,inds(n,2)),'==');
        thisSal = [];
        for ii = 1:numel(allPostSal30)
            rng((ii+numel(allPostLSD30))*count)
            thisSal = [thisSal;allPostSal30{ii}(randperm(size(allPostSal30{ii},1),...
                minSamp30),:)];
        end
        salTest = allPostSal30{testSal}(randperm(size(allPostSal30{testSal},1),...
            minSamp30),:);
        % Combine and build models
        trainX = [thisLSD24;thisSal];
        trainY = [ones(size(thisLSD24,1),1);zeros(size(thisSal,1),1)];
        testX = [lsdTest;salTest];
        testY = [ones(minSamp30,1);zeros(minSamp30,1)];

        mdl = fitglm(zscore(trainX(:,subInds)),trainY,'distribution','binomial','binomialSize',numel(trainY));
        prob = predict(mdl,zscore(testX(:,subInds)));
        [x(count,:),y(count,:),~,a(count)] = perfcurve(testY,prob,1);
        beta(:,:,count) = table2array(mdl.Coefficients);
        % Single feature models
        for jj = 1:60
            theseTrain = trainX(:,subInds);
            theseTest = testX(:,subInds);
            mdl = fitglm(zscore(theseTrain(:,jj)),trainY,'distribution','binomial','binomialSize',numel(trainY));
            acuteBeta(jj,count,:) = table2array(mdl.Coefficients(2,:));
            prob = predict(mdl,zscore(theseTest(:,jj)));
            [xS(count,jj,:),yS(count,jj,:),~,aS(count,jj)] = perfcurve(testY,prob,1);
            [xSP(count,jj,:),ySP(count,jj,:),~,aSP(count,jj)] = perfcurve(testY(randperm(numel(testY),numel(testY))),prob,1);
        end
        count = count+1;
end
% Permuted
lsdInds = nchoosek(1:6,3);
salInds = nchoosek(1:5,3);
c = 1;
for ii = 1:size(lsdInds,1)
    for jj = 1:size(salInds,1)
        permInds(c,:) = cat(2,lsdInds(ii,:),salInds(jj,:));
        c = c+1;
    end
end
thisLSD = cell(1,numel(allPostLSD30));
thisSal = cell(1,numel(allPostSal30));
for ii = 1:size(permInds,1)
    disp(ii)
    lsdInd = permInds(ii,randperm(3,3));
    salInd = permInds(ii,randperm(3,3)+3);
    lsdOtherInd = logicFind(1,~ismember(1:6,lsdInd),'==');
    salOtherInd = logicFind(1,~ismember(1:5,salInd),'==');
    for jj = 1:numel(allPostLSD30)
        thisLSD{jj} = allPostLSD30{jj}(randperm(size(allPostLSD30{jj},1),...
            minSamp30),:);
    end
    for jj = 1:numel(allPostSal30)
        thisSal{jj} = allPostSal30{jj}(randperm(size(allPostSal30{jj},1),...
            minSamp30),:);
    end
    testX = cat(1,thisLSD{lsdInd(1)},thisSal{salInd(1)});
    testY = cat(1,ones(minSamp30,1),zeros(minSamp30,1)); 
    trainX = cat(1,thisLSD{lsdInd(2:3)},thisSal{salInd(2:3)},thisLSD{lsdOtherInd},thisSal{salOtherInd});
    trainY = cat(1,ones(minSamp30*4,1),zeros(minSamp30*5,1));
    mdl = fitglm(zscore(trainX(:,subInds)),trainY,'Distribution','binomial');
    prob = predict(mdl,zscore(testX(:,subInds)));
    [xP(ii,:),yP(ii,:),~,aP(ii)] = perfcurve(testY,prob,1);
    for jj = 1:numel(subInds)
        mdl = fitglm(zscore(trainX(:,subInds(jj))),trainY,'Distribution','binomial');
        prob = predict(mdl,testX(:,subInds(jj)));
        [~,~,~,aSP(ii,jj)] = perfcurve(testY,prob,1);
    end
end
%%
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
feat = feat(subInds);
%
figure
plot(mean(acuteBeta(:,:,1),2),mean(aS,1),'.k')
xlabel('mean beta')
ylabel('mean auc')
title('single feature')
%
figure
subplot(1,2,1)
plot(mean(acuteBeta(:,:,1),2),mean(acuteBeta(:,:,4),2),'.k')
xlabel('mean beta'); ylabel('mean p')
title('single feature')
subplot(1,2,2)
plot(mean(beta(2:61,1,:),3),mean(beta(2:61,4,:),3),'.k')
xlabel('mean beta'); ylabel('mean p')
title('full')
%
figure
plot(mean(acuteBeta(:,:,1),2),mean(beta(2:61,1,:),3),'.k')
xlabel('single feature beta'); ylabel('full feature beta')
% save('acuteLSDvSalineLogLOOImag30.mat','a','aP','aS','aSP','acuteBeta',...
%     'beta','x','xP','xS','xSP','y','yP','yS','yP')
%% Plot Log LOO
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
figure 
subplot(2,1,1)
hold on
[f,xi,bw] = ksdensity(a); 
fill(xi,f*bw,'w')
plotSpread(a','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(a) mean(a)],[0 0.25],'-')
yticks(0:0.05:0.25)
xlim([0 1])
ylim([0 0.25])
subplot(2,1,2)
hold on
[f,xi,bw] = ksdensity(aP); 
fill(xi,f*bw,'w')
plotSpread(aP','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aP) mean(aP)],[0 0.2],'-')
yticks(0:0.05:0.2)
xlim([0 1])
ylim([0 0.2])
%% Single feature analysis - grid
mASingle = mean(aS,1).*sign(mean(acuteBeta(:,:,1),2))';
[~,sortInd] = sort(abs(mASingle),'descend');
mASingleSort = mASingle(sortInd)';
sortFeat = feat(sortInd)';
for ii = 1:numel(mASingle)
    p(ii) = (1+sum(aSP(:,ii)>=mean(aS(:,ii))))/(1+size(aSP,1));
end
pAdj = p.*60;
mAS(pAdj>=0.05) = NaN;
pow = reshape(mAS(1:24),6,4);
coh = reshape(mAS(25:end),6,6);
figure
pcolor(padarray(pow,[1 1],NaN,'post')')
colormap('viridis')
% caxis([-1 1])
set(gca,'xtick',1.5:6.5,'xticklabel',...
    {'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'},...
    'ytick',1.5:4.5,'yticklabel',...
    {'ILl','ILr','NAl','NAr'})

figure
pcolor(padarray(coh,[1 1],NaN,'post')')
colormap('viridis')
% caxis([-1 1])
set(gca,'xtick',1.5:6.5,'xticklabel',...
    {'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'},...
    'ytick',1.5:6.5,'yticklabel',...
    {'ILl-ILr','ILl-NAl','ILl-NAr','NAl-ILr','ILr-NAr','NAl-NAr'})
%% Single feature peformance
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
feat = feat(subInds);
figure
hold on
for ii = 1:60
    plotSpread(aS(:,ii),'xValues',ii-0.25,'binWidth',0.1)
    plot(ii,mean(aS(:,ii)),'.k')
    plotSpread(aSP(:,ii),'binWidth',0.1,'xValues',ii+0.25,'distributionColors','r')
    plot(ii,mean(aSP(:,ii)),'.k')
end
set(gca,'xtick',1:60,'xticklabel',feat)
xtickangle(45)
%% Single feature analysis - bar plot
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
mAS = mean(aS,1);
for ii = 1:60
    [~,p(ii)] = ttest2(aS(:,ii),aSP(:,ii));
end
pAdj = p*60;
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
feat = feat(subInds);
sigFeat = feat(pAdj<=0.05);
pow = sum(pAdj(1:24)<=0.05);
coh = sum(pAdj(25:end)<=0.05);
for jj = 1:6
    powFreq(jj) = sum(pAdj(1+(jj-1):6:24)<=0.05);
    cohFreq(jj) = sum(pAdj(25+(jj-1):6:end)<=0.05);
%     powFreq(jj) = sum(mAS(1+(jj-1):6:48)>=0.6);
%     cohFreq(jj) = sum(pAdj(49+(jj-1):6:end)>=0.6);
end
figure
x = categorical({'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
x = reordercats(x,{'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
bar(x,[powFreq;cohFreq]','stacked')
box off
title('Acute features by frequency')
legend({'Power','Coherence'})
%% Frequency networks
chan = 4;
sites = {'ILL','ILR','NAcL','NAcR'};
powInd = chan*6;
cmbs = nchoosek(1:chan,2);
freqs = {'d','t','a','b','lg','hg'};
thisAS = mean(aS,1).*sign(mean(acuteBeta,2))';
[~,p] = ttest2(aS,aSP);
pAdj = p.*60;
for jj = 1:6
    adjMat = [];
    count = powInd+jj;
    for ii = 1:size(cmbs,1)
        % Make sure beta coherence is significant
        if pAdj(count) <= 0.05
            adjMat(cmbs(ii,1),cmbs(ii,2)) =  thisAS(count);
        else
            adjMat(cmbs(ii,1),cmbs(ii,2)) = 0;
        end
        count = count+6;
    end
    adjMat = [adjMat;zeros(1,chan)];
    % Find scaling factor based on minimum
    s = 1-min(min(abs(adjMat(adjMat~=0))));
    % Compute graph with scale factor
    g = graph(adjMat,sites,'upper');
    if chan == 8
    % Set up node coordinates based on regular octagon
    x = [1,sqrt(2)/2,0,-sqrt(2)/2,-1,-sqrt(2)/2,0,sqrt(2)/2];
    y = [0,sqrt(2)/2,1,sqrt(2)/2,0,-sqrt(2)/2,-1,-sqrt(2)/2];
    elseif chan == 4
       % Square
       x = [1,1,-1,-1];
       y = [1,-1,1,-1];
    end
    % Set edge color based on edge sign (+ = red; - = blue)
    edgeSign = g.Edges.Weight>=0;
    edgeColor = zeros(numel(edgeSign),3);
    edgeColor(logicFind(1,edgeSign,'=='),:) = repmat([1 0 0],sum(edgeSign),1);
    edgeColor(logicFind(0,edgeSign,'=='),:) = repmat([0 0 1],sum(~edgeSign),1);
    % Set node color based on node sign (+ = red; - = blue)
    nodeSign = thisAS(jj:6:powInd)>=0;
    nodeColor = zeros(numel(nodeSign),3);
    nodeColor(logicFind(1,nodeSign,'=='),:) = repmat([1 0 0],sum(nodeSign),1);
    nodeColor(logicFind(0,nodeSign,'=='),:) = repmat([0 0 1],sum(~nodeSign),1);
    % Check node significance
    sigNode = pAdj(jj:6:powInd) <= 0.05;
    % If not significant, set to white and size 10
    nodeColor(~sigNode,:) = repmat([1 1 1],sum(~sigNode),1);
    theseNodes = abs(thisAS(jj:6:powInd));
    theseNodes(~sigNode) = 10;
    % Plot
    figure('position',[895,400,667,571])
    h = plot(g,'XData',x,'YData',y,'markersize',theseNodes,'linewidth',abs(g.Edges.Weight)+s,'edgecolor',edgeColor,'nodecolor',nodeColor);
    title(freqs{jj})
    axis off
%     print(['LSD',freqs{jj},'.eps'],'-dwinc','-painters')
end
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
%% base vs. 24 hours later (LSD and saline separate)
[baseData,~,baseFiles] = collateData(['G:\Greenlab\data\lsd\processed\',...
    'imagWater\'],{'lsd';'sal'},{'pow','coh'},'trl','rel');
[preData,~,preFiles] = collateData(['G:\Greenlab\data\lsd\processed\',...
    'imagAcute\'],{'lsd';'sal'},{'pow','coh'},'trl','rel');
[postData,~,postFiles] = collateData(['G:\Greenlab\data\lsd\processed\',...
    'imag24HR\'],{'lsd';'sal'},{'pow','coh'},'trl','rel');
% Combine data from water (baseData) and pre injection (preData)
ids = {'_26_','_29_','_30_','_37_','_38_','_39_','_80_','_81_','_82_',...
    '_88_','_90_'};
allBaseData = cell(numel(ids),1);
for ii = 1:numel(ids)
    for jj = 1:2
        inds = logicFind(1,contains(baseFiles{jj},ids{ii}),'==');
        if ~isempty(inds)
            allBaseData{ii} = cat(1,allBaseData{ii},baseData{jj}{inds});
        end
        inds = logicFind(1,contains(preFiles{jj},ids{ii}),'==');
        if ~isempty(inds)
            allBaseData{ii} = cat(1,allBaseData{ii},preData{jj}{inds,1});
        end
    end
end
% get post data
allPostData = cell(numel(ids),1);
for ii = 1:numel(ids)
    for jj = 1:2
        inds = logicFind(1,contains(postFiles{jj},ids{ii}),'==');
        if ~isempty(inds)
            allPostData{ii} = cat(1,allPostData{ii},postData{jj}{inds});
        end
        inds = logicFind(1,contains(postFiles{jj},ids{ii}),'==');
        if ~isempty(inds)
            allPostData{ii} = cat(1,allPostData{ii},postData{jj}{inds,1});
        end
    end
end
group = [0,1,0,0,0,1,1,0,1,1,1];
LSD = logicFind(1,group,'==');
SAL = logicFind(0,group,'==');
% save('G:\GreenLab\data\lsd\preVpost.mat','allBaseData','allPostData','ids')
%% LSD pre vs. post - 20 iterations of each animal left out
for n = 1:20
    disp(n)
    thisLSDpre = cell(numel(LSD),1);
    thisLSDpost = cell(numel(LSD),1);
    for ii = 1:numel(LSD)
        thisLSDpre{ii} = allBaseData{LSD(ii)}(randperm(size(...
            allBaseData{LSD(ii)},1),900),:);
        thisLSDpost{ii} = allPostData{LSD(ii)}(randperm(size(...
            allPostData{LSD(ii)},1),900),:);
    end
    for jj = 1:numel(LSD)
        theseInds = 1:numel(LSD);
        trainX = cat(1,thisLSDpre{theseInds(~ismember(theseInds,jj))},...
            thisLSDpost{theseInds(~ismember(theseInds,jj))});
        trainY = cat(1,zeros(900*5,1),ones(900*5,1));
        testX = cat(1,thisLSDpre{theseInds(jj)},thisLSDpost{theseInds(jj)});
        testY = cat(1,zeros(900,1),ones(900,1));
        mdl = fitglm(zscore(trainX(:,subInds)),trainY,...
            'distribution','binomial');
        pred = predict(mdl,zscore(testX(:,subInds)));
        [x,y,~,aLSD(n,jj)] = perfcurve(testY,pred,1);
        xLSD(n,jj,:) = interp1(linspace(0,1,numel(x)),x,...
            linspace(0,1,1800));
        yLSD(n,jj,:) = interp1(linspace(0,1,numel(y)),y,...
            linspace(0,1,1800));
    end
end
%% permuted LSD pre vs. post
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
for ii = 1:6
    disp(ii)
    testX = cat(1,allBaseData{ii,1}(randperm(size(allBaseData{ii,1},1),...
        900),subInds),...
        allPostData{ii}(randperm(size(allPostData{ii},1),900),...
        subInds));
    testY = cat(1,zeros(900,1),ones(900,1));
    otherInd = logicFind(1,~ismember(1:6,ii),'==');
    otherInd1 = nchoosek(otherInd,3);
    for jj = 1:size(otherInd1,1)
        otherInd2(jj,:) = otherInd(~ismember(otherInd,otherInd1(jj,:)));
    end
    for jj = 1:size(otherInd1,1)
        trainX = [];
        trainY = [];
        for k = 1:3
            trainX = [trainX;...
                allPostData{otherInd1(jj,k)}(randperm(size(allPostData{otherInd1(jj,k)},1),900),subInds);...
                allBaseData{otherInd1(jj,k),1}(randperm(size(allBaseData{otherInd1(jj,k),1},1),900),subInds)];
            trainY = [trainY;ones(900,1);zeros(900,1)];
        end
        for k = 1:2
            trainX = [trainX;...
                allPostData{otherInd1(jj,k)}(randperm(size(allPostData{otherInd1(jj,k)},1),900),subInds);...
                allBaseData{otherInd1(jj,k)}(randperm(size(allBaseData{otherInd1(jj,k)},1),900),subInds)];
            trainY = [trainY;zeros(900,1);ones(900,1)];
        end
        mdl = fitglm(zscore(trainX),trainY,'Distribution','binomial');
        prob = predict(mdl,zscore(testX));
        [~,~,~,aLSDp(ii,jj,1)] = perfcurve(testY,prob,1);
        trainX = [];
        trainY = [];
        for k = 1:3
            trainX = [trainX;...
                allPostData{otherInd1(jj,k)}(randperm(size(allPostData{otherInd1(jj,k)},1),900),subInds);...
                allBaseData{otherInd1(jj,k)}(randperm(size(allBaseData{otherInd1(jj,k)},1),900),subInds)];
            trainY = [trainY;zeros(900,1);ones(900,1)];
        end
        for k = 1:2
            trainX = [trainX;...
                allPostData{otherInd1(jj,k)}(randperm(size(allPostData{otherInd1(jj,k)},1),900),subInds);...
                allBaseData{otherInd1(jj,k)}(randperm(size(allBaseData{otherInd1(jj,k)},1),900),subInds)];
            trainY = [trainY;ones(900,1);zeros(900,1)];
        end
        mdl = fitglm(zscore(trainX),trainY,'Distribution','binomial');
        prob = predict(mdl,zscore(testX));
        [~,~,~,aLSDp(ii,jj,2)] = perfcurve(testY,prob,1);
    end
end
%% SAL pre vs. post  - 20 iterations of each animal left out
for n = 1:20
    for jj = 1:numel(SAL)
        thisSALpre = cell(numel(SAL),1);
        thisSALpost = cell(numel(SAL),1);
        for ii = 1:numel(SAL)
            thisSALpre{ii} = allBaseData{SAL(ii)}(randperm(size(...
                allBaseData{SAL(ii)},1),900),:);
            thisSALpost{ii} = allPostData{SAL(ii)}(randperm(size(...
                allPostData{SAL(ii)},1),900),:);
        end
        theseInds = 1:numel(SAL);
        trainX = cat(1,thisSALpre{theseInds(~ismember(theseInds,jj))},...
            thisSALpost{theseInds(~ismember(theseInds,jj))});
        trainY = cat(1,zeros(900*4,1),ones(900*4,1));
        testX = cat(1,thisSALpre{theseInds(jj)},thisSALpost{theseInds(jj)});
        testY = cat(1,zeros(900,1),ones(900,1));
        mdl = fitglm(zscore(trainX),trainY,'distribution','binomial');
        pred = predict(mdl,zscore(testX));
        [x,y,~,aSAL(n,jj)] = perfcurve(testY,pred,1);
        xSAL(n,jj,:) = interp1(linspace(0,1,numel(x)),x,...
            linspace(0,1,1800));
        ySAL(n,jj,:) = interp1(linspace(0,1,numel(y)),y,...
            linspace(0,1,1800));
    end
end
% save('G:\GreenLab\data\lsd\preVpostModels.mat','aLSD','xLSD','yLSD',...
%     'aSAL','xSAL','ySAL','aLSDp')
%%
figure
subplot(3,1,1)
hold on
[f,xi,bw] = ksdensity(mean(aLSD,1)); 
fill(xi,f*bw,'w')
plotSpread(mean(aLSD,1)','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aLSD,[1,2]) mean(aLSD,[1,2])],[0 0.25],'-')
yticks(0:0.05:0.25)
xlim([0 1])
ylim([0 0.25])
subplot(3,1,2)
hold on
[f,xi,bw] = ksdensity(reshape(aLSDp,1,120)); 
fill(xi,f*bw,'w')
plotSpread(reshape(aLSDp,1,120)','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(reshape(aLSDp,1,120)) mean(reshape(aLSDp,1,120))],[0 0.25],'-')
yticks(0:0.05:0.25)
xlim([0 1])
ylim([0 0.25])
subplot(3,1,3)
hold on
[f,xi,bw] = ksdensity(mean(aSAL,1)); 
fill(xi,f*bw,'w')
plotSpread(mean(aSAL,1)','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aSAL,[1,2]) mean(aSAL,[1,2])],[0 0.25],'-')
yticks(0:0.05:0.25)
xlim([0 1])
ylim([0 0.25])
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
load 'G:\GreenLab\data\lsd\lsdSal24Hr-BaseImag.mat'
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
    for k = 1:20
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
        mdl = fitglm(zscore(trainX(:,subInds)),trainY,...
            'distribution','binomial','binomialSize',numel(trainY));
        prob = predict(mdl,zscore(testX(:,subInds)));
        [x24(jj,k,:),y24(jj,k,:),~,a24(jj,k)] = perfcurve(testY,prob,1);
    end
end
% Permuted
lsdInds = nchoosek(1:6,3);
salInds = nchoosek(1:5,3);
c = 1;
for ii = 1:size(lsdInds,1)
    for jj = 1:size(salInds,1)
        permInds(c,:) = cat(2,lsdInds(ii,:),salInds(jj,:));
        c = c+1;
    end
end
thisLSD = cell(1,numel(lsdDiff));
thisSal = cell(1,numel(salDiff));
for ii = 1:size(permInds,1)
    lsdInd = permInds(ii,randperm(3,3));
    salInd = permInds(ii,randperm(3,3)+3);
    lsdOtherInd = logicFind(1,~ismember(1:6,lsdInd),'==');
    salOtherInd = logicFind(1,~ismember(1:5,salInd),'==');
    for jj = 1:numel(lsdDiff)
        thisLSD{jj} = lsdDiff{jj}(randperm(size(lsdDiff{jj},1),...
            minSamp),:);
    end
    for jj = 1:numel(salDiff)
        thisSal{jj} = salDiff{jj}(randperm(size(salDiff{jj},1),...
            minSamp),:);
    end
    testX = cat(1,thisLSD{lsdInd(1)},thisSal{salInd(1)});
    testY = cat(1,ones(minSamp,1),zeros(minSamp,1)); 
    trainX = cat(1,thisLSD{lsdInd(2:3)},thisSal{salInd(2:3)},thisLSD{lsdOtherInd},thisSal{salOtherInd});
    trainY = cat(1,ones(minSamp*4,1),zeros(minSamp*5,1));
    mdl = fitglm(zscore(trainX(:,subInds)),trainY,'Distribution','binomial');
    prob = predict(mdl,zscore(testX(:,subInds)));
    [xP24(ii,:),yP24(ii,:),~,aP24(ii)] = perfcurve(testY,prob,1);
end
% save('G:\GreenLab\data\lsd\LSDvSaline24HrLogLOOImag.mat','a24','x24','y24','aP24','xP24','yP24')
%%
% load('LSDvSaline24HrLogLOOImag.mat')
% Plot
figure
subplot(2,1,1)
hold on
[f,xi,bw] = ksdensity(mean(a24,2)); 
fill(xi,f*bw,'w')
plotSpread(mean(a24,2),'xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(a24,[1,2]) mean(a24,[1,2])],[0 0.25],'-')
yticks(0:0.05:0.25)
xlim([0 1])
ylim([0 0.25])
subplot(2,1,2)
hold on
[f,xi,bw] = ksdensity(aP24); 
fill(xi,f*bw,'w')
plotSpread(aP24','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aP24) mean(aP24)],[0 0.25],'-')
yticks(0:0.05:0.25)
xlim([0 1])
ylim([0 0.25])
%% Plot in 3D feature space
% Get ordered list of single features
load('G:\GreenLab\data\lsd\acuteLSDvSalineLogLOOImag.mat','aS')
mAS = mean(aS,1);
[sorted,indSort] = sort(mAS,'descend');
% Subset features
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,...
    175:180];
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
feat = feat(subInds);
% Get acute data ready
load('G:\GreenLab\data\lsd\normAcuteImag.mat')

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
[data2,samp2,files2,time2] = collateData(['G:\GreenLab\data\lsdStim\'...
    'processed\baseImag\'],{'mPFC','in'},{'pow','coh'},'trl','');

% Combine
for ii = 1%:3
    allData{ii} = [data{1,ii};data2{1,ii}];
    allFiles{ii} = [files{ii,1};files2{ii,1}];
%     relTime{ii} = [time1.rel{1,ii};time2.rel{1,ii}];
%     absTime{ii} = [time1.abs{1,ii};time2.abs{1,ii}];   
    % check for any with fewer than 216 features
    feats = cellfun(@(x) size(x,2),allData{ii});
    allData{ii}(feats(:,1)<216,:) = [];
    allFiles{ii}(feats(:,1)<216,:) = [];
%     relTime{ii}(feats(:,1)<216,:) = [];
%     absTime{ii}(feats(:,1)<216,:) = [];
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
[lsdData,lsdSamp,lsdFiles,lsdTime] = collateData(['G:\GreenLab\data\'...
    'lsdStim\processed\postLSDImag\'],{'.mat'},{'pow','coh'},'trl','');
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
    theseFiles = logicFind(1,~cellfun(@isempty,strfind(allFiles{1},...
        ids{ii})),'==');
    for jj = theseFiles
        if size(allData{1}{jj,1},1)>1
            thisBaseDiff{ii} = [thisBaseDiff{ii};(allData{1}{jj,3}-...
                mean(allData{1}{jj,1},1))./std(allData{1}{jj,1},[],1)];
        end
    end
end
thisLSDDiff = cell(1,numel(ids));
for ii = 1:numel(ids)
    thisLSDDiff{ii} = (lsdData{1}{ii,3}-mean(lsdData{1}{ii,1},1))./...
        std(lsdData{1}{ii,1},[],1);
end
% Remove IRDM2 (animal 5) due to low samples and IRDM14 (animal 1) for too
% few features
thisBaseDiff = thisBaseDiff([2:4,6:7]);
thisLSDDiff = thisLSDDiff([2:4,6:7]);
minSamp = min([cellfun(@(x) size(x,1),thisBaseDiff),cellfun(@(x) ...
    size(x,1),thisLSDDiff)]);
% save('G:\GreenLab\data\lsdStim\lsd-baseVsal-base_zscore_stimImag_all-216feat.mat','allData',...
%     'allFiles','lsdData','lsdFiles','lsdSamp','lsdTime','minSamp',...
%     'minSamps','thisBaseDiff','thisLSDDiff')
%%
x = []; y = []; a = []; 
xS = []; yS = []; aS = []; 
xP = []; yP = []; aP = [];
xSP = []; ySP = []; aSP = [];
betaStim = [];
% Subset indices to match other electrode array
subInds = [1:12,37:54,79:90,115:126,211:216];
% subInds = [5,10:12,29,40,52,58,94,99,192,196];
subInds = 1:216;
c = 1;
for ii = 1:5
    for jj = 1:5
        cmbs(c,:) = [ii,jj];
        c = c+1;
    end
end
for ii = 1:size(cmbs,1)
    thisBase = []; thisLSD = [];
    for jj = 1:size(thisBaseDiff,2)
        thisBase{jj} = thisBaseDiff{jj}(randperm(size(thisBaseDiff{jj},1),minSamp),:);
        thisLSD{jj} = thisLSDDiff{jj}(randperm(size(thisLSDDiff{jj},1),minSamp),:);
    end
    otherLSD = logicFind(1,~ismember(1:5,cmbs(ii,1)),'==');
    otherSAL = logicFind(1,~ismember(1:5,cmbs(ii,2)),'==');
    testX = cat(1,thisBase{cmbs(ii,2)},thisLSD{cmbs(ii,1)});
    testY = cat(1,zeros(minSamp,1),ones(minSamp,1));
    trainX = cat(1,thisBase{otherSAL},thisLSD{otherLSD});
    trainY = cat(1,zeros(minSamp*4,1),ones(minSamp*4,1));
    mdl = fitglm(trainX(:,subInds),trainY,'distribution','binomial');
    prob = predict(mdl,testX(:,subInds));
    [~,~,~,aLOO(ii)] = perfcurve(testY,prob,1);
end
%%
for n = 1:20
    disp(n)
    thisBase = []; thisLSD = [];
    for ii = 1:size(thisBaseDiff,2)
        thisBase{ii} = thisBaseDiff{ii}(randperm(size(thisBaseDiff{ii},1),minSamp),:);
        thisLSD{ii} = thisLSDDiff{ii}(randperm(size(thisLSDDiff{ii},1),minSamp),:);
    end
    % 80:20
    [trainX,trainY,testX,testY] = trainTest(cat(1,thisBase{:},thisLSD{:}),...
        [zeros(minSamp*numel(thisBase),1);ones(minSamp*numel(thisLSD),1)],0.20);
    mdl = fitglm(trainX(:,subInds),trainY,'distribution','binomial','binomialSize',numel(trainY));
%     betas(:,:,n) = table2array(mdl.Coefficients);
    prob = predict(mdl,zscore(testX(:,subInds)));
    [x{n},y{n},~,a(n)] = perfcurve(testY,prob,1);
%     % Permuted
%     mdl = fitglm(trainX(:,subInds),trainY(randperm(numel(trainY),numel(trainY)),:),'distribution','binomial');
%     prob = predict(mdl,testX(:,subInds));
%     [xP{n},yP{n},~,aP(n)] = perfcurve(testY,prob,1);
%     c = 1;
%     for ii = subInds
%         mdl = fitglm(trainX(:,ii),trainY,'distribution','binomial');
%         betaStim(c,n) = table2array(mdl.Coefficients(2,1));
%         prob = predict(mdl,testX(:,ii));
%         [xS{n,c},yS{n,c},~,aS(n,c)] = perfcurve(testY,prob,1);
%         [xSP{n,c},ySP{n,c},~,aSP(n,c)]= perfcurve(testY(randperm(numel(...
%             testY),numel(testY))),prob,1);
%         c  = c+1;
%     end
end
% % Single feature analysis
% feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
%     'rNAcC'},{'d','t','a','b','lg','hg'})';
% % feat = names({'lmPFC','rmPFc','lNAcC','rNAcC'},...
% %     {'d','t','a','b','lg','hg'})';
% [sortA,sortInd] = sort(mean(aS,1)','descend');
% sortFeat = feat(sortInd);
% sortBeta = mean(betaStim(sortInd,:),2);
% save(['G:\GreenLab\data\lsdStim\'...
%     'lsdStim-base_v_salStim-base_imag_all_216feat.mat'],'x','y','a',...
%     'xP','yP','aP','xS','yS','aS','betaStim','xSP','ySP','aSP',...
%     'sortFeat','sortBeta')
%% Plot
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
figure
hold on
violin({a',aP'},'facecolor',[1,1,1])
plotSpread({a',aP'},'distributionColors',{'k','k'})
title('LSD+stim vs. SAL+stim')

figure
subplot(2,1,1)
hold on
[f,xi] = ksdensity(a); fill(xi,f/100,'w')
plotSpread(a','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(a) mean(a)],[0 0.15],'-')
yticks(0:0.05:0.15)
xlim([0 1])
ylim([0 0.15])
subplot(2,1,2)
hold on
[f,xi] = ksdensity(aP); fill(xi,f/100,'w')
plotSpread(aP','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aP) mean(aP)],[0 0.15],'-')
yticks(0:0.05:0.15)
xlim([0 1])
ylim([0 0.15])
sgtitle('LSD+stim vs. SAL+stim')
%%
mAS = mean(aS,1).*sign(mean(betaStim,2))';
[h,p] = ttest2(aS,aSP);
pAdj = p.*numel(p);
mAS(pAdj>=0.05) = NaN;
pow = reshape(mAS(1:48),6,8);
coh = reshape(mAS(49:end),6,28);

figure
pcolor(padarray(pow,[1 1],NaN,'post')')
colormap('viridis')
caxis([-1 1])
set(gca,'xtick',1.5:6.5,'xticklabel',...
    {'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'},...
    'ytick',1.5:8.5,'yticklabel',...
    {'ILl','ILr','OFCl','OFCr','NAcSl','NAcSr','NAcCl','NAcCr'})

figure
pcolor(padarray(coh,[1 1],NaN,'post')')
colormap('viridis')
caxis([-1 1])
set(gca,'xtick',1.5:6.5,'xticklabel',...
    {'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'},...
    'ytick',1.5:6.5,'yticklabel',...
    {'ILl-ILr','ILl-NAl','ILl-NAr','NAl-ILr','ILr-NAr','NAl-NAr'})
%% Keep data sets separate; saline base v stim & lsd base v stim
% Raw
[data,samp1,files,time1] = collateData(['D:\dualSite\processed\'...
    'toUseSingleImag\'],{'IL','in'},{'pow','coh'},'trl','');
% Raw
[data2,samp2,files2,time2] = collateData(['G:\GreenLab\data\lsdStim\'...
    'processed\baseImag\'],{'mPFC','in'},{'pow','coh'},'trl','');
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
%     relTime{ii}(feats(:,1)<216,:) = [];
%     absTime{ii}(feats(:,1)<216,:) = [];
    samps{ii} = cellfun(@(x) size(x,1),allData{ii});
    % minSamps from base and stim
    minSamps{ii} = min(samps{ii}(:,[1,3]),[],2);
end
% Raw
[lsdData,lsdSamp,lsdFiles,lsdTime] = collateData(['G:\GreenLab\data\'...
    'lsdStim\processed\postLSDImag\'],{'.mat'},{'pow','coh'},'trl','');
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
lsd = lsdData([2:4,6:7],:);
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
%     thisTest = randi(5,1);
%     testX = zscore(cat(1,thisBase{thisTest},thisStim{thisTest}));
%     testY = [zeros(80,1);ones(80,1)];
%     others = logicFind(1,~ismember(1:5,thisTest),'==');
%     trainX = zscore([cat(1,thisBase{others});cat(1,thisStim{others})]);
%     trainY = [zeros(80*4,1);ones(80*4,1)];
    % Normalize x values
    trainX = zscore(cat(1,thisTrainX{:}));
    trainY = cat(1,thisTrainY{:});
    testX = zscore(cat(1,thisTestX{:}));
    testY = cat(1,thisTestY{:});
%     [trainX,trainY,testX,testY] = trainTest([cat(1,saline{:,1});cat(1,saline{:,2})],[zeros(sum(cellfun(@(x) size(x,1),saline(:,1))),1);ones(sum(cellfun(@(x) size(x,1),saline(:,2))),1)],0.2);

    
    mdl = fitglm(trainX,trainY,'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX);
    [xSaline{ii},ySaline{ii},tSaline{ii},aSaline(ii)] = perfcurve(testY,prob,1);
%     betaSaline(:,:,ii) = table2array(mdl.Coefficients);
%     for n = 1:216
%         mdl = fitglm(trainX(:,n),trainY,'distribution','binomial','binomialSize',numel(trainY));
%         betaSaline(n,ii) = table2array(mdl.Coefficients(2,1));
%         prob = predict(mdl,testX(:,n));
%         [xSalineSingle{n,ii},ySalineSingle{n,ii},tSalineSingle{n,ii},aSalineSingle(n,ii)] = perfcurve(testY,prob,1);
%     end
end
%%
% LSD: base v. stim model - 80:20
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
%     betaLSD(:,:,ii) = table2array(mdl.Coefficients);
%     for n = 1:216
%         mdl = fitglm(trainX(:,n),trainY,'distribution','binomial','binomialSize',numel(trainY));
%         betaLSD(n,ii) = table2array(mdl.Coefficients(2,1));
%         prob = predict(mdl,testX(:,n));
%         [xLSDsingle{n,ii},yLSDsingle{n,ii},~,aLSDsingle(n,ii)] = perfcurve(testY,prob,1);
%         
%         prob = predict(mdl,testX(randperm(size(testX,1)),n));
%         [xLSDsingleP{n,ii},yLSDsingleP{n,ii},~,aLSDsingleP(n,ii)] = perfcurve(testY,prob,1);
%     end
end
%% lsd LOO
[thisTrainX,thisTrainY,thisTestX,thisTestY] = deal(cell(1,6));
[xLSD,yLSD,tLSD,xLSDP,yLSDP,tLSDP] = deal(cell(1,100));
[aLSD,aLSDP] = deal(zeros(1,100));
[betaLSD,aLSDsingle,aLSDsingleP] = deal(zeros(216,100));
[xLSDsingle,yLSDsingle,tLSDsingle,xLSDsingleP,yLSDsingleP,tLSDsingleP] = deal(cell(216,100));
thisBase = cell(1,5);
thisStim  = cell(1,5);
for ii = 1:200
    for jj = 1:5
        thisBase{jj} = lsd{jj,1}(randperm(size(lsd{jj,1},1),minSamp),:);
        thisStim{jj} = lsd{jj,2}(randperm(size(lsd{jj,2},1),minSamp),:);
    end
    thisTest(ii) = randi(5,1);
    testX = zscore(cat(1,thisBase{thisTest(ii)},thisStim{thisTest(ii)}));
    testY = [zeros(80,1);ones(80,1)];
    others = logicFind(1,~ismember(1:5,thisTest(ii)),'==');
    trainX = zscore([cat(1,thisBase{others});cat(1,thisStim{others})]);
    trainY = [zeros(80*4,1);ones(80*4,1)];

    mdl = fitglm(trainX,trainY,'distribution','binomial','binomialSize',numel(trainY));
    prob = predict(mdl,testX);
    [xLSD{ii},yLSD{ii},tLSD{ii},aLSD(ii)] = perfcurve(testY,prob,1);
end
%% Permuted - just shuffle group assignments (equally)
% All possible groups or 3 LSD and 2 SAL
three = nchoosek(1:5,3);
two = nchoosek(1:5,2);
c = 1;
for ii = 1:size(three,1)
    for jj = 1:size(two,1)
        cmbs(c,:) = [three(ii,:),two(jj,:)];
        c = c+1;
    end
end
%%
for ii = 1:size(cmbs,1)
    disp(ii)
    for k = 1:100
        group1 = []; group2 = [];
        for jj = 1:size(cmbs,2)
            if jj <4
                thisBase = {lsd{cmbs(ii,jj),1}(randperm(size(lsd{cmbs(ii,jj),1},1),minSamp),:)};
                thisStim = {lsd{cmbs(ii,jj),2}(randperm(size(lsd{cmbs(ii,jj),2},1),minSamp),:)};
                group1 = [group1;thisBase,thisStim];
            else
                thisBase = {saline{cmbs(ii,jj),1}(randperm(size(saline{cmbs(ii,jj),1},1),minSamp),:)};
                thisStim = {saline{cmbs(ii,jj),2}(randperm(size(saline{cmbs(ii,jj),2},1),minSamp),:)};
                group1 = [group1;thisBase,thisStim];
            end
        end
        otherInds3 = logicFind(1,~ismember(1:5,cmbs(ii,1:3)),'==');
        for jj = otherInds3
            thisBase = {lsd{jj,1}(randperm(size(lsd{jj,1},1),minSamp),:)};
            thisStim = {lsd{jj,2}(randperm(size(lsd{jj,2},1),minSamp),:)};
            group2 = [group2;thisBase,thisStim];
        end
        otherInds2 = logicFind(1,~ismember(1:5,cmbs(ii,4:5)),'==');
        for jj = otherInds2
            thisBase = {saline{jj,1}(randperm(size(saline{jj,1},1),minSamp),:)};
            thisStim = {saline{jj,2}(randperm(size(saline{jj,2},1),minSamp),:)};
            group2 = [group2;thisBase,thisStim];
        end
        group1base = cat(1,group1{:,1});
        group1stim = cat(1,group1{:,2});
        [trainX,trainY,testX,testY] = trainTest([group1base;group1stim],[zeros(400,1);ones(400,1)],0.2);
        %     thisTest = randi(5,1);
        %     testX = cat(1,group1{thisTest,1},group1{thisTest,2});
        %     testY = [zeros(80,1);ones(80,1)];
        %     others = logicFind(1,~ismember(1:5,thisTest),'==');
        %     trainX = [cat(1,group1{others,1});cat(1,group1{others,2})];
        %     trainY = [zeros(80*4,1);ones(80*4,1)];
        mdl = fitglm(trainX,trainY,'distribution','binomial','binomialSize',numel(trainY));
        pred = predict(mdl,testX);
        [~,~,~,a1(ii,k)] = perfcurve(testY,pred,1);

        group2base = cat(1,group2{:,1});
        group2stim = cat(1,group2{:,2});
        [trainX,trainY,testX,testY] = trainTest([group2base;group2stim],[zeros(400,1);ones(400,1)],0.2);
        %     thisTest = randi(5,1);
        %     testX = cat(1,group2{thisTest,1},group2{thisTest,2});
        %     testY = [zeros(80,1);ones(80,1)];
        %     others = logicFind(1,~ismember(1:5,thisTest),'==');
        %     trainX = [cat(1,group2{others,1});cat(1,group2{others,2})];
        %     trainY = [zeros(80*4,1);ones(80*4,1)];
        mdl = fitglm(trainX,trainY,'distribution','binomial','binomialSize',numel(trainY));
        pred = predict(mdl,testX);
        [~,~,~,a2(ii,k)] = perfcurve(testY,pred,1);
    end
end
%% get distribution of differences
aDiff = mean(a1,2)-mean(a2,2);
figure
[f,xi] = ksdensity(aDiff); fill(xi,f/101,'w')
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
figure
hold on
violin({aLSD',aSaline',aLSDP'},'facecolor',[1,1,1])
plotSpread({aLSD',aSaline',aLSDP'},'distributionColors',{'k','k','k'})
title('LSD+stim and SAL+stim')

figure
subplot(3,1,1)
hold on
[f,xi] = ksdensity(aLSD); fill(xi,f/100,'w')
plotSpread(aLSD','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aLSD) mean(aLSD)],[0 0.15],'-')
yticks(0:0.05:0.15)
xlim([0 1])
ylim([0 0.16])
subplot(3,1,2)
hold on
[f,xi] = ksdensity(aSaline); fill(xi,f/100,'w')
plotSpread(aSaline','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aSaline) mean(aSaline)],[0 0.15],'-')
yticks(0:0.05:0.15)
xlim([0 1])
ylim([0 0.16])
subplot(3,1,3)
hold on
[f,xi] = ksdensity(aLSDP); fill(xi,f/100,'w')
plotSpread(aLSDP','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(aLSDP) mean(aLSDP)],[0 0.15],'-')
yticks(0:0.05:0.15)
xlim([0 1])
ylim([0 0.16])
sgtitle('LSD+stim vs baseline and SAL+stim vs baseline')
%%
% Single feature analysis
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'})';

lsdMeanA = mean(aLSDsingle,2).*sign(mean(betaLSD,2));
[lsdSortA,lsdSortInd] = sort(lsdMeanA,'descend');
lsdFeat = feat(lsdSortInd)';
betaLSDsort = mean(betaLSD(lsdSortInd,:),2);
[~,lsdP] = ttest2(aLSDsingle',aLSDsingleP'); 
lsdPadj = lsdP.*216;
lsdPadjSort = lsdPadj(lsdSortInd)';
lsdMeanA(lsdPadj>0.05) = NaN;
% lsdMeanA(lsdMeanA>0.6) = NaN;

salMeanA = mean(aSalineSingle,2).*sign(mean(betaSaline,2));
[salineSortA,salineSortInd] = sort(salMeanA,'descend');
salineFeat = feat(salineSortInd)';
betaSalineSort = betaSaline(salineSortInd);
[~,salP] = ttest(aSalineSingle'-0.5);
salPadj = salP.*216;
salPadjSort = salPadj(salineSortInd)';
salMeanA(salPadj>0.05) = NaN;
% salMeanA(salMeanA>0.6) = NaN;

salPow = reshape(salMeanA(1:48),6,8);
lsdPow = reshape(lsdMeanA(1:48),6,8);
salCoh = reshape(salMeanA(49:end),6,28);
lsdCoh = reshape(lsdMeanA(49:end),6,28);


figure
pcolor(padarray(lsdPow,[1 1],NaN,'post')')
colormap('viridis')
caxis([-1 1])
set(gca,'xtick',1.5:6.5,'xticklabel',...
    {'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'},...
    'ytick',1.5:8.5,'yticklabel',...
    {'ILl','ILr','OFCl','OFCr','NAcSl','NAcSr','NAcCl','NAcCr'})

figure
pcolor(padarray(lsdCoh,[1 1],NaN,'post')')
colormap('viridis')
caxis([-1 1])
set(gca,'xtick',1.5:6.5,'xticklabel',...
    {'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'},...
    'ytick',1.5:6.5,'yticklabel',...
    {'ILl-ILr','ILl-NAl','ILl-NAr','NAl-ILr','ILr-NAr','NAl-NAr'})

figure
pcolor(padarray(salPow,[1 1],NaN,'post')')
colormap('viridis')
caxis([-1 1])
set(gca,'xtick',1.5:6.5,'xticklabel',...
    {'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'},...
    'ytick',1.5:8.5,'yticklabel',...
    {'ILl','ILr','OFCl','OFCr','NAcSl','NAcSr','NAcCl','NAcCr'})

figure
pcolor(padarray(salCoh,[1 1],NaN,'post')')
colormap('viridis')
caxis([-1 1])
set(gca,'xtick',1.5:6.5,'xticklabel',...
    {'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'},...
    'ytick',1.5:6.5,'yticklabel',...
    {'ILl-ILr','ILl-NAl','ILl-NAr','NAl-ILr','ILr-NAr','NAl-NAr'})
%%
allLSDdiff = cat(1,thisLSDDiff{1});
allSalDiff = cat(1,thisBaseDiff{1});
[mASsort,mASsortInd] = sort(abs(mAS),'descend');
figure
hold on
c = 1;
for ii = mASsortInd
    if ~isnan(mAS(ii))
        errorbar(mean(allLSDdiff(:,ii)),mAS(ii),std(allLSDdiff(:,ii),[],1),'.b','horizontal')
        errorbar(mean(allSalDiff(:,ii)),mAS(ii),std(allSalDiff(:,ii),[],1),'.r','horizontal')
%         plot([mean(allLSDdiff(:,ii)) mean(allSalDiff(:,ii))],[mAS(ii) mAS(ii)],'-k')
        c = c+1;
    end
end
%%
load('lsdStim-base_v_salStim-base_imag_all_216feat.mat')
load('lsdStimvBase_salStimvBase_imag_all_216feat.mat')

lsdMeanA = mean(aLSDsingle,2).*sign(mean(betaLSD,2));
salMeanA = mean(aSalineSingle,2).*sign(mean(betaSaline,2));
diffMean = mean(aS,1).*sign(mean(betaStim,2))';
%%
figure
plot(lsdMeanA(1:48),salMeanA(1:48),'.b')
hold on
plot(lsdMeanA(49:end),salMeanA(49:end),'.r')
plot(lsdMeanA(abs(diffMean)>0.6),salMeanA(abs(diffMean)>0.6),'ok')

inds = logicFind(1,abs(diffMean)>0.6,'==');
above = diffMean(inds);
[~,sortInd] = sort(abs(above),'descend');
sorted = above(sortInd);
sortedInd = inds(sortInd);
for ii = 1:numel(sortedInd)
    text(lsdMeanA(sortedInd(ii)),salMeanA(sortedInd(ii)),num2str(ii))
end
plot([1 -1],[1 -1],'--k')
plot([-1 1],[1 -1],'--k')
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
xlim([-0.81 0.81])
ylim([-0.81 0.81])
xlabel('LSD+stim')
ylabel('SAL+stim')
%%
lims = [-1 1;1 1;-1 -1;1 -1];
figure
for ii = 1:4
    subplot(2,2,ii)
    plot(lsdMeanA(1:48),salMeanA(1:48),'.b')
    hold on
    plot(lsdMeanA(49:end),salMeanA(49:end),'.r')
    plot(lsdMeanA(abs(diffMean)>0.6),salMeanA(abs(diffMean)>0.6),'ok')
    
    inds = logicFind(1,abs(diffMean)>0.6,'==');
    above = diffMean(inds);
    [~,sortInd] = sort(abs(above),'descend');
    sorted = above(sortInd);
    sortedInd = inds(sortInd);
    for jj = 1:numel(sortedInd)
        text(lsdMeanA(sortedInd(jj)),salMeanA(sortedInd(jj)),num2str(jj))
    end
    plot([1 -1],[1 -1],'--k')
    plot([-1 1],[1 -1],'--k')
    set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
    xlim(sort([0.46 0.81].*lims(ii,1),'ascend'))
    ylim(sort([0.46 0.81].*lims(ii,2),'ascend'))
    xlabel('LSD+stim')
    ylabel('SAL+stim')    
end
%%
for ii = 1:10
doubleHist(allLSDdiff(:,mASsortInd(ii)),allSalDiff(:,mASsortInd(ii)))
end
%%
% Count number of incresed effect; decreased effect; and reversal of effect
reversalInd = sign(mean(allLSDdiff,1))~=sign(mean(allSalDiff,1));
reversal = sum(reversalInd);
closerInd = abs(mean(allLSDdiff(:,~reversalInd),1))<abs(mean(allSalDiff(:,~reversalInd),1));
closer = sum(closerInd);
furtherInd = abs(mean(allLSDdiff(:,~reversalInd),1))>abs(mean(allSalDiff(:,~reversalInd),1));
further = sum(furtherInd);
increaseInd = mean(allLSDdiff(:,~reversalInd),1)>mean(allSalDiff(:,~reversalInd),1);
increase = sum(increaseInd);
decreaseInd = mean(allLSDdiff(:,~reversalInd),1)<mean(allSalDiff(:,~reversalInd),1);
decrease = sum(decreaseInd);
%% Single feature analysis


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
%% LSD stim + IRDM; look at features through time
[lsdIRDM,lsdIRDMsamps,lsdIRDMfiles] = collateData('D:\lsdIRDM\processed\',{'27';'34';'37'},{'pow','coh'},'trl','raw');
[baseIRDM,baseIRDMsamps,baseIRDMfiles] = collateData('D:\lsdIRDM\processed\control\',{'27';'34';'37'},{'pow','coh'},'trl','raw');
%%
% Use first baseline to z-score all subsequent baselines
for ii = 1:3
    mLSD = mean(lsdIRDM{ii}{1,1},1);
    sLSD = std(lsdIRDM{ii}{1,1},[],1);
    mBase = mean(baseIRDM{ii}{1,1},1);
    sBase = std(baseIRDM{ii}{1,1},[],1);
    c = 1;
    for jj = 2:3
        zLSD(ii,c,:) = (mean(lsdIRDM{ii}{jj,1})-mLSD)./sLSD; 
        zBase(ii,c,:) = (mean(baseIRDM{ii}{jj,1})-mBase)./sBase;
        c = c+1;
    end
end
%%
for ii = 1:216
   for jj = 1:2
      [~,p(ii,jj)] = ttest2(zLSD(:,jj,ii),zBase(:,jj,ii));
   end
end
%%
for ii = [5,10:12,29,40,52,58,94,99,192,196]
    figure
    hold on
    plot([1:2],zBase(:,:,ii),'k')
    plot([1:2],zLSD(:,:,ii),'r')
end
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

