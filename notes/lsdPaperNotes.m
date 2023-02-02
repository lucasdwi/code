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
% save('G:\GreenLab\data\lsd\normAcuteImagRest.mat','allPostLSD','allPostSal','minSamp',...
%     'data','files')
%% Get rest times
files1 = fileSearch('G:\Greenlab\data\lsd\restProcessed\','_Sal_');
files2 = fileSearch('G:\Greenlab\data\lsd\restProcessed\','_LSD_');
files = [files1,files2];
time = linspace(0,1,1000);
rests = zeros(size(files,2),1000);
restsPost = nan(size(files,2),5500);
for ii = 1:size(files,2)
    cd 'G:\Greenlab\data\lsd\restProcessed\'
    load(files{ii},'hist')
    totalTime(ii) = hist.eventTs.t{1}(end);
    starts{ii} = hist.eventTs.t{4};
    stops{ii} = hist.eventTs.t{5};
    startRel{ii} = round(starts{ii}./totalTime(ii),3);
    stopRel{ii} = round(stops{ii}./totalTime(ii),3);
    restsPost(ii,1:ceil(totalTime(ii)-inj(ii))) = 0;
    for jj = 1:size(startRel{ii})
        rests(ii,nearest_idx3(startRel{ii}(jj),time):...
            nearest_idx3(stopRel{ii}(jj),time)) = 1;
        load(['G:\GreenLab\data\lsd\processed\imagAcute\',...
            files{ii}(1:end-9),'_pre_vs_post.mat'],'hist')
        inj(ii) = hist.eventTs.t{3};
        if starts{ii}(jj) >= inj(ii)
            restsPost(ii,ceil(starts{ii}(jj)-inj(ii)):...
                ceil(stops{ii}(jj)-inj(ii))) = 1;
        end
    end
end
%%
salInd = contains(files,'_Sal_');
lsdInd = contains(files,'_LSD_');
figure
h = pcolor(padarray([restsPost(salInd,:);...
    nan(1,5500);...
    mean(restsPost(salInd,:));...
    nan(1,5500);...
    mean(restsPost(lsdInd,:));...
    nan(1,5500);...
    restsPost(lsdInd,:)],...
    [1 1],NaN,'post'));
set(h,'EdgeColor','none')
%%
figure
plot(smooth(mean(restsPost(salInd,:)),120),'b')
hold on
plot(smooth(mean(restsPost(lsdInd,:)),120),'r')
set(gca,'xtick',0:600:5000,'xticklabel',[0:600:5500]./60,'ytick',0:0.25:1)
xlim([0 5000])
xlabel('mins post injection')
ylabel('% of animals resting')
legend({'SAL','LSD'})
box off
%%
figure
imagesc([rests(salInd,:);rests(lsdInd,:)]);
hold on
for ii = 1:size(inj,2)
    plot([nearest_idx3(inj(ii),time) nearest_idx3(inj(ii),time)],...
        [ii-.5 ii+.5],'r')
end
%% Collate  data and split into pre and post injection
% [data,samp,files] = collateData('G:\Greenlab\data\lsd\restProcessed\',...
%     {'_lsd';'_sal'},{'pow','coh'},'trl','');
[data,samp,files] = collateData('G:\GreenLab\data\lsd\processed\imagAcute\',...
    {'_lsd';'_sal'},{'pow','coh'},'trl','');
%%
for c = 1:2
    for ii = 1:numel(files{c})
%         load(['G:\GreenLab\data\lsd\processed\imagAcute\',...
%             files{c}{ii}(1:end-9),'_pre_vs_post.mat'],'hist')
        load(['G:\GreenLab\data\lsd\processed\imagAcute\',...
            files{c}{ii}],'hist')
        pre{c,ii} = data{c}{ii}(samp{c}{ii}(:,1)<hist.eventTs.t{2}*2000,:);
        post{c,ii} = data{c}{ii}(samp{c}{ii}(:,1)>hist.eventTs.t{3}*2000,:);
    end
end
%% Only use data from the last 30 minutes
c = 1;
for ii = 1:2
    for jj = 1:size(data{1,ii},1)
        allData30(c,:) = [data{1,ii}(jj,1),...
            data{1,ii}{jj,2}(nearest_idx3(samp{1,ii}{jj,2}(end,1)...
            -30*60*400,samp{1,ii}{jj,2}(:,1)):end,:)];
        c = c+1;
    end
end
mBase = cellfun(@(x) mean(x,1),allData30(:,1),'UniformOutput',0);
sBase = cellfun(@(x) std(x,[],1),allData30(:,1),'UniformOutput',0);
for ii = 1:numel(mBase)
    dif30{ii} = allData30{ii,2}-repmat(mBase{ii},size(allData30{ii,2},1),1);
    dif30z{ii} = (allData30{ii,2}-...
        repmat(mBase{ii},size(allData30{ii,2},1),1))./...
        repmat(sBase{ii},size(allData30{ii,2},1),1);
end
allPostLSD30 = dif30(1:numel(files{1,1}));
allPostSal30 = dif30(numel(files{1,1})+1:end);
allPostLSD30z = dif30z(1:numel(files{1,1}));
allPostSal30z = dif30z(numel(files{1,1})+1:end);
minSamp30 = min(cellfun(@(x) size(x,1),dif30));
% save('G:\GreenLab\data\lsd\normAcuteImag30.mat','allPostLSD30',...
%     'allPostSal30','minSamp30','data','files','allData30',...
%     'allPostSal30z','allPostLSD30z')
%% Compare rest times between groups on injection day and post day
salFiles = fileSearch('G:\GreenLab\data\lsd\scoredMat\','_Sal');
salPostFiles = fileSearch('G:\GreenLab\data\lsd\scoredMat\','_PostSal');
lsdFiles = fileSearch('G:\GreenLab\data\lsd\scoredMat\','_LSD');
lsdPostFiles = fileSearch('G:\GreenLab\data\lsd\scoredMat\','_PostLSD');
for ii = 1:numel(salFiles)
    load(salFiles{ii},'eventTs')
    inds = eventInd(eventTs,{'rest'});
    salRest(ii) = sum(eventTs.t{inds(2)}-eventTs.t{inds(1)});
    salRest30(ii) = sum(eventTs.t{inds(2)}(eventTs.t{inds(2)}>...
        eventTs.t{1}(end)-30*60)-eventTs.t{inds(1)}(eventTs.t{inds(2)}>...
        eventTs.t{1}(end)-30*60));
end
for ii = 1:numel(salPostFiles)
    load(salPostFiles{ii},'eventTs')
    inds = eventInd(eventTs,{'rest'});
    salPostRest(ii) = sum(eventTs.t{inds(2)}-eventTs.t{inds(1)});
    salPostRest30(ii) = sum(eventTs.t{inds(2)}(eventTs.t{inds(2)}>...
        eventTs.t{1}(end)-30*60)-eventTs.t{inds(1)}(eventTs.t{inds(2)}>...
        eventTs.t{1}(end)-30*60));
end
for ii = 1:numel(lsdFiles)
    load(lsdFiles{ii},'eventTs')
    inds = eventInd(eventTs,{'rest'});
    lsdRest(ii) = sum(eventTs.t{inds(2)}-eventTs.t{inds(1)});
    lsdRest30(ii) = sum(eventTs.t{inds(2)}(eventTs.t{inds(2)}>...
        eventTs.t{1}(end)-30*60)-eventTs.t{inds(1)}(eventTs.t{inds(2)}>...
        eventTs.t{1}(end)-30*60));
end
for ii = 1:numel(lsdPostFiles)
    load(lsdPostFiles{ii},'eventTs')
    inds = eventInd(eventTs,{'rest'});
    lsdPostRest(ii) = sum(eventTs.t{inds(2)}-eventTs.t{inds(1)});
    lsdPostRest30(ii) = sum(eventTs.t{inds(2)}(eventTs.t{inds(2)}>...
        eventTs.t{1}(end)-30*60)-eventTs.t{inds(1)}(eventTs.t{inds(2)}>...
        eventTs.t{1}(end)-30*60));
end
figure
subplot(1,2,1)
plotSpread({salRest,lsdRest,salPostRest,lsdPostRest},...
    'DistributionMarkers','o');
ylabel('time (sec)')
set(gca,'xticklabel',{'sal acute','lsd acute','sal 24 hours',...
    'lsd 24 hours'})
title('all data')
subplot(1,2,2)
plotSpread({salRest30,lsdRest30,salPostRest30,lsdPostRest30},...
    'DistributionMarkers','o');
ylabel('time (sec)')
set(gca,'xticklabel',{'sal acute','lsd acute','sal 24 hours',...
    'lsd 24 hours'})
title('last 30 minutes')
%% power and coherence changes in LSD relative to SAL
salFiles = fileSearch('G:\GreenLab\data\lsd\processed\imagAcute\','Sal');
[salPow,salCoh] = deal(cell(numel(salFiles),2));
for ii = 1:numel(salFiles)
    load(salFiles{ii},'psdTrls','coh','trls','LFPTs')
    inds30 = nearest_idx3(nearest_idx3(LFPTs.tvec(end)-30*60,...
        LFPTs.tvec),trls{1,2}.sampleinfo(:,1));
    salPow{ii,1} = psdTrls{1}.Pow;
    salPow{ii,2} = psdTrls{2}.Pow;
    salPow{ii,3} = psdTrls{2}.Pow(:,:,inds30:end);

    [b,c,t] = size(psdTrls{1}.bandPow);
    salBand{ii,1} = reshape(psdTrls{1}.bandPow,b*c,t)';
    [b,c,t] = size(psdTrls{2}.bandPow);
    salBand{ii,2} = reshape(psdTrls{2}.bandPow,b*c,t)';

    salCoh{ii,1} = coh{1}.Cxy;
    salCoh{ii,2} = coh{2}.Cxy;
    salCoh{ii,3} = coh{2}.Cxy(:,:,inds30:end);
end
lsdFiles = fileSearch('G:\GreenLab\data\lsd\processed\imagAcute\','LSD');
[lsdPow,lsdCoh] = deal(cell(numel(lsdFiles),2));
for ii = 1:numel(lsdFiles)
    load(lsdFiles{ii},'psdTrls','coh','trls','LFPTs')
    inds30 = nearest_idx3(nearest_idx3(LFPTs.tvec(end)-30*60,...
        LFPTs.tvec),trls{1,2}.sampleinfo(:,1));
    lsdPow{ii,1} = psdTrls{1}.Pow;
    lsdPow{ii,2} = psdTrls{2}.Pow;
    lsdPow{ii,3} = psdTrls{2}.Pow(:,:,inds30);

    [b,c,t] = size(psdTrls{1}.bandPow);
    lsdBand{ii,1} = reshape(psdTrls{1}.bandPow,b*c,t)';
    [b,c,t] = size(psdTrls{2}.bandPow);
    lsdBand{ii,2} = reshape(psdTrls{2}.bandPow,b*c,t)';

    lsdCoh{ii,1} = coh{1}.Cxy;
    lsdCoh{ii,2} = coh{2}.Cxy;
    lsdCoh{ii,3} = coh{2}.Cxy(:,:,inds30:end);
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
LSDpost30 = cat(3,lsdPow{:,3});
allLSDpost = cat(3,lsdPow{:,2});
mLSDpre = mean(allLSDpre,3);
sLSDpre = std(allLSDpre,[],3);
sLSDpost = std(allLSDpost,[],3);
mLSDpost = mean(allLSDpost,3);
mLSDpost30 = mean(LSDpost30,3);
sLSDpost30 = std(LSDpost30,[],3);
lsdDiff = mLSDpost-mLSDpre;
lsdDiffS = (abs(std(allLSDpost,[],3)+std(allLSDpre,[],3)-2*...
    std(cat(3,allLSDpre,allLSDpost),[],3))).^0.5;

allSALpre = cat(3,salPow{:,1});
allSALpost = cat(3,salPow{:,2});
SALpost30 = cat(3,salPow{:,3});
mSALpre = mean(allSALpre,3);
mSALpost = mean(allSALpost,3);
mSALpost30 = mean(SALpost30,3);
sSALpre = std(allSALpre,[],3);
sSALpost = std(allSALpost,[],3);
sSALpost30 = std(SALpost30,[],3);
salDiff = mSALpost-mSALpre;
salDiffS = (abs(std(allSALpost,[],3)+std(allSALpre,[],3)-2*...
    std(cat(3,allSALpre,allSALpost),[],3))).^0.5;
x = (abs(std(cat(3,allSALpost,allLSDpost),[],3)+std(cat(3,allSALpost,allLSDpost),[],3)-2*...
    std(cat(3,allLSDpre,allLSDpost,allSALpre,allSALpost),[],3))).^0.5;
lsdSalDiff = lsdDiff-salDiff;
lsdSalDiffS = (abs(lsdDiffS+salDiffS-2*(x))).^0.5;
sites = {'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'};
%%
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
%%
figure
for ii = 1:8
subplot(4,2,ii)
    hold on
    shadedErrorBar(1:100,mLSDpre(ii,:),sLSDpre(ii,:),'k',1)
    shadedErrorBar(1:100,mLSDpost(ii,:),sLSDpost(ii,:),'r',1)
    plot([1 1],[-80 -20],':k')
    plot([4 4],[-80 -20],':k')
    plot([5 5],[-80 -20],':k')
    plot([10 10],[-80 -20],':k')
    plot([11 11],[-80 -20],':k')
    plot([14 14],[-80 -20],':k')
    plot([15 15],[-80 -20],':k')
    plot([30 30],[-80 -20],':k')
    plot([45 45],[-80 -20],':k')
    plot([65 65],[-80 -20],':k')
    plot([70 70],[-80 -20],':k')
    plot([90 90],[-80 -20],':k')
    set(gca,'xtick',0:10:100)
    xlabel('frequency (Hz)')
    ylabel('a.u.')
    title(sites{ii})
end
figure
for ii = 1:8
subplot(4,2,ii)
    hold on
    shadedErrorBar(1:100,mSALpre(ii,:),sSALpre(ii,:),'k',1)
    shadedErrorBar(1:100,mSALpost(ii,:),sSALpost(ii,:),'b:',1)
    plot([1 1],[-80 -20],':k')
    plot([4 4],[-80 -20],':k')
    plot([5 5],[-80 -20],':k')
    plot([10 10],[-80 -20],':k')
    plot([11 11],[-80 -20],':k')
    plot([14 14],[-80 -20],':k')
    plot([15 15],[-80 -20],':k')
    plot([30 30],[-80 -20],':k')
    plot([45 45],[-80 -20],':k')
    plot([65 65],[-80 -20],':k')
    plot([70 70],[-80 -20],':k')
    plot([90 90],[-80 -20],':k')
    set(gca,'xtick',0:10:100)
    xlabel('frequency (Hz)')
    ylabel('a.u.')
    title(sites{ii})
end
%%
figure
for ii = 1:8
subplot(4,2,ii)
    hold on
    shadedErrorBar(1:100,mLSDpre(ii,:),sLSDpre(ii,:),'k',1)
    shadedErrorBar(1:100,mLSDpost30(ii,:),sLSDpost30(ii,:),'r',1)
    plot([1 1],[-80 -20],':k')
    plot([4 4],[-80 -20],':k')
    plot([5 5],[-80 -20],':k')
    plot([10 10],[-80 -20],':k')
    plot([11 11],[-80 -20],':k')
    plot([14 14],[-80 -20],':k')
    plot([15 15],[-80 -20],':k')
    plot([30 30],[-80 -20],':k')
    plot([45 45],[-80 -20],':k')
    plot([65 65],[-80 -20],':k')
    plot([70 70],[-80 -20],':k')
    plot([90 90],[-80 -20],':k')
    set(gca,'xtick',0:10:100)
    xlabel('frequency (Hz)')
    ylabel('a.u.')
    title(sites{ii})
end
figure
for ii = 1:8
subplot(4,2,ii)
    hold on
    shadedErrorBar(1:100,mSALpre(ii,:),sSALpre(ii,:),'k',1)
    shadedErrorBar(1:100,mSALpost30(ii,:),sSALpost30(ii,:),'b:',1)
    plot([1 1],[-80 -20],':k')
    plot([4 4],[-80 -20],':k')
    plot([5 5],[-80 -20],':k')
    plot([10 10],[-80 -20],':k')
    plot([11 11],[-80 -20],':k')
    plot([14 14],[-80 -20],':k')
    plot([15 15],[-80 -20],':k')
    plot([30 30],[-80 -20],':k')
    plot([45 45],[-80 -20],':k')
    plot([65 65],[-80 -20],':k')
    plot([70 70],[-80 -20],':k')
    plot([90 90],[-80 -20],':k')
    set(gca,'xtick',0:10:100)
    xlabel('frequency (Hz)')
    ylabel('a.u.')
    title(sites{ii})
end
%% coh - coherogram differences
allLSDpre = cat(3,lsdCoh{:,1});
allLSDpost = cat(3,lsdCoh{:,2});
LSDpost30 = cat(3,lsdCoh{:,3});
% allLSDpre = cellfun(@(x) mean(x,3),lsdCoh(:,1),'UniformOutput',0);
% allLSDpre = cat(3,allLSDpre{:});
%
% allLSDpost = cellfun(@(x) mean(x,3),lsdCoh(:,2),'UniformOutput',0);
% allLSDpost = cat(3,allLSDpost{:});
mLSDpre = mean(allLSDpre,3);
mLSDpost = mean(allLSDpost,3);
mLSDpost30 = mean(LSDpost30,3);
sLSDpost30 = std(LSDpost30,[],3);
sLSDpre = std(allLSDpre,[],3);
lsdDiff = mLSDpost-mLSDpre;
% Linear propogation of error
lsdDiffS = (abs(std(allLSDpost,[],3)+std(allLSDpre,[],3)-2*...
    std(cat(3,allLSDpre,allLSDpost),[],3))).^0.5;

allSALpre = cat(3,salCoh{:,1});
allSALpost = cat(3,salCoh{:,2});
SALpost30 = cat(3,salCoh{:,3});
% allSALpre = cellfun(@(x) mean(x,3),salCoh(:,1),'UniformOutput',0);
% allSALpre = cat(3,allSALpre{:});
%
% allSALpost = cellfun(@(x) mean(x,3),salCoh(:,2),'UniformOutput',0);
% allSALpost = cat(3,allSALpost{:});
mSALpre = mean(allSALpre,3);
mSALpost = mean(allSALpost,3);
mSALpost30 = mean(SALpost30,3);
sSALpost30 = std(SALpost30,[],3);
sSALpre = std(allSALpre,[],3);
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
%%
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
%%
figure
for ii = 1:size(cmbs,1)
    subplot(3,2,ii)
    hold on
    shadedErrorBar(1:100,decimate(mLSDpost30(ii,:),5),...
        decimate(sLSDpost30(ii,:),5),'r',1)
    shadedErrorBar(1:100,decimate(mSALpost30(ii,:),5),...
        decimate(sSALpost30(ii,:),5),'b',1)
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
%%
figure
for ii = 1:size(cmbs,1)
    subplot(3,2,ii)
    hold on
%     shadedErrorBar(1:100,decimate(mLSDpost30(ii,:),5),...
%         decimate(sLSDpost30(ii,:),5),'r',1)
%     shadedErrorBar(1:100,decimate(mLSDpre(ii,:),5),...
%         decimate(sLSDpre(ii,:),5),'k',1)
    plot(1:100,decimate(mLSDpost30(ii,:),5),'r')
    plot(1:100,decimate(mLSDpre(ii,:),5),'k')
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
    ylim([-0.2 0.2])
    xlabel('frequency (Hz)')
    ylabel('a.u.')
    title([sites{cmbs(ii,1)},'-',sites{cmbs(ii,2)}])
end
%%
figure
for ii = 1:size(cmbs,1)
    subplot(3,2,ii)
    hold on
%     shadedErrorBar(1:100,decimate(mSALpost30(ii,:),5),...
%         decimate(sLSDpost30(ii,:),5),'b',1)
%     shadedErrorBar(1:100,decimate(mSALpre(ii,:),5),...
%         decimate(sLSDpre(ii,:),5),'k',1)
        plot(1:100,decimate(mSALpost30(ii,:),5),'b')
    plot(1:100,decimate(mSALpre(ii,:),5),'k')
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
    ylim([-0.2 0.2])
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
% load('G:\GreenLab\data\lsd\normAcuteImag30.mat')
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
        thisLSD24 = [thisLSD24;allPostLSD30z{ii}(randperm(size(allPostLSD30z{ii},1),...
            minSamp30),:)];
    end
    lsdTest = allPostLSD30z{testLSD}(randperm(size(allPostLSD30z{testLSD},1),...
        minSamp30),:);
    % Saline data
    trainSal = logicFind(1,~ismember(1:salN,inds(n,2)),'==');
    testSal = logicFind(0,~ismember(1:salN,inds(n,2)),'==');
    thisSal = [];
    for ii = 1:numel(allPostSal30z)
        rng((ii+numel(allPostLSD30z))*count)
        thisSal = [thisSal;allPostSal30z{ii}(randperm(size(allPostSal30z{ii},1),...
            minSamp30),:)];
    end
    salTest = allPostSal30z{testSal}(randperm(size(allPostSal30z{testSal},1),...
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
thisLSD = cell(1,numel(allPostLSD30z));
thisSal = cell(1,numel(allPostSal30z));
for ii = 1:size(permInds,1)
    disp(ii)
    lsdInd = permInds(ii,randperm(3,3));
    salInd = permInds(ii,randperm(3,3)+3);
    lsdOtherInd = logicFind(1,~ismember(1:6,lsdInd),'==');
    salOtherInd = logicFind(1,~ismember(1:5,salInd),'==');
    for jj = 1:numel(allPostLSD30z)
        thisLSD{jj} = allPostLSD30z{jj}(randperm(size(allPostLSD30z{jj},1),...
            minSamp30),:);
    end
    for jj = 1:numel(allPostSal30z)
        thisSal{jj} = allPostSal30z{jj}(randperm(size(allPostSal30z{jj},1),...
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
% save('acuteLSDvSalineLogLOOImag30z.mat','a','aP','aS','aSP','acuteBeta',...
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
%% Plot best and worst features
load('normAcuteImag30.mat')
load('acuteLSDvSalineLogLOOImag30.mat')
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
feat = feat(subInds);
mASingle = mean(aS,1).*sign(mean(acuteBeta(:,:,1),2))';
[~,sortInd] = sort(abs(mASingle),'descend');
mASingleSort = mASingle(sortInd)';
sortFeat = feat(sortInd)';
allLSD = cat(1,allPostLSD30{:});
allLSD = allLSD(:,subInds);
allSAL = cat(1,allPostSal30{:});
allSAL = allSAL(:,subInds);
mSAL = cellfun(@(x) mean(x,1),allPostSal30,'UniformOutput',false);
mSAL = cat(1,mSAL{:});
mLSD = cellfun(@(x) mean(x,1),allPostLSD30,'UniformOutput',false);
mLSD = cat(1,mLSD{:});
%%
figure
subplot(2,2,1)
hold on
plotSpread(allSAL(:,2),'distributionColors','b','xvalues',1,'binWidth',.001);
plotSpread(mSAL(:,2),'distributionColors','#00ffff','xvalues',1,'binWidth',0.25)
plotSpread(allLSD(:,2),'distributionColors','r','xvalues',2,'binWidth',.001)
p = plotSpread(mLSD(:,2),'distributionColors','#ff9cff','xvalues',2,'binWidth',0.25);
% ylim([-75 175])
xlim([0.5 2.5])
subplot(2,2,2)
hold on
plot(squeeze(mean(xS(:,2,:),1)),squeeze(mean(yS(:,2,:),1)),'k')
plot(squeeze(mean(xSP(:,2,:),1)),squeeze(mean(yP,1)),'k--')
subplot(2,2,3)
hold on
plotSpread(allSAL(:,24),'distributionColors','b','xvalues',1,'binWidth',0.001)
plotSpread(mSAL(:,24),'distributionColors','#00ffff','xvalues',1,'binWidth',0.5)
plotSpread(allLSD(:,24),'distributionColors','r','xvalues',2,'binWidth',0.001)
plotSpread(mLSD(:,24),'distributionColors','#ff9cff','xvalues',2,'binWidth',0.5)
xlim([0.5 2.5])
% ylim([-200 400])
subplot(2,2,4)
hold on
plot(squeeze(mean(xS(:,24,:),1)),squeeze(mean(yS(:,24,:),1)),'k')
plot(squeeze(mean(xSP(:,24,:),1)),squeeze(mean(yP,1)),'k--')
%% Plot all features
% load('normAcuteImag30.mat')
subInds = [1:6,37:42,19:24,43:48,79:84,61:66,85:90,169:174,211:216,175:180];
feat = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},...
    {'d','t','a','b','lg','hg'});
feat = feat(subInds);

allSAL = cat(1,allPostSal30{:});
allLSD = cat(1,allPostLSD30{:});

these = cell(1);
for ii = subInds
    these = {these{:},allSAL(:,ii),allLSD(:,ii)};
end
%%
figure
subplot(2,2,1)
violin(these(:,2:13));
ylim([-250 375])
subplot(2,2,2)
violin(these(:,14:25));
ylim([-250 375])
subplot(2,2,3)
violin(these(:,26:37));
ylim([-250 375])
subplot(2,2,4)
violin(these(:,38:49));
ylim([-250 375])
%%
figure
subplot(3,2,1)
violin(these(:,50:61));
set(gca,'xtick',1:2:12,'xticklabel',feat(25:30))
ylim([-1 1])
subplot(3,2,2)
violin(these(:,62:73));
set(gca,'xtick',1:2:12,'xticklabel',feat(31:36))
ylim([-1 1])
subplot(3,2,3)
violin(these(:,74:85));
set(gca,'xtick',1:2:12,'xticklabel',feat(37:42))
ylim([-1 1])
subplot(3,2,4)
violin(these(:,86:97));
set(gca,'xtick',1:2:12,'xticklabel',feat(43:48))
ylim([-1 1])
subplot(3,2,5)
violin(these(:,98:109));
set(gca,'xtick',1:2:12,'xticklabel',feat(49:54))
ylim([-1 1])
subplot(3,2,6)
violin(these(:,110:121));
set(gca,'xtick',1:2:12,'xticklabel',feat(55:60))
ylim([-1 1])
%%
figure
subplot(2,1,1)
violin(allSAL(:,subInds(1:24)),'facecolor','b');
set(gca,'xtick',1:24,'xticklabel',feat(1:24))
ylim([-200 375])
subplot(2,1,2)
violin(allLSD(:,subInds(1:24)),'facecolor','r');
set(gca,'xtick',1:24,'xticklabel',feat(1:24))
ylim([-200 375])
figure
subplot(2,1,1)
violin(allSAL(:,subInds(25:end)),'facecolor','b');
set(gca,'xtick',1:36,'xticklabel',feat(25:36))
ylim([-1 1])
subplot(2,1,2)
violin(allLSD(:,subInds(25:end)),'facecolor','r');
set(gca,'xtick',1:36,'xticklabel',feat(25:36))
ylim([-1 1])
%%
figure
c = 1;
for ii = 1:24
    hold on
    plotSpread(allSAL(:,subInds(ii)),'binWidth',0.001,...
        'xvalues',c,'distributionColors','b')
    [f,xi,bw] = ksdensity(allSAL(:,subInds(ii)));
    plot((f*bw)*4.5+c,xi,'k','lineWidth',2)
    c = c+1;
    plotSpread(allLSD(:,subInds(ii)),'binWidth',0.001,...
        'xvalues',c,'distributionColors','r');
    [f,xi,bw] = ksdensity(allLSD(:,subInds(ii)));
    plot((f*bw)*4.5+c,xi,'k','lineWidth',2)
    c = c+1;
end
set(gca,'xtick',1.5:2:48,'xticklabel',feat)
figure
c = 1;
for ii = 25:60
    hold on
    plotSpread(allSAL(:,subInds(ii)),'binWidth',0.001,...
        'xvalues',c,'distributionColors','b')
    [f,xi,bw] = ksdensity(allSAL(:,subInds(ii)));
    plot((f*bw)*4.5+c,xi,'k','lineWidth',2)
    c = c+1;
    plotSpread(allLSD(:,subInds(ii)),'binWidth',0.001,...
        'xvalues',c,'distributionColors','r');
    [f,xi,bw] = ksdensity(allLSD(:,subInds(ii)));
    plot((f*bw)*4.5+c,xi,'k','lineWidth',2)
    c = c+1;
end
set(gca,'xtick',1.5:2:72,'xticklabel',feat(25:60))
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
%% Modeling under run24HrLSDvSaline.m
% cd G:\GreenLab\data\lsd\LSDvSal24Hr\
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
load('LSDvSaline24HrLogLOOImag.mat')
% Plot
figure
subplot(2,1,1)
hold on
[f,xi,bw] = ksdensity(mean(a24,2));
fill(xi,f*bw,'w')
plotSpread(mean(a24,2),'xyori','flipped','xvalues',0.05,'spreadWidth',0.1,'distributionMarkers','.')
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
% [data,samp1,files,time1] = collateData(['G:\GreenLab\data\dualSite\',...
%     'processed\'],{'IL','in'},{'pow','coh'},'trl','');
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
%% Load pow and coh data
cd 'G:\GreenLab\data\lsdStim\processed\baseImag\'
for ii = 1:size(files2{1},1)
    load(files2{1}{ii},'psdTrls','coh')
    allSALPowPre{ii} = psdTrls{1}.Pow;
    allSALCohPre{ii} = coh{1}.Cxy;
    allSALPowStim{ii} = psdTrls{3}.Pow;
    allSALCohStim{ii} = coh{3}.Cxy;
    allSALPowPost{ii} = psdTrls{2}.Pow;
    allSALCohPost{ii} = coh{2}.Cxy;
end
cd 'G:\GreenLab\data\lsdStim\processed\postLSDImag\'
for ii = 1:size(lsdFiles{1},1)
    load(lsdFiles{1}{ii},'psdTrls','coh')
    allLSDPowPre{ii} = psdTrls{1}.Pow;
    allLSDCohPre{ii} = coh{1}.Cxy;
    allLSDPowStim{ii} = psdTrls{3}.Pow;
    allLSDCohStim{ii} = coh{3}.Cxy;
    allLSDPowPost{ii} = psdTrls{2}.Pow;
    allLSDCohPost{ii} = coh{2}.Cxy;
end
%%
ids = {'IRDM15';'IRDM16';'IRDM21';'IRDM5';'IRDM6'};
for ii = 1:numel(ids)
    salInds = contains(files2{1},ids{ii});
    lsdInds = contains(lsdFiles{1},ids{ii});
    
    salPowPre{ii} = cat(3,allSALPowPre{salInds});
    lsdPowPre{ii} = cat(3,allLSDPowPre{lsdInds});
    salPowStim{ii} = cat(3,allSALPowStim{salInds});
    lsdPowStim{ii} = cat(3,allLSDPowStim{lsdInds});
    salPowPost{ii} = cat(3,allSALPowPost{salInds});
    lsdPowPost{ii} = cat(3,allLSDPowPost{lsdInds});

    salCohPre{ii} = cat(3,allSALCohPre{salInds});
    lsdCohPre{ii} = cat(3,allLSDCohPre{lsdInds});
    salCohStim{ii} = cat(3,allSALCohStim{salInds});
    lsdCohStim{ii} = cat(3,allLSDCohStim{lsdInds});
    salCohPost{ii} = cat(3,allSALCohPost{salInds});
    lsdCohPost{ii} = cat(3,allLSDCohPost{lsdInds});
end
%%
salPowPreM = mean(cat(3,salPowPre{:}),3);
salPowPreS = std(cat(3,salPowPre{:}),[],3);
salPowPostM = mean(cat(3,salPowPost{:}),3);
salPowPostS = std(cat(3,salPowPost{:}),[],3);
salPowStimM = mean(cat(3,salPowStim{:}),3);
salPowStimS = std(cat(3,salPowStim{:}),[],3);

lsdPowPreM = mean(cat(3,lsdPowPre{:}),3);
lsdPowPreS = std(cat(3,lsdPowPre{:}),[],3);
lsdPowPostM = mean(cat(3,lsdPowPost{:}),3);
lsdPowPostS = std(cat(3,lsdPowPost{:}),[],3);
lsdPowStimM = mean(cat(3,lsdPowStim{:}),3);
lsdPowStimS = std(cat(3,lsdPowStim{:}),[],3);

salCohPreM = mean(cat(3,salCohPre{:}),3);
salCohPreS = std(cat(3,salCohPre{:}),[],3);
salCohPostM = mean(cat(3,salCohPost{:}),3);
salCohPostS = std(cat(3,salCohPost{:}),[],3);
salCohStimM = mean(cat(3,salCohStim{:}),3);
salCohStimS = std(cat(3,salCohStim{:}),[],3);

lsdCohPreM = mean(cat(3,lsdCohPre{:}),3);
lsdCohPreS = std(cat(3,lsdCohPre{:}),[],3);
lsdCohPostM = mean(cat(3,lsdCohPost{:}),3);
lsdCohPostS = std(cat(3,lsdCohPost{:}),[],3);
lsdCohStimM = mean(cat(3,lsdCohStim{:}),3);
lsdCohStimS = std(cat(3,lsdCohStim{:}),[],3);

sites = {'lIL','rIL','lOFC','rOFC','lNAcS','rNAcS','lNAcC','rNAcC'};
cmbs = nchoosek(1:8,2);
% cInds = [1,4,5,10,11,23];
% cmbs = cmbs(cInds,:);
%% raw data
figure
for ii = 1:8
    subplot(4,2,ii)
    hold on
    plot(salPowPreM(ii,:),'k')
    plot(salPowStimM(ii,:),'b')
    plot([1 1],[-70 -25],':k')
    plot([5 5],[-70 -25],':k')
    plot([10 10],[-70 -25],':k')
    plot([15 15],[-70 -25],':k')
    plot([30 30],[-70 -25],':k')
    plot([45 45],[-70 -25],':k')
    plot([65 65],[-70 -25],':k')
    plot([70 70],[-70 -25],':k')
    plot([90 90],[-70 -25],':k')
%     plot(salPowPostM(ii,:),'b--')
    ylim([-70 -25])
    title(sites{ii})
end
figure
for ii = 1:size(cmbs,1)
    subplot(7,4,ii)
    hold on
    plot(decimate(salCohPreM(ii,:),5),'k')
    plot(decimate(salCohStimM(ii,:),5),'b')
    plot([1 1],[-0.2 0.15],':k')
    plot([5 5],[-0.2 0.15],':k')
    plot([10 10],[-0.2 0.15],':k')
    plot([15 15],[-0.2 0.15],':k')
    plot([30 30],[-0.2 0.15],':k')
    plot([45 45],[-0.2 0.15],':k')
    plot([65 65],[-0.2 0.15],':k')
    plot([70 70],[-0.2 0.15],':k')
    plot([90 90],[-0.2 0.15],':k')
%     plot(decimate(salCohPostM(cInds(ii),:),5),'b--')
    title([sites{cmbs(ii,1)},'-',sites{cmbs(ii,2)}])
    ylim([-0.2 0.15])
end
figure
for ii = 1:8
    subplot(4,2,ii)
    hold on
    plot(salPowPreM(ii,:),'k')
    plot(salPowStimM(ii,:),'r')
    plot([1 1],[-70 -25],':k')
    plot([5 5],[-70 -25],':k')
    plot([10 10],[-70 -25],':k')
    plot([15 15],[-70 -25],':k')
    plot([30 30],[-70 -25],':k')
    plot([45 45],[-70 -25],':k')
    plot([65 65],[-70 -25],':k')
    plot([70 70],[-70 -25],':k')
    plot([90 90],[-70 -25],':k')
%     plot(salPowPostM(ii,:),'r--')
    title(sites{ii})
    ylim([-70 -25])
end
figure
for ii = 1:size(cmbs,1)
    subplot(7,4,ii)
    hold on
    plot(decimate(lsdCohPreM(ii,:),5),'k')
    plot(decimate(lsdCohStimM(ii,:),5),'r')
    plot([1 1],[-0.2 0.15],':k')
    plot([5 5],[-0.2 0.15],':k')
    plot([10 10],[-0.2 0.15],':k')
    plot([15 15],[-0.2 0.15],':k')
    plot([30 30],[-0.2 0.15],':k')
    plot([45 45],[-0.2 0.15],':k')
    plot([65 65],[-0.2 0.15],':k')
    plot([70 70],[-0.2 0.15],':k')
    plot([90 90],[-0.2 0.15],':k')
%     plot(decimate(lsdCohPostM(cInds(ii),:),5),'r--')
    title([sites{cmbs(ii,1)},'-',sites{cmbs(ii,2)}])
    ylim([-0.2 0.15])
end
%% subtraction figure
salPowDiff = salPowStimM-salPowPreM;
salPowDiffS = (abs(std(cat(3,salPowStim{:}),[],3)+std(cat(3,salPowPre{:}),[],3)-2*...
    std(cat(3,salPowStim{:},salPowPre{:}),[],3))).^0.5;
lsdPowDiff = lsdPowStimM-lsdPowPreM;
lsdPowDiffS = (abs(std(cat(3,lsdPowStim{:}),[],3)+std(cat(3,lsdPowPre{:}),[],3)-2*...
    std(cat(3,lsdPowStim{:},lsdPowPre{:}),[],3))).^0.5;
x = (abs(std(cat(3,salPowStim{:},lsdPowStim{:}),[],3)+std(cat(3,salPowStim{:},lsdPowStim{:}),[],3)-2*...
    std(cat(3,lsdPowPre{:},lsdPowStim{:},salPowPre{:},salPowStim{:}),[],3))).^0.5;
lsdSalPowDiff = lsdPowDiff-salPowDiff;
lsdSalPowDiffS = (abs(lsdPowDiffS+salPowDiffS-2*(x))).^0.5;

salCohDiff = salCohStimM-salCohPreM;
salCohDiffS = (abs(std(cat(3,salCohStim{:}),[],3)+std(cat(3,salCohPre{:}),[],3)-2*...
    std(cat(3,salCohStim{:},salCohPre{:}),[],3))).^0.5;
lsdCohDiff = lsdCohStimM-lsdCohPreM;
lsdCohDiffS = (abs(std(cat(3,lsdCohStim{:}),[],3)+std(cat(3,lsdCohPre{:}),[],3)-2*...
    std(cat(3,lsdCohStim{:},lsdCohPre{:}),[],3))).^0.5;
x = (abs(std(cat(3,salCohStim{:},lsdCohStim{:}),[],3)+std(cat(3,salCohStim{:},lsdCohStim{:}),[],3)-2*...
    std(cat(3,lsdCohPre{:},lsdCohStim{:},salCohPre{:},salCohStim{:}),[],3))).^0.5;
lsdSalCohDiff = lsdCohDiff-salCohDiff;
lsdSalCohDiffS = (abs(lsdCohDiffS+salCohDiffS-2*(x))).^0.5;
%%
figure
for ii = 1:8
    subplot(2,4,ii)
    hold on
    shadedErrorBar(1:100,lsdPowDiff(ii,:),lsdPowDiffS(ii,:),{'r','linewidth',1.5},1)
    shadedErrorBar(1:100,salPowDiff(ii,:),salPowDiffS(ii,:),{'b','linewidth',1.5},1)
    plot([1 1],[-10.5 1],':k')
    plot([5 5],[-10.5 1],':k')
    plot([10 10],[-10.5 1],':k')
    plot([15 15],[-10.5 1],':k')
    plot([30 30],[-10.5 1],':k')
    plot([45 45],[-10.5 1],':k')
    plot([65 65],[-10.5 1],':k')
    plot([70 70],[-10.5 1],':k')
    plot([90 90],[-10.5 1],':k')
    set(gca,'xtick',0:10:100)
%     xlabel('frequency (Hz)')
%     ylabel('a.u.')
    ylim([-10.5 1])
    title(sites{ii})
    set(gca,'xticklabel',[],'yticklabel',[],'xcolor','k','ycolor','k',...
        'linewidth',0.75)
end
% for ii = 1:8
%     subplot(4,4,8+ii)
%     hold on
%     shadedErrorBar(1:100,lsdSalPowDiff(ii,:),lsdSalPowDiffS(ii,:),'k',1)
%     set(gca,'xtick',0:10:100)
%     ylim([-4 4])
%     plot([0 100],[0 0],'--k')
% end
%%
sites = {'lIL','rIL','lOFC','rOFC','lNAcS','rNAcS','lNAcC','rNAcC'};
cmbs = nchoosek(1:8,2);
% cInds = [1,4,5,10,11,23];
% cmbs = cmbs(cInds,:);
figure
for ii = 1:size(cmbs,1)
    subplot(7,4,ii)
    hold on
    shadedErrorBar(1:100,decimate(lsdCohDiff(ii,:),5),...
        decimate(lsdCohDiffS(ii,:),5),{'r','linewidth',1.5},1)
    shadedErrorBar(1:100,decimate(salCohDiff(ii,:),5),...
        decimate(salCohDiffS(ii,:),5),{'b','linewidth',1.5},1)
    plot([1 1],[-6.5 2],':k')
    plot([5 5],[-6.5 2],':k')
    plot([10 10],[-6.5 2],':k')
    plot([15 15],[-6.5 2],':k')
    plot([30 30],[-6.5 2],':k')
    plot([45 45],[-6.5 2],':k')
    plot([65 65],[-6.5 2],':k')
    plot([70 70],[-6.5 2],':k')
    plot([90 90],[-6.5 2],':k')
    set(gca,'xtick',0:10:100)
%     xlabel('frequency (Hz)')
%     ylabel('a.u.')
    ylim([-0.4 0.3])
    title([sites{cmbs(ii,1)},'-',sites{cmbs(ii,2)}])
%     set(gca,'xticklabel',[],'yticklabel',[],'xcolor','k','ycolor','k',...
%         'linewidth',0.75)
end
% for ii = 1:6
%     subplot(4,3,6+ii)
%     hold on
% %     shadedErrorBar(1:100,decimate(lsdSalCohDiff(cInds(ii),:),5),...
% %         decimate(lsdSalCohDiffS(cInds(ii),:),5),'k',1)
%     plot(1:100,decimate(lsdSalCohDiff(cInds(ii),:),5),'k')
%     set(gca,'xtick',0:10:100)
%     ylim([-0.1 0.05])
%     plot([0 100],[0 0],'--k')
% end
%%
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
%% plot violins of feature changes
% load('lsd-baseVsal-base_zscore_stimImag_all-216feat.mat')
% subInds = [1:12,37:54,79:90,115:126,211:216];
subInds = 1:216;
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'})';
feat = feat(subInds);
allLSD = cat(1,thisLSDDiff{:});
allSAL = cat(1,thisBaseDiff{:});

these = cell(1);
for ii = subInds
    these = {these{:},allSAL(:,ii),allLSD(:,ii)};
end
%% significance
x = [];
ID = [];
group = [];
for ii = 1:5
    x = [x;thisBaseDiff{ii};thisLSDDiff{ii}];
    ID = [ID;repmat((ii-1)*2+1,size(thisBaseDiff{ii},1),1);...
        repmat(ii*2,size(thisLSDDiff{ii},1),1)];
    group = [group;repmat(0,size(thisBaseDiff{ii},1),1);...
        repmat(1,size(thisLSDDiff{ii},1),1)];
end
data = table(x(:,2),group,ID);
fitglme(data,'group ~ Var1+ID','Distribution','Binomial');
%% band pow
figure
for jj = 1:8
subplot(4,2,jj)
h = violin(these(:,(jj-1)*12+2:(jj-1)*12+13),'medc',[]);
for ii = 1:12
    if iseven(ii)
        h(ii).EdgeColor = 'r';
        h(ii).FaceColor = '#ffd2ed';
    else
        h(ii).EdgeColor = 'b';
        h(ii).FaceColor = '#e0dfff';
    end
    h(ii).FaceAlpha = 1;
    h(ii).LineWidth = 1;
end
legend off
ylim([-400 400])
end
%% band coh 
figure
for jj = 1:28
subplot(4,7,jj)
h = violin(these(:,(jj-1)*12+98:(jj-1)*12+109),'medc',[]);
for ii = 1:12
    if iseven(ii)
        h(ii).EdgeColor = 'r';
        h(ii).FaceColor = '#ffd2ed';
    else
        h(ii).EdgeColor = 'b';
        h(ii).FaceColor = '#e0dfff';
    end
    h(ii).FaceAlpha = 1;
    h(ii).LineWidth = 1;
end
title(feat((jj-1)*6+49))
ylim([-1.25 1.25])
legend off
end
%% power and coherence changes in LSD+stim relative to SAL+stim
salFiles = fileSearch('H:\Shared drives\dwielDoucetteLab\data\dualSite\processed\toUseSingleImag\','IL');
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
%%
x = []; y = []; a = [];
xS = []; yS = []; aS = [];
xP = []; yP = []; aP = [];
xSP = []; ySP = []; aSP = [];
betaStim = [];
% Subset indices to match other electrode array
% subInds = [1:12,37:54,79:90,115:126,211:216];
% subInds = [5,10:12,29,40,52,58,94,99,192,196];
subInds = 1:216;
% c = 1;
% for ii = 1:5
%     for jj = 1:5
%         cmbs(c,:) = [ii,jj];
%         c = c+1;
%     end
% end
% for ii = 1:size(cmbs,1)
%     thisBase = []; thisLSD = [];
%     for jj = 1:size(thisBaseDiff,2)
%         thisBase{jj} = thisBaseDiff{jj}(randperm(size(thisBaseDiff{jj},1),minSamp),:);
%         thisLSD{jj} = thisLSDDiff{jj}(randperm(size(thisLSDDiff{jj},1),minSamp),:);
%     end
%     otherLSD = logicFind(1,~ismember(1:5,cmbs(ii,1)),'==');
%     otherSAL = logicFind(1,~ismember(1:5,cmbs(ii,2)),'==');
%     testX = cat(1,thisBase{cmbs(ii,2)},thisLSD{cmbs(ii,1)});
%     testY = cat(1,zeros(minSamp,1),ones(minSamp,1));
%     trainX = cat(1,thisBase{otherSAL},thisLSD{otherLSD});
%     trainY = cat(1,zeros(minSamp*4,1),ones(minSamp*4,1));
%     mdl = fitglm(trainX(:,subInds),trainY,'distribution','binomial');
%     prob = predict(mdl,testX(:,subInds));
%     [~,~,~,aLOO(ii)] = perfcurve(testY,prob,1);
% end
%%
for n = 1:100
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
        % Permuted
        mdl = fitglm(trainX(:,subInds),trainY(randperm(numel(trainY),numel(trainY)),:),'distribution','binomial');
        prob = predict(mdl,testX(:,subInds));
        [xP{n},yP{n},~,aP(n)] = perfcurve(testY,prob,1);
        c = 1;
        for ii = subInds
            mdl = fitglm(trainX(:,ii),trainY,'distribution','binomial');
            betaStim(c,n) = table2array(mdl.Coefficients(2,1));
            prob = predict(mdl,testX(:,ii));
            [xS{n,c},yS{n,c},~,aS(n,c)] = perfcurve(testY,prob,1);
            [xSP{n,c},ySP{n,c},~,aSP(n,c)]= perfcurve(testY(randperm(numel(...
                testY),numel(testY))),prob,1);
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
sortA = sortA.*sign(sortBeta(sortInd));
% save(['G:\GreenLab\data\lsdStim\'...
%     'lsdStim-base_v_salStim-base_imag_all_216Feat.mat'],'x','y','a',...
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
%% single feature performance spreads
% subInds = [1:12,37:54,79:90,115:126,211:216];
% this = aS(:,subInds);
% thisP = aSP(:,subInds);
this = aS;
thisP = aSP;
figure
for ii = 1:8
    subplot(4,2,ii)
    hold on
    plotSpread(this(:,(ii-1)*6+1:ii*6),'xvalues',1-0.75:2*3:12*3-0.75,...
        'distributionColors','k','binWidth',0.005);
    for jj = 1:6
        plot([(jj-1)*6-1.25 (jj-1)*6+1.75],[mean(this(:,(ii-1)*6+jj)) ...
            mean(this(:,(ii-1)*6+jj))],'k','lineWidth',3)
        plot([(jj-1)*6+1.5 (jj-1)*6+4.25],[mean(thisP(:,(ii-1)*6+jj)) ...
            mean(thisP(:,(ii-1)*6+jj))],'color','#73898e','lineWidth',3)
    end
    plotSpread(thisP(:,(ii-1)*6+1:ii*6),'xvalues',2+0.75:2*3:12*3+0.75,...
        'distributionColors','#73898e','binWidth',0.005)
    ylim([0.3 0.8])
    title(feat((ii-1)*6+1))
%     set(gca,'xtick',1.5:2:12.5,'xticklabel',feat((ii-1)*6+1:ii*6))
end
%%
starts = (1:10:55)-2;
stops = (1:10:55)+2;
startsP = (5:10:55)-2;
stopsP = (5:10:55)+2;
figure
for ii = 1:28
    subplot(4,7,ii)
    hold on
    plotSpread(this(:,(ii-1)*6+1:ii*6),'xvalues',1:10:55,...
        'distributionColors','k','binWidth',0.005);
    plotSpread(thisP(:,(ii-1)*6+1:ii*6),'xvalues',5:10:55,...
        'distributionColors','#73898e','binWidth',0.005)
    for jj = 1:6
        plot([starts(jj) stops(jj)],[mean(this(:,(ii-1)*6+jj)) ...
            mean(this(:,(ii-1)*6+jj))],'k','lineWidth',3)
        plot([startsP(jj) stopsP(jj)],[mean(thisP(:,(ii-1)*6+jj)) ...
            mean(thisP(:,(ii-1)*6+jj))],'color','#73898e','lineWidth',3)
    end
%     set(gca,'xtick',1.5:2:12.5,'xticklabel',feat((ii-1)*6+1:(ii-1)*6))
end
%%
figure
hold on
plotSpread(this(:,1:24),'xvalues',1:2:48,'distributionColors','k')
plotSpread(thisP(:,1:24),'xvalues',2:2:48)
set(gca,'xtick',1.5:2:48.5,'xticklabel',featSub(1:24))
for ii = 1:24
        plot([(ii-1)*2+1-.25 (ii-1)*2+1+.25],[mean(this(:,ii)) mean(this(:,ii))],'b')
        plot([ii*2-.25 ii*2+.25],[mean(thisP(:,ii)) mean(thisP(:,ii))],'k')
end
figure 
hold on
plotSpread(this(:,25:60),'xvalues',1:2:72,'distributionColors','k')
plotSpread(thisP(:,25:60),'xvalues',2:2:72)
set(gca,'xtick',1.5:2:72.5,'xticklabel',featSub(25:60))
for ii = 1:36
        plot([(ii-1)*2+1-.25 (ii-1)*2+1+.25],[mean(this(:,ii+24)) ...
            mean(this(:,ii+24))],'b')
        plot([ii*2-.25 ii*2+.25],[mean(thisP(:,ii+24)) ...
            mean(thisP(:,ii+24))],'k')
end
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
% [data,samp1,files,time1] = collateData(['D:\dualSite\processed\'...
%     'toUseSingleImag\'],{'IL','in'},{'pow','coh'},'trl','');
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
feat = names({'lIL','rIL','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
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
plot(lsdMeanA(abs(diffMean)>=0.6),salMeanA(abs(diffMean)>=0.6),'ok')

inds = logicFind(1,abs(diffMean)>=0.6,'==');
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
%% New LSD + stim data
% load baselines
[baseData,baseSamps,baseFiles] = collateData('F:\LSD+stim\processed\baselines\',...
    {'.mat'},{'pow','coh'},'trl','rel');
% get average and std of each animal
for ii = 1:8
    these = contains(baseFiles{1},['LSD',num2str(ii)]);
    baseCat{ii} = cat(1,baseData{1}{these});
    baseM{ii} = mean(cat(1,baseData{1}{these}));
    baseS{ii} = std(cat(1,baseData{1}{these}));
end
% load last stim days
[lastData,lastSamps,lastFiles] = collateData('F:\LSD+stim\processed\lastStim\',...
    {'.mat'},{'pow','coh'},'trl','rel');
%%
% Find last stim epoch and split data into stim data and postStim data
for ii = 1:numel(lastFiles{1})
    load(['F:\LSD+stim\processed\lastStim\',lastFiles{1}{ii}],'hist')
    % Get animal ID
    parts = strsplit(lastFiles{1}{ii},'-');
    ID = str2double(parts{1}(end));
    % Get stim times
    if contains(lastFiles{1}{ii},'sham')
        stimStart = 600;
        stimStop = 6100;
    else
        stimStart = hist.eventTs.t{9}(nearest_idx3(600,hist.eventTs.t{9}));
        stimStop = hist.eventTs.t{9}(end);
    end
    % Stim data
    stimData{ii} = lastData{1}{ii}(...
        lastSamps{1}{ii}(:,1)>(stimStart*2000) & ...
        lastSamps{1}{ii}(:,2)<(stimStop*2000),:);%-baseM{ID})./baseS{ID};
    % Post stim data
    postStimData{ii} = lastData{1}{ii}(...
        lastSamps{1}{ii}(:,2)>(stimStop*2000),:);%-baseM{ID})./baseS{ID};
end
%% First build models differentiating between stim and baseline
% code in runLSDstim.mat

for ii = 1:100
    clear lsdStimA lsdStimAP lsdSham salStimA salStimAP salShamA salShamAP
    load(['F:\LSD+stim\LSDvSAL_acute\LSDstim',num2str(ii),'.mat'])
    lsdStimAs(ii) = lsdStimA{1}.acc;
    lsdStimAPs(ii) = lsdStimAP{1}.acc;
    lsdShamAs(ii) = lsdShamA{1}.acc;
    lsdShamAPs(ii) = lsdShamAP{1}.acc;
    if exist('salStimA')
        salStimAs(ii) = salStimA{1}.acc;
    else
        salStimAs(ii) = NaN;
    end
    if exist('salStimAP')
        salStimAPs(ii) = salStimAP{1}.acc;
    else
        salStimAPs(ii) = NaN;
    end
    if exist('salShamA')
        salShamAs(ii) = salShamA{1}.acc;
    else
        salShamAs(ii) = NaN;
    end
    if exist('salShamAP')
        salShamAPs(ii) = salShamAP{1}.acc;
    else
        salShamAPs(ii) = NaN;
    end
end
%%
figure
subplot(2,2,1)
[f,xi,bw] = ksdensity(lsdStimAPs);
fill(xi,f*bw,'w')
plotSpread(lsdStimAPs','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(lsdStimAs) mean(lsdStimAs)],[0 0.25],'-')
subplot(2,2,2)
[f,xi,bw] = ksdensity(lsdShamAPs);
fill(xi,f*bw,'w')
plotSpread(lsdShamAPs','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(lsdShamAs) mean(lsdShamAs)],[0 0.25],'-')
subplot(2,2,3)
[f,xi,bw] = ksdensity(salStimAPs);
fill(xi,f*bw,'w')
plotSpread(salStimAPs','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(salStimAs) mean(salStimAs)],[0 0.25],'-')
subplot(2,2,4)
[f,xi,bw] = ksdensity(salShamAPs);
fill(xi,f*bw,'w')
plotSpread(salShamAPs','xyori','flipped','xvalues',0.05,'spreadWidth',0.1)
plot([mean(salShamAs) mean(salShamAs)],[0 0.25],'-')
%%
thisDiff = [];
for ii = 1:60
    thisDiff = [thisDiff;lsdStimAPs(ii)-lsdShamAPs(ii);...
        lsdStimAPs(ii)-salStimAPs(ii);...
        lsdStimAPs(ii)-salShamAPs(ii);...
        lsdShamAPs(ii)-salStimAPs(ii);...
        lsdShamAPs(ii)-salShamAPs(ii);...
        salStimAPs(ii)-salShamAPs(ii)];
end
%% Then use post stim effects
% code in runLSDpostSim.mat
[lsdPostStimAs,lsdPostStimAPs,lsdPostShamAs,lsdPostShamAPs,...
    salPostStimAs,salPostStimAPs,salPostShamAs,salPostShamAPs] = deal([]);
for ii = 1:60
    load(['F:\LSD+stim\LSDvsSAL_post\LSDpostStim',num2str(ii),'.mat'])
    lsdPostStimAs(ii) = lsdStimPostA{1}.acc;
    lsdPostStimAPs(ii) = lsdStimPostAP{1}.acc;
    lsdPostShamAs(ii) = lsdShamPostA{1}.acc;
    lsdPostShamAPs(ii) = lsdShamPostAP{1}.acc;
    salPostStimAs(ii) = salStimPostA{1}.acc;
    salPostStimAPs(ii) = salStimPostAP{1}.acc;
    salPostShamAs(ii) = salShamPostA{1}.acc;
    salPostShamAPs(ii) = salShamPostAP{1}.acc;
end
%% Then use baseline normalized post days in X minute bins
[postData,postSamps,postFiles] = collateData(...
    'F:\LSD+stim\processed\postStim\',{'.mat'},{'pow','coh'},'trl','');
% Number of minutes per bin
bin = 5;
binSamps = bin*60*2000;
% Total length of data to look at (minutes)
total = 30;
totalSamps = total*60*2000;
for ii = 1:numel(postFiles{1})
    for k = 1:total/bin
        thisStart = (k-1)*binSamps+1;
        thisStop = k*binSamps;
        inds = postSamps{1}{ii}(:,1)>=thisStart & ...
            postSamps{1}{ii}(:,2)<=thisStop;
        post{ii,k} = postData{1}{ii}(inds,:);
    end
end
%% Build models for each bin
stimInd = logical([1,0,0,1,1,1,0,0,1,0,0,0,0]);
shamInd = ~stimInd;
lsdInd = logical([1,0,1,0,0,1,0,1,1,0,1,0,1]);
salInd = ~lsdInd;
n = 50;
for ii = 1:100
    for k = 1:size(post,2)
        for jj = 1:size(post,1)
            % Grab n samples, or weight as needed
            if size(post{jj,k},1) < n
                thisPost{jj,k} = post{jj,k};
                thisWeight{jj,k} = repmat(n/size(post{jj,k},1),...
                    size(post{jj,k},1),1);
            else
                thisPost{jj,k} = post{jj,k}(randperm(size(post{jj,k},1),...
                    n),:);
                thisWeight{jj,k} = ones(n,1);
            end
        end
        % Combine into stim and sham groups
        allStim = cat(1,thisPost{stimInd,k});
        allSham = cat(1,thisPost{shamInd,k});
        allW = cat(1,thisWeight{stimInd,k},thisWeight{shamInd,k});
        % 80:20 split
        [trainX,trainY,testX,testY,trainInd] = trainTest([allStim;allSham],...
            [ones(size(allStim,1),1);zeros(size(allSham,1),1)],0.2);
        cfg = lassoNetCfg({testX,testY},[],'n','n','n',100,'1se',...
            allW(trainInd));
        [~,~,~,~,accPostStim{ii,k},~] = lassoNet(trainX,trainY,'binomial',...
            'class',1,10,1,cfg);
        % Permuted
        stimIndPerm = randperm(numel(postFiles{1}),sum(stimInd));
        shamIndPerm = randperm(numel(postFiles{1}),sum(shamInd));
        allStimPerm = cat(1,thisPost{stimIndPerm,k});
        allShamPerm = cat(1,thisPost{shamIndPerm,k});
        allWPerm = cat(1,thisWeight{stimIndPerm,k},...
            thisWeight{shamIndPerm,k});
        [trainX,trainY,testX,testY,trainInd] = trainTest([allStimPerm;...
            allShamPerm],[ones(size(allStimPerm,1),1);...
            zeros(size(allShamPerm,1),1)],0.2);
        cfg = lassoNetCfg({testX,testY},[],'n','n','n',100,'1se',[]);
        [~,~,~,~,accPostStimP{ii,k},~] = lassoNet(trainX,trainY,...
            'binomial','class',1,10,1,cfg);
    end
end
%% Models of every stim day, using that day's baseline and post stim
[stimData,stimSamps,stimFiles] = collateData(...
    'F:\LSD+stim\processed\allStim\',{'.mat'},{'pow','coh'},'trl','rel');
%%
[lsdSham,lsdStim,salSham,salStim] = deal(cell(1,2));
for ii = 1:numel(stimFiles{1})
    load(stimFiles{1}{ii},'hist')
    if contains(stimFiles{1}{ii},'sham')
        firstStim = 600;
        lastStim = 6100;
    else
        firstStim = hist.eventTs.t{9}(nearest_idx3(600,hist.eventTs.t{9}));
        lastStim = hist.eventTs.t{9}(end);
    end
    thisPre = stimData{1}{ii}(stimSamps{1}{ii}(:,1)<firstStim*2000,:);
    thisPost = stimData{1}{ii}(stimSamps{1}{ii}(:,2)>lastStim*2000,:);
    if contains(stimFiles{1}{ii},'sham') && contains(stimFiles{1}{ii},'_LSD_')
        lsdSham = [lsdSham;{thisPre},{thisPost}];
    end
    if contains(stimFiles{1}{ii},'sham') && contains(stimFiles{1}{ii},'_SAL_')
        salSham = [salSham;{thisPre},{thisPost}];
    end
    if contains(stimFiles{1}{ii},'sIL') && contains(stimFiles{1}{ii},'_LSD_')
        lsdStim = [lsdStim;{thisPre},{thisPost}];
    end
    if contains(stimFiles{1}{ii},'sIL') && contains(stimFiles{1}{ii},'_SAL_')
        salStim = [salStim;{thisPre},{thisPost}];
    end
end
lsdSham(1,:) = [];
lsdStim(1,:) = [];
salSham(1,:) = [];
salStim(1,:) = [];
%%
n = 50;
% LSD stim
for k = 1:100
    disp(k)
    this = [];
    for ii = 1:size(lsdStim,1)
        this{ii,1} = lsdStim{ii,1}(randperm(size(lsdStim{ii,1},1),n),:);
%         this{ii,2} = lsdStim{ii,2}(randperm(size(lsdStim{ii,2},1),n),:);
        this{ii,2} = lsdStim{ii,2}(1:50,:);
    end
    [trainX,trainY,testX,testY] = trainTest(cat(1,this{:,1},this{:,2}),...
        cat(1,ones(size(lsdStim,1)*n,1),zeros(size(lsdStim,1)*n,1)),0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~, lsdStimA(k)] = perfcurve(testY,prob,1);
    % Permuted
    these = logical(round(rand(size(lsdStim,1),1)));
    [trainX,trainY,testX,testY] = trainTest(cat(1,this{these,1},...
        this{~these,2},this{~these,1},this{these,2}),[ones(sum(these)*n*2,1);...
        zeros(sum(~these)*n*2,1)],0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,lsdStimAP(k)] = perfcurve(testY,prob,1);
    % LSD sham
    this = [];
    for ii = 1:size(lsdSham,1)
        this{ii,1} = lsdSham{ii,1}(randperm(size(lsdSham{ii,1},1),n),:);
        this{ii,2} = lsdSham{ii,2}(randperm(size(lsdSham{ii,2},1),n),:);
    end
    [trainX,trainY,testX,testY] = trainTest(cat(1,this{:,1},this{:,2}),...
        cat(1,ones(size(lsdSham,1)*n,1),zeros(size(lsdSham,1)*n,1)),0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,lsdShamA(k)] = perfcurve(testY,prob,1);
    % Permuted
    these = logical(round(rand(size(lsdSham,1),1)));
    [trainX,trainY,testX,testY] = trainTest(cat(1,this{these,1},...
        this{~these,2},this{~these,1},this{these,2}),[ones(sum(these)*n*2,1);...
        zeros(sum(~these)*n*2,1)],0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,lsdShamAP(k)] = perfcurve(testY,prob,1);
    % SAL stim
    this = [];
    for ii = 1:size(salStim,1)
        this{ii,1} = salStim{ii,1}(randperm(size(salStim{ii,1},1),n),:);
        this{ii,2} = salStim{ii,2}(randperm(size(salStim{ii,2},1),n),:);
    end
    [trainX,trainY,testX,testY] = trainTest(cat(1,this{:,1},this{:,2}),...
        cat(1,ones(size(salStim,1)*n,1),zeros(size(salStim,1)*n,1)),0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,salStimA(k)] = perfcurve(testY,prob,1);
    % Permuted
    these = logical(round(rand(size(salStim,1),1)));
    [trainX,trainY,testX,testY] = trainTest(cat(1,this{these,1},...
        this{~these,2},this{~these,1},this{these,2}),[ones(sum(these)*n*2,1);...
        zeros(sum(~these)*n*2,1)],0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,salStimAP(k)] = perfcurve(testY,prob,1);
    % SAL sham
    this = [];
    for ii = 1:size(salSham,1)
        this{ii,1} = salSham{ii,1}(randperm(size(salSham{ii,1},1),n),:);
        this{ii,2} = salSham{ii,2}(randperm(size(salSham{ii,2},1),n),:);
    end
    [trainX,trainY,testX,testY] = trainTest(cat(1,this{:,1},this{:,2}),...
        cat(1,ones(size(salSham,1)*n,1),zeros(size(salSham,1)*n,1)),0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,salShamA(k)] = perfcurve(testY,prob,1);
    % Permuted
    these = logical(round(rand(size(salSham,1),1)));
    [trainX,trainY,testX,testY] = trainTest(cat(1,this{these,1},...
        this{~these,2},this{~these,1},this{these,2}),[ones(sum(these)*n*2,1);...
        zeros(sum(~these)*n*2,1)],0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [~,~,~,salShamAP(k)] = perfcurve(testY,prob,1);
end
%% Build model of LSD+stim vs base; apply over post
