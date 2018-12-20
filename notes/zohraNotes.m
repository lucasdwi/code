data = [reshape(pow{1}.relPow,48,size(pow{1}.Pow,3));reshape(coh{1}.normBandCoh,28*6,size(pow{1}.Pow,3))]';
for ii = 3:98
    if ~isempty(pow{ii})
        this = [reshape(pow{ii}.relPow(:,:,1),48,1);reshape(coh{ii}.normBandCoh(:,:,1),28*6,1)]';
        data = [data;this];
    end
end
x{1} = data(1:size(pow{1}.Pow,3),:);
x{2} = data(size(pow{1}.Pow,3)+1:end,:);
[all,~,rnd,~] = evenDataSplit(x,0.8,0.2,'',100);
%%
for ii = 1:100
   load(['C:\Users\Pythia\Documents\GreenLab\data\Zohra\model\150v250_noImput_v2\ZA1_150v250_',num2str(ii),'.mat']) 
   A(ii) = a;
   X(:,ii) = x;
   Y(:,ii) = y;
   rA(ii) = rndA;
   rX{ii} = rndX;
   rY{ii} = rndY;
%    betas(:,:,ii) = model.glmnet_fit.beta;
end
clear rndX rndY
[rndX,rndY] = randRoc(rX,rY);
figure
hold on
plot(mean(X,2),mean(Y,2),'-k')
plot(mean(rndX,2),mean(rndY,2),'--','Color',[0.5 0.5 0.5])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
legend({['Real: ',num2str(round(mean(A),3)),'\pm',num2str(round(std(A),3))],['Permuted: ',num2str(round(mean(rA),2)),'\pm',num2str(round(std(rA),2))]},'Location','se')
%%
for ii = 1:100
    load(['C:\Users\Pythia\Documents\GreenLab\data\Zohra\model\150v250_noImput_log_v3\ZA1_250v150_log_',num2str(ii),'.mat'])
    allA(ii,:) = a;
    allOR(ii,:) = or;
end
[sortMA,sInd] = sort(mean(allA,1),'descend');
sortSA = std(allA(:,sInd),[],1);
sortOR = mean(allOR(:,sInd),1);
% scatterErr(1:216,sortMA,sortSA,1,'line','-');
%% Compare 150 and 250
load('C:\Users\Pythia\Documents\GreenLab\data\Zohra\analyzed\ZA1_NAc150Burst_2018_05_02_analyzed.mat')
pow150 = [];
powZ150 = [];
for ii = 3:98
    if ~isempty(pow{ii})
        pow150(:,:,ii) = pow{ii}.Overall(:,:,1);
        powZ150(:,:,ii) = (pow{ii}.Overall(:,:,1)-pow{1}.Overall)./pow{1}.OverallStd;
    end
end
base150Pow = pow{1}.Overall;

load('C:\Users\Pythia\Documents\GreenLab\data\Zohra\analyzed\ZA1_NAcBurst250A_2018_04_16_analyzed.mat')
pow250 = [];
powZ250 = [];
for ii = 3:98
    if ~isempty(pow{ii})
        pow250(:,:,ii) = pow{ii}.Overall(:,:,1);
        powZ250(:,:,ii) = (pow{ii}.Overall(:,:,1)-pow{1}.Overall)./pow{1}.OverallStd;
    end
end
base250Pow = pow{1}.Overall;
%% Plot
chan = 2;
names = {'mPFCr','mPFCl','OFCr','OFCl','NAcr','NAcl','LHr','LHl'};
figure
plot(pow{1}.f,mean(pow150(chan,:,:),3))
hold on
plot(pow{1}.f,mean(pow250(chan,:,:),3))
legend({'150A','250A'})
title([names{chan}, ' Power'])
xlabel('Frequency (Hz)') 
ylabel('Power (dB)')

figure
hold on
plot(pow{1}.f,mean(powZ150(chan,:,:),3))
plot(pow{1}.f,mean(powZ250(chan,:,:),3))
legend({'150A','250A'},'Location','se')
title([names{chan},' Power'])
xlabel('Frequency (Hz)') 
ylabel('Power (z-score)')
%% Indices rank from betas
for jj = 1:size(betas,3)
    for ii = 1:size(betas,1)
        thisInd = logicFind(0,betas(ii,:,jj),'~=','first');
        if ~isempty(thisInd)
            ind(ii,jj) = thisInd;
        end
    end
end
ind(ind==0) = NaN;
for ii = 1:size(ind,2)
    [b(ii,:),sInd(ii,:)]=sort(ind(:,ii));
end
nameVect = names({'mPFCr','mPFCl','OFCr','OFCl','NAcr','NAcl','LHr','LHl'},{'d','t','a','b','lg','hg'});
topFeat = nameVect(sInd);
ranks = [];
for ii = 1:size(sInd,1)
    for jj = 1:size(sInd,2)
        if isnan(ind(sInd(ii,jj),ii))
            ranks(ii,sInd(ii,jj)) = NaN;
        else
            ranks(ii,sInd(ii,jj)) = jj;
        end
    end
end
mRank = mean(ranks,1,'omitnan');
%% Plot average power and coherence for NAc L of 250A comparing interstim 
% intervals and baseline
load('ZA1_NAcBurst250A_2018_04_16_analyzed.mat')
% load('ZA1_NAcBurst_2018_04_11_analyzed.mat')
interP = []; interC = [];
interBandP = []; interBandC = [];
for ii = 3:98
    if ~isempty(pow{ii})
        interP = cat(3,interP,pow{ii}.Pow);
%         interP = cat(3,interP,zeros(8,100).*NaN);
        interC = cat(3,interC,coh{ii}.Cxy);
%         interC = cat(3,interC,zeros(28,500).*NaN);
        
        interBandP = cat(3,interBandP,pow{ii}.bandPow);
%         interBandP = cat(3,interBandP,zeros(6,8).*NaN);
        interBandC = cat(3,interBandC,coh{ii}.mBandCoh);
%         interBandC = cat(3,interBandC,zeros(28,6).*NaN);
    end
end
%%
% NAc L
figure
hold on
h(1) = shadedErrorBar(1:100,pow{1}.Overall(6,:),pow{1}.OverallStd(6,:),'-k');
% h(2) = shadedErrorBar(1:100,mean(interP(6,:,:),3,'omitnan'),nanstd(interP(6,:,:),[],3),'-b',1);
title('Left NAc Power')
legend([h(1).mainLine,h(2).mainLine],{'Base','Inter'})
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
% NAc L - mPFC L
figure
hold on
h(1) = shadedErrorBar(coh{1}.f,mean(coh{1}.Cxy(11,:,:),3),std(coh{1}.Cxy(11,:,:),[],3),'-k');
% h(2) = shadedErrorBar(coh{1}.f,mean(interC(11,:,:),3,'omitnan'),nanstd(interC(11,:,:),[],3),'-b',1);
title('Left NAc - Left mPFC Coherence')
% legend([h(1).mainLine,h(2).mainLine],{'Base','Inter'},'location','se')
xlabel('Frequency (Hz)')
ylabel('Magnitude Squared Coherence')
%% Plot all bands in z-score space
zInterBandP = (interBandP-mean(pow{1}.bandPow,3))./std(pow{1}.bandPow,[],3);
zInterBandC = (interBandC-mean(coh{1}.mBandCoh,3))./std(coh{1}.mBandCoh,[],3);
feats = [reshape(zInterBandP,48,size(zInterBandP,3));reshape(zInterBandC,168,size(zInterBandC,3))];
figure
h = pcolor(feats);
set(h,'edgecolor','none');
colormap('viridis')
xlabel('Interstim Interval')
ylabel('Features')
h = colorbar;
ylabel(h,'Z-Score')
%% Plot all features (frequencies) in z-score space
zInterP = (interP-pow{1}.Overall)./pow{1}.OverallStd;
zInterC = (interC-mean(coh{1}.Cxy,3))./std(coh{1}.Cxy,[],3);

figure
h = pcolor(squeeze(zInterP(6,:,:)));
set(h,'edgecolor','none');
colormap('viridis')
xlabel('Interstim Interval')
ylabel('Frequency (Hz)')
title('Interstim Left NAcPower')
h = colorbar;
ylabel(h,'Z-Score')


figure
h = pcolor(squeeze(zInterC(4,:,:)));
set(h,'edgecolor','none','YData',coh{1}.f);
ylim([0 100])
colormap('viridis')
xlabel('Interstim Interval')
ylabel('Frequency (Hz)')
title('Interstim Left NAc - Left mPFC Coherence')
h = colorbar;
ylabel(h,'Z-Score')
%% 3 Hz Stim
load('ZA1_StimNAc3Hz_2018_04_18_base_vs_post.mat')
powZ = (psdTrls{2}.Pow-psdTrls{1}.Overall)./psdTrls{1}.OverallStd;
cohZ = (coh{2}.Cxy-mean(coh{1}.Cxy,3))./std(coh{1}.Cxy,[],3);
%%
figure
for ii = 1:8
    subplot(2,4,ii)
    imagesc(squeeze(powZ(ii,:,1:36)))
    caxis([-3 5])
    colormap viridis
end
%%
figure
c = 1;
for ii = [4,5,10,11,15,16,19,20,25]
   subplot(3,3,c)
   plot(coh{1}.f(1:100),mean(coh{1}.Cxy(ii,1:100,:),3))
   hold on
   plot(coh{1}.f(1:100),mean(coh{2}.Cxy(ii,1:100,:),3))   
   yLim = ylim;
   plot([3 3],[0 1],':k')
   ylim(yLim)
   c = c+1;
end
tightfig
%%
figure
for ii = 1:28
    subplot(4,7,ii)
    imagesc(squeeze(cohZ(ii,:,:)))
    colormap viridis
end
%%
pre = [reshape(pow{1}.bandPow,48,182);reshape(coh{1}.mBandCoh,168,182)];
inter = [reshape(interBandP,48,96);reshape(interBandC,168,96)];
for ii = 1:size(pre,1)
   p(ii) = ranksum(pre(ii,:),inter(ii,:)); 
end
%%
fileName = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\Zohra\','.mat','in');
for ii = 1:size(fileName,2)
    load(fileName{ii})
    
    LFPTs.data = filter60(LFPTs.data,adfreq,0);
    [LFPTspre,chk_nan,zeroedChannel] = threshFilt(LFPTs,0.75,0.125,0.5,...
        5,adfreq);
    preS = []; postS = [];
    for jj = 1:8
        [preS{ii}(:,:,jj),f,t] = spectrogram(LFPTspre.data(jj,1:nearest_idx3(eventTs.t{1}(1),LFPTspre.tvec)),1024,512,1:100,adfreq);
        [postS{ii}(:,:,jj),f,t] = spectrogram(LFPTs.data(jj,nearest_idx3(eventTs.t{2}(end),LFPTs.tvec):end),1024,512,1:100,adfreq);
    end
    
    preP{ii} = 10*log10(abs(preS{ii}).^2);
    mPreP{ii} = mean(preP{ii},2,'omitnan');
    sPreP{ii} = nanstd(preP{ii},[],2);
    postP{ii} = 10*log10(abs(postS{ii}).^2);
    zPostP{ii} = (postP{ii}-repmat(mPreP{ii},1,size(postP{ii},2),1))./repmat(sPreP{ii},1,size(postP{ii},2),1);
    percPostP{ii} = (postP{ii}-repmat(mPreP{ii},1,size(postP{ii},2),1))./repmat(mPreP{ii},1,size(postP{ii},2),1);
    
    zDeltaP(ii,:,:) = mean(zPostP{ii}(1:4,1:100,:));
    percDeltaP(ii,:,:) = mean(percPostP{ii}(1:4,1:100,:));
    zThetaP(ii,:,:) = mean(zPostP{ii}(5:10,1:100,:));
    percThetaP(ii,:,:) = mean(percPostP{ii}(5:10,1:100,:));
end
labels = {'OFC-R','OFC-L','mPFC-R','mPFC-L','NAc-R','NAc-L','LH-R','LH-L'};
%%
for jj = 1:4
    figure
    for ii = 1:8
        subplot(2,4,ii)
        h = imagesc(zPostP{jj}(:,1:30,ii));
        set(gca,'XTickLabel',round(t([10,20,30])),'ydir','normal')
        if ii == 5
            xlabel('Time (s)')
            ylabel('Frequency (Hz)')
        end
        colormap viridis
        caxis([-3 3])
    end
end
%%
for jj = 1:4
    figure
    for ii = 1:8
        subplot(2,4,ii)
        h = imagesc(percPostP{jj}(:,:,ii).*100);
%         set(gca,'XTickLabel',round(t([10,20,30])),'ydir','normal')
        if ii == 5
            xlabel('Time (s)')
            ylabel('Frequency (Hz)')
        end
        colormap viridis
        caxis([-200 200])
    end
suptitle(fileName{jj})    
end

%% Line Delta
figure
for ii = 1:4
    for jj = 1:8
        subplot(2,4,jj)
        hold on
        %         h = plot(t(1:30),zDeltaP(ii,1:30,jj),'linewidth',1);
        h = plot(t(1:30),percDeltaP(ii,1:30,jj),'linewidth',1);
        if ii == 1
            set(h,'Color','b','linestyle','--');
        elseif ii == 2
            set(h,'Color','r','linestyle','--');
        elseif ii == 3
            set(h,'color','b');
        elseif ii == 4
            set(h,'color','r');
        end
%         plot([0 15],[1.96 1.96],':k','linewidth',1)
%         plot([0 15],[-1.96 -1.96],':k','linewidth',1)
        title(labels{jj})
        ylim([-3.5 2])
        xlim([0 15])
    end
end
tightfig
suptitle('Delta')
%% Bar Delta
figure
bar(squeeze([zDeltaP(1,1,:);zDeltaP(3,1,:);zDeltaP(2,1,:);zDeltaP(4,1,:)])')
set(gca,'XTickLabel',labels,'xticklabelrotation',45)
legend({'Burst250','Cont250','Burst150','Cont150'},'location','se')
box off
title('Delta: First Sample')
%% Line Theta
figure
for ii = 1:4
    for jj = 1:8
        subplot(2,4,jj)
        hold on
        h = plot(t(1:30),percThetaP(ii,1:30,jj),'linewidth',1);
        if ii == 1
            set(h,'Color','b','linestyle','--');
        elseif ii == 2
            set(h,'Color','r','linestyle','--');
        elseif ii == 3
            set(h,'color','b');
        elseif ii == 4
            set(h,'color','r');
        end
%         plot([0 15],[1.96 1.96],':k','linewidth',1)
%         plot([0 15],[-1.96 -1.96],':k','linewidth',1)
        title(labels{jj})
%         ylim([-3 2])
        xlim([0 15])
    end
end
tightfig
suptitle('Theta')
%% Bar Theta
figure
bar(squeeze([zThetaP(1,1,:);zThetaP(3,1,:);zThetaP(2,1,:);zThetaP(4,1,:)])')
set(gca,'XTickLabel',labels,'xticklabelrotation',45)
legend({'Burst250','Cont250','Burst150','Cont150'},'location','se')
box off
title('Theta: First Sample')
%% Get mean zPostP for delta and theta
mDeltaPreP = mean(sum(preP(1:4,:,:)),2,'omitnan');
sDetlaPreP = nanstd(sum(preP(1:4,:,:)),[],2);
mThetaPreP = mean(sum(preP(5:10,:,:)),2,'omitnan');
sThetaPreP = nanstd(sum(preP(5:10,:,:)),[],2);

deltaPostP = sum(postP(1:4,:,:));
thetaPostP = sum(postP(5:10,:,:));
zDelta = (deltaPostP-repmat(mDeltaPreP,1,size(deltaPostP,2),1))./repmat(mDeltaPreP,1,size(deltaPostP,2),1);

%% Grab baseline and post - set eventTs
start = nearest_idx3(eventTs.t{1,1}(1),LFPTs.tvec);
eventTs.t{1,3} = LFPTs.tvec(1);
eventTs.t{1,4} = LFPTs.tvec(start);
eventTs.label{1,3} = 'Base (Start)';
eventTs.label{1,4} = 'Base (End)';

stop = nearest_idx3(eventTs.t{1}(end),LFPTs.tvec);
eventTs.t{1,5} = LFPTs.tvec(stop);
eventTs.t{1,6} = LFPTs.tvec(end);
eventTs.label{1,5} = 'Post (Start)';
eventTs.label{1,6} = 'Post (End)';

% Separate stim and clean data
stimInd = []; stim = []; stimStop = []; stimStart = [];
% stimInd = logicFind(1,~cellfun(@isempty,eventTs.t),'==');
stim = eventTs.t{1,13};
% Get all stops
stimStop = stim(logicFind(0,round(diff(stim)),'~='));
% Add last stop since above algorithm skips
stimStop = [stimStop;stim(end)];
% Go 0.5 seconds before each stop
stimStop = stimStop;
% Get starts using the next index after stops from stim
stimStart = stim(logicFind(0,round(diff(stim)),'~=')+1);
% Add first start since above algorithm skips
stimStart = [stim(1);stimStart];
% Use stim start and stops to determine data start and stops; also include
% first/last data point - assumes that no stim starts or ends a recording
dataStart = [LFPTs.tvec(1);LFPTs.tvec(nearest_idx3(stimStop,LFPTs.tvec)+1)'];
dataStop = [LFPTs.tvec(nearest_idx3(stimStart,LFPTs.tvec)-1)';LFPTs.tvec(end)];
struct.t = [];
struct.label = [];
for ii = 2:size(stimStart,1)
    struct.t = [struct.t,{dataStart(ii)},{dataStop(ii)}];
    struct.label = [struct.label,{['Int',num2str(ii-1),' (Start)']},{['Int',num2str(ii-1),' (End)']}];
end
struct.t = [{dataStart(1)},{dataStop(1)},struct.t,{LFPTs.tvec(stop)},{LFPTs.tvec(end)}];
struct.label = [{'Base (Start)'},{'Base (End)'},struct.label,{'Post (Start)'},{'Post (End)'}];
oldEventTs = eventTs;
eventTs = struct;
eoi = {'base',[0 5];'post',[0 5]};
for ii = 1:size(stimStart,1)-1
    eoi = [eoi;{['Int',num2str(ii)],[0 5]}];
end
%%
[LFPTs,~,zeroedChannel,~,trls,adfreq] = preProcess(LFPTs,adfreq,1,2,0.125,0.5,eoi,eventTs);
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];...
    'lgamma',[45,65];'hgamma',[70,90]};
[pow,~] = powerComp(trls,adfreq,bands,[55 65],[1 1 100],eoi,'n');
[coh,~] = cohComp(trls,adfreq,eoi,bands,zeroedChannel,[1 1 100],...
    [55 65],'n','mtm','NW',8);
%%
fileName = 'ZA1_Stim22250_2018_04_16.mat';
load(fileName)
eventTs.t(:,[1:14,17:end]) = [];
eventTs.label(:,[1:14,17:end]) = [];
LFPTs.data = filter60(LFPTs.data,adfreq,0);
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,1,0.125,0.5,...
    5,adfreq);
%% Grab baseline and post - set eventTs
start = nearest_idx3(eventTs.t{1,1}(1),LFPTs.tvec);
eventTs.t{1,3} = LFPTs.tvec(1);
eventTs.t{1,4} = LFPTs.tvec(start);
eventTs.label{1,3} = 'Base (Start)';
eventTs.label{1,4} = 'Base (End)';

stop = nearest_idx3(eventTs.t{2}(end),LFPTs.tvec);
eventTs.t{1,5} = LFPTs.tvec(stop);
eventTs.t{1,6} = LFPTs.tvec(end);
eventTs.label{1,5} = 'Post (Start)';
eventTs.label{1,6} = 'Post (End)';
%% Calculate power and coherence
eoi = {'base',[0 5];'post',[0 5]};

LFPTs.data(isnan(LFPTs.data)) = eps;
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];...
    'lgamma',[45,65];'hgamma',[70,90]};
[pow,powerPlots] = powerComp(trls,adfreq,bands,'y',[1 1 100],eoi,'n');
[coh,cohPlots] = cohComp(trls,adfreq,eoi,bands,zeroedChannel,[1 1 100],...
    'y','n','mtm','NW',8);
%% Reshape raw power and coherence
prePow = pow{1}.Pow;
prePow(abs(prePow)>200) = NaN;
prePowM = mean(prePow,3,'omitnan');
prePowS = nanstd(prePow,[],3);
postPow = pow{2}.Pow;
postPow(abs(postPow)>200) = NaN;
zPostPow = (postPow - repmat(prePowM,1,1,size(postPow,3)))./repmat(prePowS,1,1,size(postPow,3));

preCohM = mean(coh{1}.Cxy,3,'omitnan');
preCohS = nanstd(coh{1}.Cxy,[],3);
postCoh = coh{2}.Cxy;
zPostCoh = (postCoh - repmat(preCohM,1,1,size(postCoh,3)))./repmat(preCohS,1,1,size(postCoh,3));
%%
figure
for ii = 1:8
    subplot(2,4,ii)
    imagesc(squeeze(zPostPow(ii,:,1:10)))
    set(gca,'ydir','normal')
    caxis([-3 3])
    colormap viridis
end
%%
figure
c = 1;
for ii = [4,5,10,11,15,16,19,20,23]
    subplot(2,5,c)
    imagesc(squeeze(zPostCoh(ii,1:100,1:10)))
    caxis([-3 3])
    colormap viridis
    c = c+1;
end
%%
load('ZA1_StimNAc_2018_04_10_8chan.mat')
LFPTs.data = filter60(LFPTs.data,adfreq,0);
[LFPTs,adfreq] = dwnSample(LFPTs,4,adfreq);
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,0.5,0.125,1,...
    5,adfreq);
%% Grab baseline and post - set eventTs
start = nearest_idx3(eventTs.t{1,3}(1),LFPTs.tvec);
eventTs.t{1,5} = LFPTs.tvec(1);
eventTs.t{1,6} = LFPTs.tvec(start);
eventTs.label{1,5} = 'Base (Start)';
eventTs.label{1,6} = 'Base (End)';

stop = nearest_idx3(eventTs.t{1,4}(end),LFPTs.tvec);
eventTs.t{1,7} = LFPTs.tvec(stop);
eventTs.t{1,8} = LFPTs.tvec(end);
eventTs.label{1,7} = 'Post (Start)';
eventTs.label{1,8} = 'Post (End)';
%%
eoi = {'base',[0 5];'post',[0 5]};
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
%%
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];...
    'lgamma',[45,65];'hgamma',[70,90]};
[psdtrlData,powerPlots] = powerComp(trls,adfreq,bands,'y',[1 1 100],eoi,'n');
%%
[coh,cohPlots] = cohComp(trls,adfreq,eoi,bands,zeroedChannel,[1 1 100],...
    'y','n','mtm','NW',8);
%% Find differences between pre and post
% Collate data
pre = [reshape(psdtrlData{1,1}.relPow,48,38);reshape(coh{1,1}.normBandCoh,168,38)];
post = [reshape(psdtrlData{1,2}.relPow,48,79);reshape(coh{1,2}.normBandCoh,168,79)];
% Get z-scores of post
zPost = (mean(pre,2)-post)./std(pre,[],2);
% Run t-tests
for ii = 1:size(pre,1)
    [~,p(ii)] = ttest2(pre(ii,:),post(ii,:));
end
p = p.*length(p);
%%
baseCoh = coh;
basePow = psdtrlData;
%%
load('ZA1_Stim2_2018_04_11_8chan.mat')
LFPTs.data = filter60(LFPTs.data,adfreq,0);
[LFPTs,adfreq] = dwnSample(LFPTs,4,adfreq);
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,0.5,0.125,1,...
    5,adfreq);
%% Grab baseline and post - set eventTs
start = nearest_idx3(eventTs.t{1,1}(1),LFPTs.tvec);
eventTs.t{1,5} = LFPTs.tvec(1);
eventTs.t{1,6} = LFPTs.tvec(start);
eventTs.label{1,5} = 'Base (Start)';
eventTs.label{1,6} = 'Base (End)';

stop = nearest_idx3(eventTs.t{1}(end),LFPTs.tvec);
eventTs.t{1,7} = LFPTs.tvec(stop);
eventTs.t{1,8} = LFPTs.tvec(end);
eventTs.label{1,7} = 'Post (Start)';
eventTs.label{1,8} = 'Post (End)';
%%
% Separate stim and clean data
stimInd = []; stim = []; stimStop = []; stimStart = [];
% stimInd = logicFind(1,~cellfun(@isempty,eventTs.t),'==');
stim = eventTs.t{1,1};
% Get all stops
stimStop = stim(logicFind(0,round(diff(stim)),'~='));
% Add last stop since above algorithm skips
stimStop = [stimStop;stim(end)];
% Go 0.5 seconds before each stop
stimStop = stimStop-0.5;
% Get starts using the next index after stops from stim
stimStart = stim(logicFind(0,round(diff(stim)),'~=')+1);
% Add first start since above algorithm skips
stimStart = [stim(1);stimStart];
stimStart = stimStart+0.5;
% Use stim start and stops to determine data start and stops; also include
% first/last data point - assumes that no stim starts or ends a recording
dataStart = [LFPTs.tvec(1);LFPTs.tvec(nearest_idx3(stimStop,LFPTs.tvec)+1)'];
dataStop = [LFPTs.tvec(nearest_idx3(stimStart,LFPTs.tvec)-1)';LFPTs.tvec(end)];
struct.t = [];
struct.label = [];
for ii = 1:size(stimStart,1)
    struct.t = [struct.t,{dataStart(ii)},{dataStop(ii)}];
    struct.label = [struct.label,{['Int',num2str(ii),' (Start)']},{['Int',num2str(ii),' (End)']}];
end
struct.t = [struct.t,{LFPTs.tvec(1)},{LFPTs.tvec(start)},{LFPTs.tvec(stop)},{LFPTs.tvec(end)}];
struct.label = [struct.label,{'Base (Start)'},{'Base (End)'},{'Post (Start)'},{'Post (End)'}];
oldEventTs = eventTs;
eventTs = struct;
%%
eoi = {'Base',[0 5]};
for ii = 1:size(stimStart,1)
    eoi = [eoi;{['Int',num2str(ii)],[0 5]}];
end
%%
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
%%
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];...
    'lgamma',[45,65];'hgamma',[70,90]};
[psdtrlData,powerPlots] = powerComp(trls,adfreq,bands,'y',[1 1 100],eoi,'n');
%%
[coh,cohPlots] = cohComp(trls,adfreq,eoi,bands,zeroedChannel,[1 1 100],...
    'y','n','mtm','NW',8);
%% Find differences between pre and post
% Collate data
preBurst = [reshape(psdtrlData{1,1}.relPow,48,62);reshape(coh{1,1}.normBandCoh,168,62)];
for ii = 1:97
    if isempty(psdtrlData{1+ii})
        interBurst(:,:,ii) = NaN;
    else
        interBurst(:,:,ii) = [reshape(psdtrlData{1,1+ii}.relPow(:,:,1),48,1);reshape(coh{1,1+ii}.normBandCoh(:,:,1),168,1)];
    end
end
interBurst = squeeze(interBurst);
% Get z-scores of post
zInterBurst = (mean(preBurst,2)-interBurst)./std(preBurst,[],2);
% Run t-tests
for ii = 1:size(preBurst,1)
    [~,pBurst(ii)] = ttest2(preBurst(ii,:),postBurst(ii,:));
end
pBurst = pBurst.*length(pBurst);
%%
load('ZA1_StimNAc250_2018_04_13.mat')
eventTs.t(:,[1:14,17:end]) = [];
eventTs.label(:,[1:14,17:end]) = [];

LFPTs.data = filter60(LFPTs.data,adfreq,0);
[LFPTs,adfreq] = dwnSample(LFPTs,4,adfreq);
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,1,0.125,0.5,...
    5,adfreq);
%% Grab baseline and post - set eventTs
start = nearest_idx3(eventTs.t{1,1}(1),LFPTs.tvec);
eventTs.t{1,3} = LFPTs.tvec(1);
eventTs.t{1,4} = LFPTs.tvec(start);
eventTs.label{1,3} = 'Base (Start)';
eventTs.label{1,4} = 'Base (End)';

stop = nearest_idx3(eventTs.t{2}(end),LFPTs.tvec);
eventTs.t{1,5} = LFPTs.tvec(stop);
eventTs.t{1,6} = LFPTs.tvec(end);
eventTs.label{1,5} = 'Post (Start)';
eventTs.label{1,6} = 'Post (End)';
%%
eoi = {'base',[0 5];'post',[0 5]};
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
bands = {'delta',[1,4];'theta',[5,10];'alpha',[11,14];'beta',[15,30];...
    'lgamma',[45,65];'hgamma',[70,90]};
[psdtrlData,powerPlots] = powerComp(trls,adfreq,bands,'y',[1 1 100],eoi,'n');
[coh,cohPlots] = cohComp(trls,adfreq,eoi,bands,zeroedChannel,[1 1 100],...
    'y','n','mtm','NW',8);
%%
pre = [reshape(psdtrlData{1,1}.relPow,48,size(psdtrlData{1}.Pow,3));reshape(coh{1,1}.normBandCoh,168,size(psdtrlData{1}.Pow,3))];
post = [reshape(psdtrlData{1,2}.relPow,48,size(psdtrlData{2}.Pow,3));reshape(coh{1,2}.normBandCoh,168,size(psdtrlData{2}.Pow,3))];
% Get z-scores of post
zPost = (post-mean(pre,2))./std(pre,[],2);
figure
imagesc(zPost)
colormap viridis
%%
for ii = 1:8
    [preS(:,:,ii),f,t] = spectrogram(LFPTs.data(ii,1:nearest_idx3(eventTs.t{1}(1),LFPTs.tvec)),1024,512,1:100,adfreq);
    [postS(:,:,ii),f,t] = spectrogram(LFPTs.data(ii,nearest_idx3(eventTs.t{2}(end),LFPTs.tvec):end),1024,512,1:100,adfreq);
end
%%
preP = 10*log10(abs(preS).^2);
mPreP = mean(preP,2,'omitnan');
sPreP = nanstd(preP,[],2);
postP = 10*log10(abs(postS).^2);
zPostP = (postP-repmat(mPreP,1,size(postP,2),1))./repmat(sPreP,1,size(postP,2),1);
%%
figure
for ii = 1:8
    subplot(2,4,ii)
    h = pcolor(zPostP(:,1:end,ii));
    set(h,'EdgeColor','none')
    colormap viridis
    caxis([-3 3])
end
%%
z = (psdtrlData{2}.Pow-repmat(psdtrlData{1}.Overall,1,1,89))./repmat(psdtrlData{1}.OverallStd,1,1,89);
figure
for ii = 1:8
    subplot(2,4,ii)
    h = pcolor(squeeze(z(ii,:,1:30)));
    set(h,'EdgeColor','none')
    colormap viridis
    caxis([-3 3])
end
%%
load('ZA1_StimNAc22250_2018_04_16.mat')
eventTs.t(:,[1:14,17:end]) = [];
eventTs.label(:,[1:14,17:end]) = [];
