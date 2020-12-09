%% Determine latencies
files = fileSearch('F:/irdmRound2/processedContinuous3/','.mat');
allTrials = [];
thesePercent = [];
for ii = 1:numel(files)
    load(files{ii},'trials')
    if ~isempty(trials)
        allTrials = [allTrials;trials];
        % Use the block percent of each first trial
        this = struct2cell(trials);
        theseInds = cell2mat(this(5,:))==1;
        thesePercent = [thesePercent;cell2mat(this(12,theseInds))'];
    end
end
allCell = struct2cell(allTrials);
% Remove doubleTones
allCell(:,cellfun(@(x) strcmp(x,'doubleTone'),allCell(1,:))) = [];
% Get tone to nose poke latency
tnpLat = cell2mat(allCell(9,:));
% Get nose poke to lever press latency
nplpLat = cell2mat(allCell(10,:));
% Combine to get tone to lever press latency
tlpLat = tnpLat + nplpLat;
% Convert block percent and blockN
percent = cell2mat(allCell(12,:));
blockN = cell2mat(allCell(4,:));
% Use trialN = 1 to grab percent from each block
trialN = cell2mat(allCell(5,:));
blockPercent = percent(trialN==1);
% Calculate ecdfs
[tlpY,tlpX] = ecdf(tlpLat);
[nplpY,nplpX] = ecdf(nplpLat);
%% Plot histogram of block percents
figure
subplot(2,3,1)
histogram(blockPercent,'normalization','probability')
box off
xlabel('block % delay choice')
ylabel('% of blocks')
title('all blocks')
% Plot block percents through blocks
for ii = 1:5
    subplot(2,3,ii+1)
    histogram(percent(blockN==ii),'normalization','probability')
    box off
    xlabel('block % delay choice')
    ylabel('% of blocks')
    title(['block ',num2str(ii)])
end
%% Plot % delay by delay time
delay = cell2mat(allCell(13,:));
theseDelay = delay(trialN==1);
figure
plot(theseDelay,blockPercent,'.','markersize',10)
%% Plot latency histograms and ecdfs
figure
subplot(1,3,1)
histogram(nplpLat)
title('nose poke to lever press')
xlim([0 20])
xlabel('latency (sec)')
ylabel('count')

subplot(1,3,2)
histogram(tlpLat)
xlabel('latency (sec)')
title('tone to lever press')
xlim([0 20])

subplot(1,3,3)
plot(nplpX,nplpY)
hold on
plot(tlpX,tlpY)
xlim([0 20])
xlabel('latency (sec)')
ylabel('cumulative probability')
legend({'nose poke','tone'})
%% Immediate vs delay
imInds = cellfun(@(x) strcmp(x,'immediate_free'),allCell(1,:));
[imNoseY,imNoseX] = ecdf(nplpLat(imInds));
[imToneY,imToneX] = ecdf(tlpLat(imInds));
deInds = cellfun(@(x) strcmp(x,'delay_free'),allCell(1,:));
[deNoseY,deNoseX] = ecdf(nplpLat(deInds));
[deToneY,deToneX] = ecdf(tlpLat(deInds));
figure
subplot(1,2,1)
title('tone to lever press')
hold on
plot(imToneX,imToneY)
plot(deToneX,deToneY)
xlabel('latency (sec)')
ylabel('cumulative probability')
xlim([0 20])
legend({'immediate','delay'})

subplot(1,2,2)
title('nose poke to lever press')
hold on
plot(imNoseX,imNoseY)
plot(deNoseX,deNoseY)
xlabel('latency (sec)')
ylabel('cumulative probability')
xlim([0 20])
legend({'immediate','delay'})
%% Break into groups based on block number
figure
title('tone to lever press by block')
% title('nose poke to lever press by block')
hold on
[blockX,blockY] = deal(cell(1,5));
ninetyfve = zeros(1,5);
for ii = 1:5
    [blockY{ii},blockX{ii}] = ecdf(tlpLat(blockN==ii));
%     [blockY{ii},blockX{ii}] = ecdf(nplpLat(blockN==ii));
    % Get 95 percentile
    ninetyfive(ii) = blockX{ii}(nearest_idx3(0.95,blockY{ii}));
    mBlock(ii) = mean(tlpLat(blockN==ii),'omitnan');
%     mBlock(ii) = mean(nplpLat(blockN==ii),'omitnan');
    plot(blockX{ii},blockY{ii})
end
plot([0 20],[0.95 0.95],':k')
xlabel('latency (sec)')
ylabel('cumulative probability')
xlim([0 20])
legend({'1','2','3','4','5'})
axes('position',[0.7 0.2 0.2 0.2])
plot(1:5,ninetyfive,'.k')
hold on
plot(1:5,mBlock,'.r')
box off
xlabel('block')
ylabel('latency (sec)')
legend({'95%','mean'},'location','nw')
set(gca,'xtick',1:5,'xlim',[0 6],'ylim',[0 7])
%% Break into groups based on block percent delay
% Round block percent delays
rBlockDelay = round(cell2mat(allCell(12,:)),1);
for ii = 1:11
        [delayY{ii},delayX{ii}] = ecdf(nplpLat(rBlockDelay==(ii-1)/10));
%     [delayY{ii},delayX{ii}] = ecdf(tlpLat(rBlockDelay==(ii-1)/10));
    ninetyfiveDelay(ii) = delayX{ii}(nearest_idx3(0.95,delayY{ii}));
        m(ii) = mean(nplpLat(rBlockDelay==(ii-1)/10),'omitnan');
%     m(ii) = mean(tlpLat(rBlockDelay==(ii-1)/10),'omitnan');
end
figure
% title('tone to lever press by % delay choice')
title('nose poke to lever press by % delay choice')
hold on
plot(0:0.1:1,ninetyfiveDelay,'.k','markersize',20)
plot(0:0.1:1,m,'.r','markersize',20)
xlim([-0.1 1.1])
legend({'95%','mean'})
set(gca,'xticklabel',0:0.2:1)
xlabel('% delay choice')
ylabel('latency (sec)')
%% Get clean data associated with free decisions
files = fileSearch('F:/irdmRound2/processedContinuous3/','DDT');
% cd('G:\GreenLab\data\irdmNew\processedContinuous\')
data = []; c = 0; aC = 1;
% Set window size (sec) desired and % of window required to keep
winSize = 10;
for ii = 1:numel(files)
    disp(ii)
    load(files{ii},'trials')
    parts = strsplit(files{ii},'_');
    if exist('trials','var') && ~isempty(trials)
        % Only load whole file if trials are good
        load(files{ii},'trls','psdTrls','coh','hist')
        % check for 8 channels
        if size(trls{1}.trial,1) == 8
            % Convert struct to cell
            this = struct2cell(trials);
            % Setup first; used to determine which is the first trial; jj=1
            % may not be used
            first = 1;
            if trials(end).blockN == 5 && trials(end).trialN > 5
                for jj = 1:size(this,1)
                    % Subtract 0.5 seconds from trial start to account for
                    % feature times given as center of 1-second window
                    % (i.e., window centered at 0.5 includes data from [0
                    % 1])
                    if strcmp('immediate_free',this(1,jj)) || strcmp(...
                            'delay_free',this(1,jj))
                        stop = nearest_idx3(cell2mat(this(2,jj))-0.5,...
                            psdTrls{1}.t,-1);
                        start = stop-winSize+1;
                        % Check for NaNs; only keep if 3/5 not Nans (using
                        % first values in relPower as proxy)
                        if sum(~isnan(psdTrls{1}.relPow(1,1,1,start:...
                                stop)))>=1
                            % Add to counter
                            c = c+1;
                            % Get starting index for this session
                            if first == 1
                                animalInds(aC,1) = c;
                            end
                            data(c).animal = parts{1};
                            data(c).date = parts{3};
                            
                            
                            data(c).rawPower = psdTrls{1}.Pow(:,:,:,...
                                start:stop);
                            data(c).relPower = psdTrls{1}.relPow(:,:,:,...
                                start:stop);
                            
                            data(c).rawCoh = coh{1}.Cxy(:,:,:,start:stop);
                            data(c).normBandCoh = coh{1}.normBandCoh(:,...
                                :,:,start:stop);
                            
                            if strcmp('immediate_free',this(1,jj))
                                data(c).delay = 0;
                            else
                                data(c).delay = 1;
                            end
                            
                            data(c).blockN = cell2mat(this(4,jj));
                            data(c).trialN = cell2mat(this(5,jj));
                            data(c).leverPress = cell2mat(this(6,jj));
                            data(c).feeder = cell2mat(this(7,jj));
                            data(c).headEntry = cell2mat(this(8,jj));
                            data(c).tone_nosePoke_latency = ...
                                cell2mat(this(9,jj));
                            data(c).nosePoke_lever_latency = ...
                                cell2mat(this(10,jj));
                            data(c).lever_headEntry_latency = ...
                                cell2mat(this(11,jj));
                            data(c).blockD = cell2mat(this(12,jj));
                            % Reset first
                            first = 0;
                        end
                    end
                end
                % Get ending index for this session
                animalInds(aC,2) = c;
                aC = aC+1;
            end
        end
    end
    clear trials
end
% Get rid of zero animalInds due to NaNed data
animalInds(animalInds(:,1)==0,:) = [];
%%
for ii = 1:size(data,2)
    clean(ii,:) = ~isnan(data(ii).rawPower(1,1,1,:));
end
figure
subplot(1,3,1)
imagesc(clean)
set(gca,'xtick',1:2:10,'xticklabel',1.5:2:10.5)
title('clean data per lever press')
colormap gray
xlabel('time before lever press (sec)')
ylabel('lever press')
% get sum of clean trials through time
sum(clean);
for ii = 1:10
    cumSum(ii) = sum(reshape(clean(:,end-ii+1:end),1,ii*size(clean,1)));
end
subplot(1,3,2)
plot(1:10,cumSum,'.k','markersize',20)
title('total clean windows through time')
set(gca,'xtick',1:2:10,'xticklabel',1.5:2:10.5)
xlabel('time before lever press (sec)')
ylabel('total data')
% percent of window available
for ii = 1:size(clean,1)
   for jj = 1:10
        percentClean(ii,jj) = sum(clean(ii,1:jj))/jj;
   end
end
% Get percent clean by animal
% Convert and split apart anmimal ID to get just IRDM# (removes things
% after -, e.g., IRDM11-stimcore would count as IRDM11)
cellData = squeeze(struct2cell(data));
allID = unique(cellData(1,:));
splitID = cellfun(@(x) strsplit(x,'-'),allID,'UniformOutput',false);
uIDs = unique(cellfun(@(x) x{1},splitID,'UniformOutput',0));
for ii = 1:numel(uIDs)
    % Find all indices of data that are from given animal
    inds = logicFind(1,cellfun(@(x) ~isempty(x),strfind(cellData(1,:),...
        uIDs{ii})),'==');
    c = 1;
    for jj = inds
        for k = 1:10
            ratPercentClean{ii}(c,k) = sum(clean(jj,1:k))/k;
        end 
        c = c+1;
    end
end
subplot(1,3,3)
hold on
for ii = 1:numel(ratPercentClean)
    plot(1:10,mean(ratPercentClean{ii}),'-','color',[0.5 0.5 0.5])%,'.','markerSize',20)
end
plot(1:10,mean(percentClean,1),'.k','markersize',20)
set(gca,'xtick',1:2:10,'xticklabel',1.5:2:10.5)
xlabel('time before lever press (sec)')
ylabel('% of clean data')
title('% of clean data through time')
%% Get trial number by length of window and # clean required
this = [];
for ii = 1:size(clean,1)
    for jj = 1:10
        for k = 1:jj
           this(ii,k,jj) = sum(clean(ii,1:jj)) >= k; 
        end
    end
end
figure
imagesc(squeeze(sum(this,1)))
set(gca,'xtick',1:2:10,'xticklabel',1.5:2:10.5)
xlabel('time before lever press (sec)')
ylabel('# clean bins required')
c = colorbar;
c.Label.String = 'trials';
colormap viridis
%% Immediate and delay lever press
allData = squeeze(struct2cell(data));
delay = []; immediate = [];
% Find minimal number of delays and immediates each session contributes
for ii = 1:size(animalInds,1)
    % Find indices of delays and immediates
    delayInd = logicFind(1,cell2mat(allData(7,animalInds(ii,1):...
        animalInds(ii,2))),'==');
    immediateInd = logicFind(0,cell2mat(allData(7,animalInds(ii,1):...
        animalInds(ii,2))),'==');
    % Count number of delays and immediates
    dN(ii) = numel(delayInd);
    iN(ii) = numel(immediateInd);
    % Find min
    thisMin = min([dN(ii),iN(ii)]);
    if ii == 1
        start = 0;
    else
        start = animalInds(ii-1,2);
    end
    % Pull that number of both trials
    delay = cat(2,delay,allData(:,start+delayInd(randperm(dN(ii),...
        thisMin))));
    immediate = cat(2,immediate,allData(:,start+immediateInd(randperm(...
        iN(ii),thisMin))));
end
% Count number of trials per animal (just counts delays, so half of real
% number)
u = unique(delay(1,:));
for ii = 1:numel(u)
    nD(ii) = sum(strcmp(u(ii),delay(1,:)));
end
%% Decisive split - decisive = >80% one lever or the other, indecisive = 
% <40% one lever or the other
for ii = 1:size(animalInds,1)
    delayPerc(ii) = data(animalInds(ii,1)).blockD;
end
allDelay = extractfield(data,'blockD');
topInd = logicFind(1,allDelay<=0.2 | allDelay>=0.8,'==');
midInd = allDelay<=0.6 & allDelay>=0.4;

top = squeeze(struct2cell(data(topInd(randperm(numel(topInd),sum(midInd))))));
mid = squeeze(struct2cell(data(midInd)));
%%
figure
set(gcf,'renderer','Painters')
hold on
histogram(delayPerc,'binwidth',.1)
plot([0.2 0.2],[0 16],'-k')
plot([0.8 0.8],[0 16],'-k')
plot([0.6 0.6],[0 16],'--k')
plot([0.4 0.4],[0 16],'--k')
box off

%% Median split AUC
auc{1} = [345,325,365,355,355,360,370,375,365,350,365,365,365];
auc{2} = [215,210,225,220,180,100,180,265,225,220];
auc{3} = [95,90,60,40,75,70,55,55,70,75,85,120];
auc{4} = [295,260,280,235,180,135,170,185,160,110,210,210,180];
auc{5} = [120,135,225,200,160,155,190,185,165,145,190,195,180];
auc{6} = [95,20,55,240,100,105,90,75,80,90,25,65];
auc{7} = [80,90,120,325,175,145,120];
allAuc= cat(2,auc{:})';
figure
set(gcf,'renderer','Painters')
histogram(allAuc,20)
xlabel('AUC')
ylabel('# of sessions')
hold on
plot([median(allAuc) median(allAuc)],[0 10],'--k')
box off
%% Load IRDM model
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\delayImmediateGen\')
naive = []; naiveP = []; cv = []; cvP = [];
for ii = 1:100
    load(['delayImmediateGen',num2str(ii),'.mat'])
    naive(ii) = accArray{1}.acc;
    naiveP(ii) = accArrayP{1}.acc;
    
    cv = cat(1,cv,allLambda{1}.allErr);
    cvP = cat(1,cvP,allLambdaP{1}.allErr);
    
    [x(ii,:),y(ii,:),~,a(ii)] = perfcurve(hist.cfg.naive.testY,accArray{1}.pred,1,'tvals',0:0.01:1,'usenearest',0);
    [xP(ii,:),yP(ii,:),~,aP(ii)] = perfcurve(histP.cfg.naive.testY,accArrayP{1}.pred,1,'tvals',0:0.01:1,'usenearest',0);
end
figure
hold on
plot(mean(x,1),mean(y,1),'-k')
plot(mean(xP,1),mean(yP,1),'--k')

mA = round(mean(a),2);
cA = round(conf(a,0.95),2);

mAp = round(mean(aP),2);
cAp = round(conf(aP,0.95),2);
legend({['Real: \mu = ',num2str(mA),'\pm',num2str(cA)],['Permuted: \mu = ',num2str(mAp),'\pm',num2str(cAp)]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR');
%%
dC = 1; iC = 1;
for ii = 1:1050
    if data(ii).delay == 1
        delayRaw(dC,:,:) = squeeze(mean(data(ii).rawPower,4,'omitnan'));
        delayRel(dC,:,:) = squeeze(mean(data(ii).relPower,4,'omitnan'));
        dC = dC+1;
    else
        immediateRaw(iC,:,:) = squeeze(mean(data(ii).rawPower,4,'omitnan'));
        immediateRel(iC,:,:) = squeeze(mean(data(ii).relPower,4,'omitnan'));
        iC = iC+1;
    end
end
%% Within animal models - find bad channels per animal
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\irdm\'...
    'toProcess'],'.mat');
parts = cellfun(@(x) strsplit(x,'_'),files,'UniformOutput',0);
parts = cat(1,parts{:});
u = unique(parts(:,1));
for ii = 1:numel(u)
    theseInds = logicFind(1,strcmp(u(ii),parts(:,1)),'==');
    c = 1;
    for jj = theseInds
        load(files{jj},'LFPTs')
        over{ii}(:,c) = sum(abs(LFPTs.data)>2,2)./length(LFPTs.data);
        c = c+1;
    end
    over50{ii} = over{ii}>=0.3;
    med(:,ii) = median(over{ii},2);
    avg(:,ii) = mean(over{ii},2);
end
%% Flexibility within block (decisiveness)
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\topMid2\')
naive = []; naiveP = []; cv = []; cvP = [];
for ii = 1:100
    load(['topMid',num2str(ii),'.mat'])
    naive(ii) = accArray{1}.acc;
    naiveP(ii) = accArrayP{1}.acc;
    
    cv = cat(1,cv,allLambda{1}.allErr);
    cvP = cat(1,cvP,allLambdaP{1}.allErr);
    
    [x(ii,:),y(ii,:),~,a(ii)] = perfcurve(hist.cfg.naive.testY,accArray{1}.pred,1,'tvals',0:0.01:1,'usenearest',0);
    [xP(ii,:),yP(ii,:),~,aP(ii)] = perfcurve(histP.cfg.naive.testY,accArrayP{1}.pred,1,'tvals',0:0.01:1,'usenearest',0);
end
figure
hold on
plot(mean(x,1),mean(y,1),'-k')
plot(mean(xP,1),mean(yP,1),'--k')

mA = round(mean(a),2);
cA = round(conf(a,0.95),2);

mAp = round(mean(aP),2);
cAp = round(conf(aP,0.95),2);
legend({['Real: \mu = ',num2str(mA),'\pm',num2str(cA)],['Permuted: \mu = ',num2str(mAp),'\pm',num2str(cAp)]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR');
%% Run individual animal models
allData = squeeze(struct2cell(data));
delay = []; immediate = [];
% Find minimal number of delays and immediates each session contributes
for jj = 1:size(animalInds,1)
    % Find indices of delays and immediates
    delayInd = logicFind(1,cell2mat(allData(7,animalInds(jj,1):...
        animalInds(jj,2))),'==');
    immediateInd = logicFind(0,cell2mat(allData(7,animalInds(jj,1):...
        animalInds(jj,2))),'==');
    % Count number of delays and immediates
    dN(jj) = numel(delayInd);
    iN(jj) = numel(immediateInd);
    % Find min
    thisMin = min([dN(jj),iN(jj)]);
    if jj == 1
        start = 0;
    else
        start = animalInds(jj-1,2);
    end
    % Pull that number of both trials
    delay = cat(2,delay,allData(:,start+delayInd(randperm(dN(jj),...
        thisMin))));
    immediate = cat(2,immediate,allData(:,start+immediateInd(randperm(...
        iN(jj),thisMin))));
end
u = unique(delay(1,:));
for ii = 1:numel(u)
    theseInds = strcmp(u(ii),delay(1,:));
    % Delay
    delayRelPower = cellfun(@(x) reshape(squeeze(mean(x,4,'omitnan')),1,48),...
        delay(4,theseInds),'UniformOutput',0);
    catDelayPower = cat(1,delayRelPower{:});
    delayNormCoh = cellfun(@(x) reshape(squeeze(mean(x,4,'omitnan')),1,168),...
        delay(6,theseInds),'UniformOutput',0);
    catDelayCoh = cat(1,delayNormCoh{:});
    % Immediate
    immediateRelPower = cellfun(@(x) reshape(squeeze(mean(x,4,'omitnan')),...
        1,48),immediate(4,theseInds),'UniformOutput',0);
    catImmediatePower = cat(1,immediateRelPower{:});
    immediateNormCoh = cellfun(@(x) reshape(squeeze(mean(x,4,'omitnan')),...
        1,168),immediate(6,theseInds),'UniformOutput',0);
    catImmediateCoh = cat(1,immediateNormCoh{:});
    % Combine
    indDataX{ii} = [catDelayPower,catDelayCoh;catImmediatePower,catImmediateCoh];
    indDataY{ii} = [ones(numel(delayRelPower),1);zeros(numel(immediateRelPower),1)];
    % Split into train and test
    x = indDataX{ii};
    y = indDataY{ii};
    trainInds = randperm(numel(y),ceil(numel(y)*0.8));
    trainX = x(trainInds,:);
    trainY = y(trainInds,:);
    testInds = 1:numel(y);
    testInds = ~ismember(testInds,trainInds);
    testX = x(testInds,:);
    % Build models
    cfg = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se',[]);
    [~,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x,y,...
        'binomial','deviance',1,10,1,cfg);
    
    cfg = lassoNetCfg({testX,testY},[],'y','y','n',100,'1se',[]);
    [~,allLambdaP,allBetaP,cvFitsArrayP,accArrayP,histP] = lassoNet(x,y,...
        'binomial','deviance',1,10,1,cfg);
    
    save(['/ihome/ldwiel/data/delayImmediate',u{ii},'_',num2str(ii+20*(k-1)),'.mat'],...
        'allLambda','allBeta','cvFitsArray','accArray','hist','allLambdaP',...
        'allBetaP','cvFitsArrayP','accArrayP','histP')
end
%% Ind models
fSearch = {'IRDM11','IRDM14','IRDM15','IRDM16','IRDM18','IRDM21','IRDM22'};
for ii = 1:numel(fSearch)
    files{ii} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\irdm\delayImmediateInd',fSearch{ii});
    for jj = 1:numel(files{ii})
        load(files{ii}{jj})
        auc{ii}(jj) = accArray{1}.acc;
        aucP{ii}(jj) = accArrayP{1}.acc;
        
        [x{ii}(jj,:),y{ii}(jj,:),~,a{ii}(jj)] = perfcurve(hist.cfg.naive.testY,accArray{1}.pred,1,'tvals',0:0.01:1,'usenearest',0);
        [xP{ii}(jj,:),yP{ii}(jj,:),~,aP{ii}(jj)] = perfcurve(histP.cfg.naive.testY,accArrayP{1}.pred,1,'tvals',0:0.01:1,'usenearest',0);
    end
end
for ii = 1:numel(fSearch)
    aM(ii) = mean(a{ii});
    apM(ii) = mean(aP{ii});
    aC(ii) = conf(a{ii},0.95);
    apC(ii) = conf(a{ii},0.95);
%     aucM(ii) = mean(auc{ii});
%     aucpM(ii) = mean(aucP{ii});
end
figure
hold on
for ii = 1:numel(fSearch)
    plot(mean(x{ii},1),mean(y{ii},1))
end
plot(mean(xP{1}),mean(yP{1}),'--k')
legend({'IRDM11: 76','IRDM14: 162','IRDM15: 114','IRDM16: 109','IRDM18: 280','IRDM21: 59','IRDM22: 22','Permuted: '})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR');
%% Ind p-values
for ii = 1:numel(fSearch)
    z(ii) = (aM(ii)-mean(aP{ii}))/std(aP{ii});
end
p = 2*normcdf(-abs(z));
%% Tone - nose poke latency
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\toneNose\')
for ii = 1:100
    load(['toneNose',num2str(ii),'.mat'])
    acc(ii) = mean(abs(accArray{1}.acc-hist.cfg.naive.testY'));
    accP(ii) = mean(abs(accArrayP{1}.acc-histP.cfg.naive.testY'));
end
doubleHist(acc,accP,'main','Tone -> Nose Poke')
% load('C:\Users\Pythia\Documents\GreenLab\data\irdm\irdmDelayImmediateDataMinus1.mat')
% allData = squeeze(struct2cell(data));
% delay = cell2mat(allData(13,cell2mat(allData(7,:))==1));
% immediate = cell2mat(allData(13,cell2mat(allData(7,:))==0));
% figure
% violin({delay',immediate'})
%% Nose poke - lever latency
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\noseLever\')
for ii = 1:100
    load(['noseLever',num2str(ii),'.mat'])
    acc(ii) = mean(abs(accArray{1}.acc-hist.cfg.naive.testY'));
    accP(ii) = mean(abs(accArrayP{1}.acc-histP.cfg.naive.testY'));
end
doubleHist(acc,accP,'main','Nose Poke -> Lever')
%%
load('C:\Users\Pythia\Documents\GreenLab\data\irdm\irdmDelayImmediateDataMinus1.mat')
allData = squeeze(struct2cell(data));
delay = cell2mat(allData(14,cell2mat(allData(7,:))==1));
immediate = cell2mat(allData(14,cell2mat(allData(7,:))==0));
%%
figure
% set(gcf,'renderer','Painters')
hold on
plotSpread({delay',immediate'})
plot([1 1],[mean(delay) mean(delay)],'ok')
plot([2 2],[mean(immediate) mean(immediate)],'ok')
% plot([1 1],[mean(delay)-conf(delay,0.95) mean(delay)+conf(delay,0.95)],'-r')
% plot([2 2],[mean(immediate)-conf(immediate,0.95) mean(immediate)+conf(immediate,0.95)],'-r')
set(gca,'xticklabel',{'Delay','Immediate'})
ylabel('Latency (sec)')
%% Median split AUC plot
figure
hold on
scatter(1:8,[343,241,63,200,146,96,177,225],200,'.')
plot([0.5 8.5],[188.5 188.5],'--')
xlim([0.5 8.5])
%% Raw AUC plot
auc{1} = [345,325,365,355,355,360,370,375,365,350,365,365,365];
auc{2} = [215,210,225,220,180,100,180,265,225,220];
auc{3} = [95,90,60,40,75,70,55,55,70,75,85,120];
auc{4} = [295,260,280,235,180,135,170,185,160,110,210,210,180];
auc{5} = [120,135,225,200,160,155,190,185,165,145,190,195,180];
auc{6} = [95,20,55,240,100,105,90,75,80,90,25,65];
auc{7} = [80,90,120,325,175,145,120];
for ii = 1:7  
    thisAucZ = zscore(auc{ii});
    aucZ{ii} = thisAucZ';
end
allAucZ = cat(1,aucZ{:});
allAuc= cat(2,auc{:})';
% Raw AUC plot
figure
pcolor(padarray(flipud(repmat(allAuc,1,2)),[1 0],0))
colormap viridis
% ylim([2 80])
figure
pcolor(padarray(flipud(repmat(allAuc(randperm(80,80)),1,2)),[1 0],0))
colormap viridis
% ylim([2 80])
%% Zscore AUC plot
figure
pcolor(padarray(flipud(repmat(allAucZ,1,2)),[1 0],0))
colormap viridis
% ylim([2 80])

figure
trainInds = [1:10,14:21,24:33,36:45,49:58,62:71,74:79];
pcolor(padarray(flipud(repmat(allAucZ(trainInds),1,2)),[1 0],0))
colormap viridis
caxis([-2.6 2.6])
% ylim([2 65])

figure
testInds = [11:13,22:23,34:35,45:48,59:61,72:73,80];
pcolor(padarray(flipud(repmat(allAucZ(testInds),1,2)),[1 0],0))
colormap viridis
caxis([-2.6 2.6])
% ylim([2 18])

permZ = allAucZ(randperm(80,80));
figure
pcolor(padarray(flipud(repmat(permZ(trainInds),1,2)),[1 0],0))
colormap viridis
caxis([-2.6 2.6])
% ylim([2 65])

figure
pcolor(padarray(flipud(repmat(permZ(testInds),1,2)),[1 0],0))
colormap viridis
caxis([-2.6 2.6])
% ylim([2 18])
%% talk figure
% load data
load('IRDM18_DDT_2019-03-15.mat')
% run ddtTrial
trial = ddtTrials(eventTs,LFPTs,1);
hold on
LFPTs.data = filter60(LFPTs,2000,0);
LFPTs = dwnSample(LFPTs,5,2000);
for ii = 1:8
plot(LFPTs.tvec,LFPTs.data(ii,:)+10+ii)
end
ylim([0 20])
xlim([1000 1090])
for ii = 1:62
    plot([trial(ii).start trial(ii).start],[0 70],'-k')
    plot([trial(ii).start-5 trial(ii).start-5],[0 70],'--k')
end
%% Data for coherence and power figures
load('C:\Users\Pythia\Documents\GreenLab\data\irdm\matNew\IRDM18_DDT_2019-03-15.mat')
LFPTs.data = filter60(LFPTs,2000,0);
LFPTs = dwnSample(LFPTs,5,2000);
%% Coherence
sigInd = 350000:352000;
% Signals
figure
plot(t,smooth(LFPTs.data(2,sigInd),20))
hold on
plot(t,smooth(LFPTs.data(3,sigInd),20))
box off
% 3 Hz filter
[b,a] = cheby2(6,40,0.015); % 0.015 = 3/200
test = filtfilt(b,a,LFPTs.data(2,sigInd));
test2 = filtfilt(b,a,LFPTs.data(3,sigInd));
t = 0:1/400:5;
h1 = hilbert(test);
h2 = hilbert(test2);
p1 = angle(h1);
p2 = angle(h2);
figure
hold on
plot(t,test)
plot(t,test2)
figure
plot(p1,p2,'.')
set(gca,'xtick',-4:2:4,'ytick',-4:2:4)
lsline
% 85 Hz filter
d = designfilt('bandpassfir','FilterOrder',20, ...
    'CutoffFrequency1',85,'CutoffFrequency2',86, ...
    'SampleRate',400);
test = filtfilt(d,LFPTs.data(2,sigInd));
test2 = filtfilt(d,LFPTs.data(3,sigInd));
t = 0:1/400:5;

h1 = hilbert(test);
h2 = hilbert(test2);
p1 = angle(h1);
p2 = angle(h2);

figure
hold on
plot(t,smooth(test))
plot(t,smooth(test2))
xlim([1 1.2])

figure
plot(p1,p2,'.')

set(gca,'xtick',-4:2:4,'ytick',-4:2:4)
lsline
corrcoef(p1,p2)
% Sine wave 
figure
t = 1.5:.01:2.5;
x = sin(2*pi*t);
plot(t,x)
set(gca,'xtick',1:1/4:5)
box off
%% Power
% sigInd = [nearest_idx3(1042,LFPTs.tvec):nearest_idx3(1047,LFPTs.tvec)];
sigInd = 350000:352000;
[Pxx,F] = pwelch(LFPTs.data(2,sigInd),256,128,1:1:100,400);

fig = figure;
set(gcf,'color',[1 1 1])
im = cell(1,1);
axis tight manual

sigBack = plot3(ones(2001,1)*20,LFPTs.tvec(sigInd),...
    flipud(smooth(LFPTs.data(2,sigInd),40)).*350,'-b');
axis([0 20 875 880 -60 60]);
xlabel('Freq (Hz)')
ylabel('Time (s)')
zlabel('Magnitude (au)')
view(90,0)
frame = getframe(fig);
im{1} = frame2im(frame);
pause(1)

hold on

inds = [6,10;11,16;21,25;41,45;81,85];
freq = [1,2,4,8,16];
t = 0:1/400:5;
y = fft(LFPTs.data(2,sigInd));
l = 2000;
f = 400*(0:(l/2))/l;

thisXtick = [];
for ii = 1:size(inds,1)
    phase(ii) = mean(angle(y(inds(ii,1):inds(ii,2))))*pi;
    mag(ii) = 10*log10(abs(Pxx(freq(ii))))+71;
    x(ii,:) = sin(2*pi*freq(ii)*t+phase(ii))*mag(ii);
    if ii > 1
        for jj = 1:ii-1
            this{jj}.Color = [0.8 0.8 0.8];
        end
    end
    this{ii} = plot3(ones(1,numel(t)).*freq(ii),LFPTs.tvec(350000:352000),x(ii,:),'color',[0.4660 0.6740 0.1880]);
    thisXtick = cat(1,thisXtick,freq(ii));
    set(gca,'xtick',thisXtick)
    thisText = text(0,877.5,40,[num2str(freq(ii)),' Hz'],'Color',[0.4660 0.6740 0.1880],'FontSize',15,'HorizontalAlignment','center');
    frame = getframe(fig);
    im = [im,{frame2im(frame)}];
    pause(1)
    set(thisText,'visible','off')
end
view(-90,0)
set(sigBack,'visible','off')
for ii = 1:size(inds,1)
    set(this{ii},'visible','off')
    this{ii} = plot3(ones(1,numel(t)).*freq(ii),LFPTs.tvec(350000:352000),fliplr(x(ii,:)),'color',[0.8 0.8 0.8]);
    set(this{ii},'visible','on')
end

sigFor = plot3(zeros(2001,1),LFPTs.tvec(sigInd),...
    smooth(LFPTs.data(2,sigInd),40).*350,'-b');
this{end}.Color = [0.8 0.8 0.8];

% Turn
for jj = 1:90
    view(-90+jj,0)
    frame = getframe(fig);
    im = [im,{frame2im(frame)}];
    pause(0.02)
end

pause(1)
set(sigFor,'visible','off')
frame = getframe(fig);
im = [im,{frame2im(frame)}];
pause(0.5)
zeroLine = plot3([0 100],[875 875],[-0.1 -0.1],'-k','linewidth',1);
frame = getframe(fig);
im = [im,{frame2im(frame)}];
pause(0.5)

for ii = 1:size(inds,1)
    set(this{ii},'visible','off')
    stem(ii) = plot3([freq(ii) freq(ii)],[875 875],[0 abs(mag(ii))],'-r');
    plot3([freq(ii) freq(ii)],[875 875],abs([mag(ii) mag(ii)]),'.r')
end
frame = getframe(fig);
im = [im,{frame2im(frame)}];
pause(1)
set(stem,'visible','off')
frame = getframe(fig);
im = [im,{frame2im(frame)}];
pause(1)

dot = plot3(F,ones(1,100)*875,10*log10(abs(Pxx))+71,'.r');
axis([0 20 875 880 -60 60]);
set(gca,'xtick',0:10:100)
for ii = 1:60%lcm(60,80)
    axis([0 20+ii/.75 875 880 -60+ii 60])
    frame = getframe(fig);
    im = [im,{frame2im(frame)}];
end
set(dot,'visible','off')
plot3(F,ones(1,100)*875,10*log10(abs(Pxx))+71,'-r')
frame = getframe(fig);
    im = [im,{frame2im(frame)}];
    % Add subplot with orignal signal and psd side by side
    % Remove numbers on magnitude axis
%%
filename = 'powerDot.gif';
for idx = 1:numel(im)
   [A,map] = rgb2ind(im{idx},256);
   if idx == 1
       imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',2);
   else
       if ismember(idx,2:6)
           imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1.5);
       end
       if ismember(idx,7:95)
           imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.03);
       end
       if ismember(idx,96:100)
           imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1.5);
       end  
       if ismember(idx,101:160)
           imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0);
       end 
       if idx == 161
           imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',2);
       end
   end
end
%%
load('exampleDDTAUC.mat')
IRDM11 = (IRDM11./10)*100;
IRDM15 = (IRDM15./10)*100;
IRDM18 = (IRDM18./10)*100;
figure
these = linspace(0.2,0.8,size(IRDM11,1));
hold on
for ii = 1:size(IRDM11,1)
    plot(IRDM11(ii,:),'-','color',repmat(these(ii),1,3))
end
plot(mean(IRDM11,1),'-b')
ylim([0 100])
set(gca,'xtick',1:5,'ytick',0:10:100);
xlabel('Block')
ylabel('% Delay')
title('IRDM11')
%%
figure
these = linspace(0.2,0.8,size(IRDM15,1));
hold on
for ii = 1:size(IRDM15,1)
    plot(IRDM15(ii,:),'-','color',repmat(these(ii),1,3))
end
plot(mean(IRDM15,1),'-b')
ylim([0 100])
set(gca,'xtick',1:5,'ytick',0:10:100);
xlabel('Block')
ylabel('% Delay')
title('IRDM15')
%%
figure
these = linspace(0.2,0.8,size(IRDM18,1));
hold on
for ii = 1:size(IRDM18,1)
    plot(IRDM18(ii,:),'-','color',repmat(these(ii),1,3))
end
plot(mean(IRDM18,1),'-b')
ylim([0 100])
set(gca,'xtick',1:5,'ytick',0:10:100);
xlabel('Block')
ylabel('% Delay')
title('IRDM18')