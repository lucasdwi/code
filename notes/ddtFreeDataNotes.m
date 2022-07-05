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
data = []; c = 0; aC = 1;
% Set window size (sec) desired and minimum number of clean data for use
winSize = 10;
notNan = 1;
ti = 1;
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
                for jj = 1:size(this,2)
                    % Subtract 1.5 seconds from lever press to account for
                    % feature times given as center of 3-second window
                    % (i.e., window centered at 1.5 includes data from [0
                    % 3])
                    if strcmp('immediate_free',this(1,jj)) || strcmp(...
                            'delay_free',this(1,jj))
                        % Use leverPress as the last possible point to
                        % include
                        stop = nearest_idx3(cell2mat(this(6,jj)),...
                            psdTrls{1}.t+1.5,-1);
                        % Subtract winSize+1.5 to get all possible data
                        % within the time of interest
                        start = nearest_idx3(cell2mat(this(6,jj))-winSize,...
                            psdTrls{1}.t-1.5,1);
                        % Check for NaNs; only keep if >= notNan (using
                        % first values in relPower as proxy)
                        if sum(~isnan(psdTrls{1}.relPow(1,1,1,start:...
                                stop)))>=notNan
                            % Add to counter
                            c = c+1;
                            % Get starting index for this session
                            if first == 1
                                animalInds(aC,1) = c;
                            end
                            % data
                            data(c).animal = parts{1};
                            data(c).date = parts{end-1};
                                                        
%                             data(c).rawPower = psdTrls{1}.Pow(:,:,:,...
%                                 start:stop);
                            data(c).relPower = psdTrls{1}.relPow(:,:,:,...
                                start:stop);
                            
%                             data(c).rawCoh = coh{1}.Cxy(:,:,:,start:stop);
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
                            % Save data window
                            data(c).start = psdTrls{1}.t(start);
                            data(c).stop = psdTrls{1}.t(stop);
                            % time
                            data(c).winCenter = cell2mat(this(6,jj))-...
                                psdTrls{1}.t(start:stop);
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
% save('F:/irdmRound2/trialData.mat','animalInds','data')
%% collate all times
% allTimes = [];
% clean = [];
% for ii = 1:numel(data)
%     % Check NaNs
%     cleanInds = ~isnan(data(ii).relPower(1,1,1,:));
%     clean(ii) = sum(cleanInds);
%     % Account for the times being center points by adding and subtracting
%     % 1.5 seconds
%    allTimes = [allTimes,data(ii).winCenter(cleanInds)-1.5,...
%        data(ii).winCenter(cleanInds)+1.5];
% end
% figure
% subplot(1,2,1)
% histogram(allTimes)
% box off
% xlabel('seconds before lever (sec)')
% ylabel('number of trials')
% subplot(1,2,2)
% [f,x] = ecdf(allTimes);
% plot(x,f.*sum(clean))
% box off
% xlabel('seconds before lever (sec)')
% ylabel('number of trials')
%% get number of trials as we allow time to go further back
% load('F:/irdmRound2/trialData.mat')
c=1;
cleanTrials = cell(1,8);
for x = [3.1,4:10]
    thresh = 0:0.1:x;
%     cleanTrials = zeros(numel(data),numel(thresh));
    for ii = 1:numel(data)
        % Check NaNs
        cleanInds = ~isnan(data(ii).relPower(1,1,1,:));
        clean{c}(ii) = sum(cleanInds);
        cleanCenters = data(ii).winCenter(cleanInds);
        % Check if any window contains thresh but does not exceed max thresh
        for jj = 1:numel(thresh)
            if any(thresh(jj) >= cleanCenters-1.5 & ...
                    thresh(jj) <= cleanCenters+1.5 & ...
                    thresh(end) >= cleanCenters+1.5)
                cleanTrials{c}(ii,jj) = 1;
            else
                cleanTrials{c}(ii,jj) = 0;
            end
        end
    end
    c = c+1;
end
cumSum = cell(1,8);
for c = 1:8
    for ii = 1:size(cleanTrials{c},2)
        cumSum{c}(ii) = sum(any(cleanTrials{c}(:,1:ii),2));
    end
end
%% plots how many trials can be gotten max 
figure
plot(3:10,cellfun(@max,cumSum),'.')
figure
thresh = [3.1,4:10];
hold on
for ii = 1:8
    plot(0:0.1:thresh(ii),sum(cleanTrials{ii},1))
end
%%
figure
subplot(1,2,1)
imagesc(cleanTrials{3})
set(gca,'xtick',0:9:51,'xticklabel',0:5)
xlabel('time before lever press (sec)')
ylabel('trial')
title('clean data per trial')
colormap gray

subplot(1,2,2)
plot(0:0.1:5,(sum(cleanTrials{3})/cumSum{3}(end))*100)
xlabel('time before lever press (sec)')
ylabel('% trials')
title('% of trials with included time')
%% plot data and behavior
figure
hold on
% plot(LFPTs.tvec,LFPTs.data(1,:))
plot(psdTrls{1}.t,squeeze(psdTrls{1}.Overall(1,1,1,:)),'ok')
for ii = 1:size(times,1)
    plot([psdTrls{1}.t(times(ii,1))-1.5 psdTrls{1}.t(times(ii,2))+1.5],[2 2],'-r')
end
%% Plot histogram of gap between data start and lever press
figure
% Account for the 1.5 seconds subtracted earlier
histogram(extractfield(data,'stop')-extractfield(data,'leverPress')+1.5);
figure
histogram(extractfield(data,'start')-extractfield(data,'leverPress')+1.5);
%% Immediate and delay lever press
% Find when session switches by comparing animal ID and date
allData = squeeze(struct2cell(data));
delay = []; immediate = [];
% Find minimal number of delays and immediates each session contributes
for ii = 1:size(animalInds,1)
    % Find indices of delays and immediates
    delayInd = logicFind(1,cell2mat(allData(5,animalInds(ii,1):...
        animalInds(ii,2))),'==');
    immediateInd = logicFind(0,cell2mat(allData(5,animalInds(ii,1):...
        animalInds(ii,2))),'==');
    if ii == 1
        start = 0;
    else
        start = animalInds(ii-1,2);
    end
%     % Count number of delays and immediates
%     dN(ii) = numel(delayInd);
%     iN(ii) = numel(immediateInd);
%     % Find min
%     thisMin = min([dN(ii),iN(ii)]);
%     % Pull that number of both trials
%     delay = cat(2,delay,allData(:,start+delayInd(randperm(dN(ii),...
%         thisMin))));
%     immediate = cat(2,immediate,allData(:,start+immediateInd(randperm(...
%         iN(ii),thisMin))));
    % Grab all delays and immediates
    delay = cat(2,delay,allData(:,start+delayInd));
    immediate = cat(2,immediate,allData(:,start+immediateInd));
end
% Count number of trials per animal (just counts delays, so half of real
% number)
u = unique(delay(1,:));
for ii = 1:numel(u)
    nD(ii) = sum(strcmp(u(ii),delay(1,:)));
    nI(ii) = sum(strcmp(u(ii),immediate(1,:)));
end
% save('F:/irdmRound2/delayImmediate.mat','delay','immediate','u','nD','nI')
%% Get trials per animal (not intervention)
% baseInd = ~contains(u,'-');
baseInd = 1:79; % Includes intervention
uBase = u(baseInd);
% Only use an animal if it has at least 20 total samples (immediate and
% delay)
theseInd = nD(baseInd)+nI(baseInd) >= 20;
[catDelay,catImmediate] = deal(cell(1,sum(theseInd)));
for ii = 1:sum(theseInd)
    thisAnimalDeInd = strcmp(delay(1,:),uBase{ii});
    thisAnimalImInd = strcmp(immediate(1,:),uBase{ii});
    % Delay
    delayRelPower = cellfun(@(x) reshape(squeeze(mean(x,4,'omitnan')),1,...
        48),delay(3,thisAnimalDeInd),'UniformOutput',0);
    catDelayPower = cat(1,delayRelPower{:});
    delayNormCoh = cellfun(@(x) reshape(permute(squeeze(mean(x,4,...
        'omitnan')),[2,1,3]),1,168),delay(4,thisAnimalDeInd),...
        'UniformOutput',0);
    catDelayCoh = cat(1,delayNormCoh{:});
    catDelay{ii} = [catDelayPower,catDelayCoh];
    % Immediate
    immediateRelPower = cellfun(@(x) reshape(squeeze(mean(x,4,'omitnan')...
        ),1,48),immediate(3,thisAnimalImInd),'UniformOutput',0);
    catImmediatePower = cat(1,immediateRelPower{:});
    immediateNormCoh = cellfun(@(x) reshape(permute(squeeze(mean(x,4,...
        'omitnan')),[2,1,3]),1,168),immediate(4,thisAnimalImInd),...
        'UniformOutput',0);
    catImmediateCoh = cat(1,immediateNormCoh{:});
    catImmediate{ii} = [catImmediatePower,catImmediateCoh];
end
%%
all = logicFind(1,theseInd,'==');
% Calculate weights
baseDelay = nD(baseInd);%cellfun(@(x) size(x,1),catDelay);
baseIm = cellfun(@(x) size(x,1),catImmediate);
delayW = 1./(baseDelay/sum(baseDelay));
imW = 1./(baseIm/sum(baseIm));
% Propogate weights
for ii = 1: numel(all)
    dWeights{ii} = repmat(delayW(ii),baseDelay(ii),1);
    iWeights{ii} = repmat(imW(ii),baseIm(ii),1);
end
% LOO
for ii = 1:numel(all)
    % Leave out this animal and build a model from all others
    others = ~ismember(all,all(ii));
    testXLOO{ii} = [catDelay{ii};catImmediate{ii}];
    testYLOO{ii} = [zeros(size(catDelay{ii},1),1);...
        ones(size(catImmediate{ii},1),1)];
    
    trainXLOO{ii} = [cat(1,catDelay{others});cat(1,catImmediate{others})];
    trainYLOO{ii} = [zeros(size(cat(1,catDelay{others}),1),1);...
        ones(size(cat(1,catImmediate{others}),1),1)];
    % Grab weights
    trainWLOO{ii} = [cat(1,dWeights{others});cat(1,iWeights{others})];
    cfg = lassoNetCfg({testXLOO{ii},testYLOO{ii}},[],'n','n','n',100,...
        '1se',trainWLOO{ii});
    [~,~,~,~,accLOO,histLOO] = lassoNet(trainXLOO{ii},trainYLOO{ii},...
        'binomial','deviance',1,10,1,cfg);
end
% 80/20
xAll = [cat(1,catDelay{:});cat(1,catImmediate{:})];
yAll = [zeros(size(cat(1,catDelay{:}),1),1);...
    ones(size(cat(1,catImmediate{:}),1),1)];
wAll = [cat(1,dWeights{:});cat(1,iWeights{:})];
rng(n)
[trainXind,trainYind,testXind,testYind] = trainTest((1:size(xAll,1))',...
    (1:size(yAll,1))',0.2);
cfg = lassoNetCfg({xAll(testXind,:),yAll(testYind)},[],'n','n','n',...
    100,'1se',wAll(trainXind));
[~,~,~,~,accAll,histAll] = lassoNet(xAll(trainXind,:),yAll(trainYind),...
    'binomial','deviance',1,10,1,cfg);
%% Plot above (delay vs. immmediate) results
for ii = 1:100
    load(['F:/irdmRound2/delayImmediate/delayIm',num2str(ii),'.mat']')
    allA(ii) = accAll{1}.acc;
    allAX(ii,:) = interp1(linspace(0,1,numel(accAll{1}.x)),accAll{1}.x,...
        linspace(0,1,600));
    allAY(ii,:) = interp1(linspace(0,1,numel(accAll{1}.y)),accAll{1}.y,...
        linspace(0,1,600));
    for jj = 1:size(accLOO,2)
        if ~isempty(accLOO{jj})
            looA(ii,jj) = accLOO{jj}{1}.acc;
            % Convert all ROCs to length 100 for plotting convenience
            looAX(ii,jj,:) = interp1(linspace(0,1,...
                numel(accLOO{jj}{1}.x)),accLOO{jj}{1}.x,linspace(0,1,100));
            looAY(ii,jj,:) = interp1(linspace(0,1,...
                numel(accLOO{jj}{1}.y)),accLOO{jj}{1}.y,linspace(0,1,100));
        else
            looA(ii,jj) = NaN;
            [looAX(ii,jj,:),looAY(ii,jj,:)] = deal(NaN);
        end
    end
    glmA(ii) = glmAllA;
end

figure
hold on
for ii = 1:size(looAX,2)
   plot(squeeze(mean(looAX(:,ii,:),1)),squeeze(mean(looAY(:,ii,:),1)))
end
title('delay vs. immediate: each LOO')
xlabel('FPR')
ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)

figure
hold on
plot(squeeze(mean(allAX,1)),squeeze(mean(allAY,1)))
plot(squeeze(mean(mean(looAX,1,'omitnan'),2,'omitnan')),...
    squeeze(mean(mean(looAY,1,'omitnan'),2,'omitnan')))
legend({['80:20 \mu = ',num2str(round(mean(allA),2)),'\pm',...
    num2str(round(conf(allA,0.95),3))],['LOO \mu = ',...
    num2str(round(mean(mean(looA,'omitnan'),'omitnan'),2)),'\pm',...
    num2str(round(conf(mean(looA(:,[1:16,18:end])),0.95),2))]},...
    'location','nw')
box off
title('delay vs. immediate')
xlabel('FPR')
ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
%% Decisive split - decisive = >80% one lever or the other, indecisive = 
% <40% one lever or the other
load('F:/irdmRound2/trialData.mat')
u = unique(extractfield(data,'animal'));
for ii = 1:size(animalInds,1)
    delayPerc(ii) = data(animalInds(ii,1)).blockD;
end
allDelay = extractfield(data,'blockD');
topInd = logicFind(1,allDelay<=0.2 | allDelay>=0.8,'==');
midInd = logicFind(1,allDelay<=0.6 & allDelay>=0.4,'==');
% Remove data with missing channels

top = squeeze(struct2cell(data(topInd)));%squeeze(struct2cell(data(topInd(randperm(numel(topInd),sum(midInd))))));
mid = squeeze(struct2cell(data(midInd)));

% Get indices of each animal
uT = unique(top(1,:));
uTLog = cellfun(@(x) strcmp(top(1,:),x),uT,'UniformOutput',0);
uM = unique(mid(1,:));
uMLog = cellfun(@(x) strcmp(mid(1,:),x),uM,'UniformOutput',0);
% Get base
uBaseLog = ~contains(u,'-');
uBaseInd = logicFind(1,uBaseLog,'==');
% Get number of trials per animal
% save('F:/irdmRound2/topMid.mat','top','mid','uBaseInd','u')
%% Plot top (im. vs delay) and mid (im vs. delay)
for ii = 1:100
    load(['F:/irdmRound2/decisive/decisiveImDe',num2str(ii),'.mat'])
    for jj = 1:numel(accTopLOO)
        if iscell(accTopLOO{jj})
            topLOO(ii,jj) = accTopLOO{jj}{1}.acc;
            topGLMLOO(ii,jj) = glmTopLOOA(jj);
            topLOOX(ii,jj,:) = interp1(linspace(0,1,numel(accTopLOO{jj}{1}.x)),...
                accTopLOO{jj}{1}.x,linspace(0,1,400));
            topLOOY(ii,jj,:) = interp1(linspace(0,1,numel(accTopLOO{jj}{1}.y)),...
                accTopLOO{jj}{1}.y,linspace(0,1,400));
        else
            topLOO(ii,jj) = NaN;
            topGLMLOO(ii,jj) = NaN;
            topLOOX(ii,jj,:) = ones(1,400).*NaN;
            topLOOY(ii,jj,:) = ones(1,400).*NaN;
        end
    end
    for jj = 1:numel(accMidLOO)
        midLOO(ii,jj) = accMidLOO{jj}{1}.acc;
        midGLMLOO(ii,jj) = glmMidLOOA(jj);
        midLOOX(ii,jj,:) = interp1(linspace(0,1,numel(accMidLOO{jj}{1}.x)),...
            accMidLOO{jj}{1}.x,linspace(0,1,400));
        midLOOY(ii,jj,:) = interp1(linspace(0,1,numel(accMidLOO{jj}{1}.y)),...
            accMidLOO{jj}{1}.y,linspace(0,1,400));
    end
    allTop(ii) = accTopAll{1}.acc;
    allTopX(ii,:) = interp1(linspace(0,1,numel(accTopAll{1}.x)),...
        accTopAll{1}.x,linspace(0,1,400));
    allTopY(ii,:) = interp1(linspace(0,1,numel(accTopAll{1}.y)),...
        accTopAll{1}.y,linspace(0,1,400));
    allMid(ii) = accMidAll{1}.acc;
    allMidX(ii,:) = interp1(linspace(0,1,numel(accMidAll{1}.x)),...
        accMidAll{1}.x,linspace(0,1,400));
    allMidY(ii,:) = interp1(linspace(0,1,numel(accMidAll{1}.y)),...
        accMidAll{1}.y,linspace(0,1,400));
end
subplot(2,3,1)
hold on
plot(mean(allTopX),mean(allTopY))
plot(mean(allMidX),mean(allMidY))
title('immediate vs. delay (top and mid)')
xlabel('FPR'); ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['top: \mu = ',num2str(mean(allTop)),'\pm',num2str(conf(allTop,0.95))],...
    ['mid: \mu = ',num2str(mean(allMid)),'\pm',num2str(conf(allMid,0.95))]})

subplot(2,3,2)
hold on
plot(squeeze(mean(mean(topLOOX,1,'omitnan'),2,'omitnan')),...
    squeeze(mean(mean(topLOOY,1,'omitnan'),2,'omitnan')))
plot(squeeze(mean(mean(midLOOX,1,'omitnan'),2,'omitnan')),...
    squeeze(mean(mean(midLOOY,1,'omitnan'),2,'omitnan')))
legend({['top: \mu = ',num2str(mean(mean(topLOO,'omitnan'),'omitnan')),...
    '\pm',num2str(conf(squeeze(mean(topLOO(:,[1:14,16]),1)),0.95))],...
    ['mid: \mu = ',num2str(mean(mean(midLOO,'omitnan'),'omitnan')),...
    '\pm',num2str(conf(squeeze(mean(midLOO,1)),0.95))]})
title('top and mid: immediate vs. delay (LOO)')
xlabel('FPR'); ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)

subplot(2,3,3)
hold on
for ii = 1:16
    plot(squeeze(mean(topLOOX(:,ii,:),1)),squeeze(mean(topLOOY(:,ii,:),1)))
end
title('immediate vs. delay (top): each LOO')
xlabel('FPR'); ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)

subplot(2,3,4)
hold on
for ii = 1:4
    plot(squeeze(mean(midLOOX(:,ii,:),1)),squeeze(mean(midLOOY(:,ii,:),1)))
end
title('immediate vs. delay (mid): each LOO')
xlabel('FPR'); ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
% Plot top vs. mid
for ii = 1:100
    load(['F:/irdmRound2/topMid/topMid',num2str(ii),'.mat'])
    for jj = 1:20
        looLassoA(ii,jj) = accLOO{jj}{1}.acc;
        looLassoX(ii,jj,:) = interp1(linspace(0,1,numel(accLOO{jj}{1}.x)),...
            accLOO{jj}{1}.x,linspace(0,1,24));
        looLassoY(ii,jj,:) = interp1(linspace(0,1,numel(accLOO{jj}{1}.y)),...
            accLOO{jj}{1}.y,linspace(0,1,24));
    end
    allA(ii) = accArray{1}.acc;
    if ~isnan(allA(ii))
        allX(ii,:) = interp1(linspace(0,1,numel(accArray{1}.x)),...
            accArray{1}.x,linspace(0,1,24));
        allY(ii,:) = interp1(linspace(0,1,numel(accArray{1}.y)),...
            accArray{1}.y,linspace(0,1,24));
    else
        allX(ii,:) = deal(NaN);
        allY(ii,:) = deal(NaN);
    end
end
subplot(2,3,5)
hold on
plot(mean(allX,'omitnan'),mean(allY,'omitnan'))
plot(squeeze(mean(mean(looLassoX,1),2)),squeeze(mean(mean(looLassoY,1),2)))
xlabel('FPR'); ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
legend({['80:20 \mu = ',num2str(round(mean(allA),2)),'\pm',...
    num2str(round(conf(allA,0.95),2))],['LOO \mu = ',...
    num2str(round(mean(mean(looLassoA,1),2),2)),'\pm',...
    num2str(round(conf(mean(looLassoA,1),0.95),2))]})

subplot(2,3,6)
hold on
for ii = 1:20
    plot(squeeze(mean(looLassoX(:,ii,:),1)),squeeze(mean(looLassoY(:,ii,:),1)))
end
xlabel('FPR'); ylabel('TPR')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
%% Plot immediate vs. delay (individual animal models)
for ii = 1:100
    load(['F:/irdmRound2/delayImmediateInd/delayImInd',num2str(ii),'.mat'])
    allA(:,:,ii) = a;
end
%%
% Get animal ID
load('F:/irdmRound2/delayImmediate.mat','u','nI','nD')
baseIndLog = ~contains(u,'-');
baseInd = logicFind(1,~contains(u,'-'),'==');
uBase = u(baseIndLog);
theseIndLog = nD(baseIndLog)+nI(baseIndLog) >= 40;
theseInd = logicFind(1,theseIndLog,'==');
theseID = uBase(theseInd);
figure
thisM = mean(allA,3,'omitnan');
pcolor(padarray(thisM([2,7,8,9,15,17,10,16,18,3,5,11,13,1,4,6,12,14],...
    [2,7,8,9,15,17,10,16,18,3,5,11,13,1,4,6,12,14]),[1,1],'post'))
caxis([0.33 0.8])
colormap viridis
set(gca,'xtick',0.5:18.5,'ytick',0.5:18.5)
xlabel('train')
ylabel('test')
set(gca,'xtick',1.5:19.5,'xticklabel',...
    theseID([2,7,8,9,15,17,10,16,18,3,5,11,13,1,4,6,12,14]),...
    'ytick',1.5:19.5,'yticklabel',...
    theseID([2,7,8,9,15,17,10,16,18,3,5,11,13,1,4,6,12,14]))
xtickangle(45)
colorbar

% Get sex; 1 = male
sex = [zeros(1,9),ones(1,9)];
% Get left or right; left = 1
leftRight = [ones(1,6),0,0,0,1,1,1,1,0,0,0,0,0];
% Set data to use; where dim = test x train
this = thisM([2,7,8,9,15,17,10,16,18,3,5,11,13,1,4,6,12,14],...
    [2,7,8,9,15,17,10,16,18,3,5,11,13,1,4,6,12,14]);
% male->female
mf = this(sex==0,sex==1);
% male->male
mm = this(sex==1,sex==1);
mm(1:size(mm,1)+1:end) = NaN;
% female->female
ff = this(sex==0,sex==0);
ff(1:size(ff,1)+1:end) = NaN;
% female->male
fm = this(sex==1,sex==0);

totalMM = sum(sum(~isnan(mm)));
totalFF = sum(sum(~isnan(ff)));
totalMF = sum(sum(~isnan(mf)));
totalFM = sum(sum(~isnan(fm)));
p = 0.7;
nMM = sum(sum(mm>p));
nFF = sum(sum(ff>p));
nMF = sum(sum(mf>p));
nFM = sum(sum(fm>p));

% left->left
ll = this(leftRight==1,leftRight==1);
ll(1:size(ll,1)+1:end) = NaN;
% left->right
lr = this(leftRight==0,leftRight==1);
% right->left 
rl = this(leftRight==1,leftRight==0);
% right->right
rr = this(leftRight==0,leftRight==0);
rr(1:size(rr,1)+1:end) = NaN;

totalLL = sum(sum(~isnan(ll)));
totalRR = sum(sum(~isnan(rr)));
totalLR = sum(sum(~isnan(lr)));
totalRL = sum(sum(~isnan(rl)));

nLL = sum(sum(ll>p));
nRR = sum(sum(rr>p));
nLR = sum(sum(lr>p));
nRL = sum(sum(rl>p));

perc(4,:) = [nMM/totalMM,nMF/totalMF,nFF/totalFF,nFM/totalFM,...
    nLL/totalLL,nLR/totalLR,nRR/totalRR,nRL/totalRL];
n(4,:) = [nMM,nMF,nFF,nFM,nLL,nLR,nRR,nRL];
total(4,:) = [totalMM,totalMF,totalFF,totalFM,totalLL,totalLR,totalRR,...
    totalRL];
%% Plot cross performance (above) vs. change in trapz AUC
% Get IDs of animals with intervention and immediate vs delay model
these = theseInter(contains(theseInter,theseBase));
baseAUC = [];
inter = {'sIL','sNAcC','guan0.3','mph3','mph0.3','mph1','mph0.1',...
    'sNAcC-LSD','sIL-LSD'};
for ii = 1:numel(these)
    for jj = 1:numel(inter)
        % Get baseline mean trapz AUC
        
    end
end
%% Plot cross performance (above) vs. effect
thisInter = contains(theseInter,theseBase);
thisBase = contains(theseBase,theseInter);
thisBaseInd = logicFind(1,thisBase,'==');
interI = 1;
for ii = [1,2,4]
    c = 1;
    
        figure
    for jj = 1:15
        subplot(3,5,c)
        plot(thisEffect(thisInter,ii),mean(allA(thisBaseInd(jj),thisBaseInd,:),3),'o')
%         thisEffect(thisInter,jj) == 1;
%         
        [~,p(interI,jj)] = ttest2(thisEffect(thisInter,ii),mean(allA(thisBaseInd(jj),thisBaseInd,:),3));
        c = c+1;
    end
    interI = interI+1;
end
%%
figure
set(gcf,'renderer','Painters')
hold onclear
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
allTimes = linspace(0.2,0.8,size(IRDM11,1));
hold on
for ii = 1:size(IRDM11,1)
    plot(IRDM11(ii,:),'-','color',repmat(allTimes(ii),1,3))
end
plot(mean(IRDM11,1),'-b')
ylim([0 100])
set(gca,'xtick',1:5,'ytick',0:10:100);
xlabel('Block')
ylabel('% Delay')
title('IRDM11')
%%
figure
allTimes = linspace(0.2,0.8,size(IRDM15,1));
hold on
for ii = 1:size(IRDM15,1)
    plot(IRDM15(ii,:),'-','color',repmat(allTimes(ii),1,3))
end
plot(mean(IRDM15,1),'-b')
ylim([0 100])
set(gca,'xtick',1:5,'ytick',0:10:100);
xlabel('Block')
ylabel('% Delay')
title('IRDM15')
%%
figure
allTimes = linspace(0.2,0.8,size(IRDM18,1));
hold on
for ii = 1:size(IRDM18,1)
    plot(IRDM18(ii,:),'-','color',repmat(allTimes(ii),1,3))
end
plot(mean(IRDM18,1),'-b')
ylim([0 100])
set(gca,'xtick',1:5,'ytick',0:10:100);
xlabel('Block')
ylabel('% Delay')
title('IRDM18')