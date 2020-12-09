[data,~,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\data\'...
    'irdm\processed\base\'],{'IRDM11';'IRDM14';'IRDM15';'IRDM16';...
    'IRDM18';'IRDM21';'IRDM22';'IRDM23'},{'pow','coh'},'avg','rel');
save('C:\Users\Pythia\Documents\GreenLab\data\irdm\irdmBaseData.mat',...
    'data','files')
%%
load('C:\Users\Pythia\Documents\GreenLab\data\irdm\irdmBaseData.mat')
% Pre-surgey median split; animals 15,16,21,23 = low (0); 11,14,18,22 =
% high (1)
catData = cat(1,data{:});
x = cat(1,catData{:});
% Pre-surgery
y = [ones(9,1);zeros(14,1);ones(4,1);zeros(4,1);ones(4,1)];
cfg = lassoNetCfg([],[],'n','y','n',100,'1se',[]);
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x,y,...
    'binomial','class',1,5,1,cfg);
save('C:\Users\Pythia\Documents\GreenLab\data\irdm\basePreReal.mat',...
    'catData','x','y','accArray','allAlpha','allBeta','allLambda',...
    'cvFitsArray','hist')
% Post/Thethered; 15,16,18,21 = low (0); 11,14,22,23 = high (1); based on
% average of last 5 tethered sessions (used post-surgery value for 23)
y = [ones(9,1);zeros(18,1);ones(8,1)];
cfg = lassoNetCfg([],[],'n','y','n',100,'1se',[]);
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x,y,...
    'binomial','class',1,5,1,cfg);
save('C:\Users\Pythia\Documents\GreenLab\data\irdm\basePostReal.mat',...
    'catData','x','y','accArray','allAlpha','allBeta','allLambda',...
    'cvFitsArray','hist')
%% Pre-Surgery Prediction - Using even split permuted
permErr = [];
for ii = 1:36
    load(['C:\Users\Pythia\Documents\GreenLab\data\irdm\basePerm\evenPre\'...
        'irdmBasePerm',num2str(ii),'.mat']) 
    permErr = [permErr;allRndLambda{1}.allErr];
end
load('C:\Users\Pythia\Documents\GreenLab\data\irdm\baseReal.mat')
doubleHist((1-allLambda{1}.allErr).*100,(1-permErr).*100,'main',...
    'Pre-Surgery','xlab','Accuracy')
%% Post-Surgery Prediction - Using even split permuted
permErr = [];
for ii = 1:36
    load(['C:\Users\Pythia\Documents\GreenLab\data\irdm\basePerm\'...
        'evenPost\irdmPostBasePerm',num2str(ii),'.mat'])
    permErr = [permErr;allRndLambda{1}.allErr];
end
load('C:\Users\Pythia\Documents\GreenLab\data\irdm\basePostReal.mat',...
    'allLambda')
doubleHist((1-allLambda{1}.allErr).*100,(1-permErr).*100,'main',...
    'Post-Surgery','xlab','Accuracy')
% Post surgery rand
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\basePostRand\')
err = [];
for ii = 1:100
    load([num2str(ii),'.mat'])
    err(ii,:) = 1-allRndLambda{1}.allErr;
end
%% DDT
[data,~,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\data\'...
    'irdm\processed\ddt\'],{'IRDM11';'IRDM14';'IRDM15';'IRDM16';'IRDM18'...
    ;'IRDM21';'IRDM22'},{'pow','coh'},'avg','rel');
auc{1} = [345,325,365,355,355,360,370,375,365,350,365,365,365];
auc{2} = [215,210,225,220,180,100,180,265,225,220];
auc{3} = [95,90,60,40,75,70,55,55,70,75,85,120];
auc{4} = [295,260,280,235,180,135,170,185,160,110,210,210,180];
auc{5} = [120,135,225,200,160,155,190,185,165,145,190,195,180];
auc{6} = [95,20,55,240,100,105,90,75,80,90,25,65];
auc{7} = [80,90,120,325,175,145,120];
for ii = 1:7
    thisData = cat(1,data{ii}{:});
    thisZ = zscore(thisData);
    allZ{ii} = thisZ;
    allD{ii} = thisData;
    
    thisAucZ = zscore(auc{ii});
    aucZ{ii} = thisAucZ';
end
% allData = cat(1,allZ{:});
allData = cat(1,allD{:});
% Get rid of three recordings that do not have AUCs
allData = allData([1:9,11:46,48:64,66:end],:);

% allAucZ = cat(1,aucZ{:});
allAuc= cat(2,auc{:})';

x = allData;
% y = allAucZ;
y = allAuc;
save('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtDataRaw.mat','x','y')
%% Z-score AUC DDT
for ii = 1:100
   cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtModel\')
   load([num2str(ii),'.mat'])
   real(ii,:) = accArray{1}.acc;
   perm(ii,:) = acc;
end
doubleHist(abs(reshape(real,1,1600)),abs(reshape(perm,1,1600)),'xlab','Mean Absolute Error','main','Z-Score AUC')
doubleHist(reshape(real,1,1600),reshape(perm,1,1600),'xlab','Mean Error','main','Z-Score AUC')
%% Raw AUC DDT
for ii = 1:100
    cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtRawModel\')
    load([num2str(ii),'.mat'])
    real(ii,:) = accArray{1}.acc;
    perm(ii,:) = acc;
end
doubleHist(abs(reshape(real,1,1600)),abs(reshape(perm,1,1600)),'xlab','Mean Absolute Error','main','Raw AUC')
doubleHist(reshape(real,1,1600),reshape(perm,1,1600),'xlab','Mean Error','main','Raw AUC')
%% Median split AUC DDT
cd('G:\GreenLab\data\irdm\ddtMedianTest14\')
% cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtMedianAnimalSplit\')
for ii = 1:100
    load([num2str(ii),'.mat'])
    real(ii) = accArray{1}.acc;
    [realX(ii,:),realY(ii,:)] = perfcurve(hist.cfg.naive.testY,accArray{1}.pred,1,...
        'tvals',0:0.05:1,'usenearest',0);
    [randX(ii,:),randY(ii,:),~,randA(ii)] = perfcurve(hist.cfg.naive.testY,accArray{1}.pred(randperm(16,16)),1,'tvals',0:0.05:1,'usenearest',0);
    if ii <= 36
       perm(ii) = accArrayPerm{1}.acc;
       [permX(ii,:),permY(ii,:)] = perfcurve(histPerm.cfg.naive.testY,...
           accArrayPerm{1}.pred,1,'tvals',0:0.05:1,'usenearest',0);
    end
%     randA(ii) = fullPerm.auc;
%     randX(ii,:) = fullPerm.x;
%     randY(ii,:) = fullPerm.y;
end
% doubleHist(real,randA,'xlab','AUC','main','Median Split AUC')
figure
hold on
plot(mean(realX,1),mean(realY,1),'-k')
plot(mean(permX,1),mean(permY,1),'--k')
plot(mean(randX,1),mean(randY,1),'-.k')
%%
real = [];
perm = [];
for ii = 1:100
   cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtMedianModel\')
   load([num2str(ii),'.mat'])
   real(ii) = accArray{1}.acc;
   perm(ii) = acc;
end
doubleHist(real,perm,'xlab','AUC','main','Median Split AUC')

perm2 = [];
for ii = 1:25
    cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtMedianAnimalSplitModel\')
    load([num2str(ii),'.mat'])
    real(ii) = accArray{1}.acc;
    perm(ii) = accArrayPerm{1}.acc;
end
doubleHist(real,perm,'xlab','AUC','main','Median Split AUC')
%%
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtMedianTest\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    base5(ii) = accArray{1}.acc;
end
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtMedianTest2\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    base10(ii) = accArray{1}.acc;
end
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtMedianTest3\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    baseRnd10(ii) = accArray{1}.acc;
end
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtMedianTest4\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    baseRndNaive10(ii) = accArray{1}.acc;
end
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtMedianTest5\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    baseRndNaive5Dev(ii) = accArray{1}.acc;
end
cd('C:\Users\Pythia\Documents\GreenLab\data\irdm\ddtMedianAnimalSplitModel\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    new(ii) = accArray{1}.acc;
end
%% Timings
% Count number of tones
nTones = size(eventTs.t{7},1);
% Replace eventTs labels
eventTs.label = {'levIM','levDEL','head','feeder','ltDEL','ltIM','tone',...
    'TrialInit'};
% Preallocate structure
trial = struct('outcome',NaN,'start',NaN,'stop',NaN,'blockN',NaN,...
    'trialN',NaN,'feeder',NaN,'leverPress',NaN,'headEntry',NaN,...
    'tone_nosePoke_latency',NaN,'nosePoke_lever_latency',NaN,...
    'lever_headEntry_latency',NaN,'blockPercentDelay',NaN);
trial = repmat(trial,nTones,1);
% Set up counters
forcedCounter = 0;
blockCounter = 1;
trialCounter = 1;
% Set skipTone to false
skipTone = 0;
% Use tone as start time for trial
for ii = 1:58%nTones
    if skipTone == 0
        % Tone
        tone = eventTs.t{7}(ii);
        % Next tone
        if ii ~= nTones
            nextTone = eventTs.t{7}(ii+1);
        else
            nextTone = NaN;
        end
        % Check for double-tone
        if nextTone - tone < 30
            % Skip next tone
            skipTone = 1;
        end
        % Find nearest (subsequent) nose poke
        nearestNosePoke = eventTs.t{8}(nearest_idx3(tone,eventTs.t{8}),1);
        % Determine if nose poke occurs within 49 seconds of tone
        if nearestNosePoke >= tone+49
            trial(ii).outcome = 'timeout';
        else
            % Calculate latency to nose poke after tone
            trial(ii).tone_nosePoke_latency = nearestNosePoke - tone;
            % Get nearest feeder
            thisFeed = eventTs.t{4}(nearest_idx3(tone,eventTs.t{4},1));
            % Get next feeder - if exists, otherwise NaN
            if nearest_idx3(tone,eventTs.t{4}) ~= numel(eventTs.t{4})
                nextFeed = eventTs.t{4}(nearest_idx3(tone,eventTs.t{4})+1);
            else
                nextFeed = NaN;
            end
            % Determine if delay or immediate reward based off pellet
            % number; if no nextFeed, then immediate
            if nextFeed - thisFeed > 1 || isnan(nextFeed)
                thisReward = 'immediate';
            else
                thisReward = 'delay';
            end
            % Determine if forced or free
            if forcedCounter < 2
                % Check for failed forced trial (i.e., no feeder before
                % next tone)
                if thisFeed >= tone+49
                    trial(ii).outcome = 'failed_forced';
                else
                    trial(ii).outcome = [thisReward,'_forced'];
                    % Add one to forcedCounter
                    forcedCounter = forcedCounter + 1;
                    % Add blockN
                    trial(ii).blockN = blockCounter;
                    % Add feeder
                    trial(ii).feeder = thisFeed;
                    % Grab lever press time, using thisReward to determine
                    % lever
                    if strcmp(thisReward,'delay')
                        trial(ii).leverPress = eventTs.t{2}(...
                            nearest_idx3(tone,eventTs.t{2},1));
                    else
                        trial(ii).leverPress = eventTs.t{1}(...
                            nearest_idx3(tone,eventTs.t{1},1));
                    end
                end
            % Free
            else
                % Determine if any lever is pressed (i.e., did feeder
                % activate)
                if thisFeed >= tone+49
                    % If feeder activates more than 49 seconds after tone,
                    % fail
                    trial(ii).outcome = 'failed_free';
                else
                    % Otherwise, success
                    trial(ii).outcome = [thisReward,'_free'];
                    trial(ii).feeder = thisFeed;
                    % Grab lever press time, using thisReward to determine
                    % lever
                    if strcmp(thisReward,'delay')
                        trial(ii).leverPress = eventTs.t{2}(...
                            nearest_idx3(tone,eventTs.t{2},1));
                    else
                        trial(ii).leverPress = eventTs.t{1}(...
                            nearest_idx3(tone,eventTs.t{1},1));
                    end
                    % Calculate latency to lever press from nose poke
                    trial(ii).nosePoke_lever_latency = ...
                        trial(ii).leverPress - nearestNosePoke;
                    % Get head entry time post lever press
                    trial(ii).headEntry = eventTs.t{3}(nearest_idx3(...
                        trial(ii).leverPress,eventTs.t{3},1));
                    % Calculate latency to head entry from lever press
                    trial(ii).lever_headEntry_latency = ...
                        trial(ii).headEntry - trial(ii).leverPress;
                end
                % Set block and trial number
                trial(ii).blockN = blockCounter;
                trial(ii).trialN = trialCounter;
                % Add one to trial counter
                trialCounter = trialCounter + 1;
            end
        end
        % Get start and stop time for trial
        trial(ii).start = LFPTs.tvec(nearest_idx3(tone,LFPTs.tvec));
        % Use next tone to find stop time, unless last tone
        if ii == nTones
            trial(ii).stop = LFPTs.tvec(end);
        else
            if skipTone == 1
                trial(ii).stop = LFPTs.tvec(nearest_idx3(...
                    eventTs.t{7}(ii+2),LFPTs.tvec)-1);
            else
                trial(ii).stop = LFPTs.tvec(nearest_idx3(...
                    eventTs.t{7}(ii+1),LFPTs.tvec)-1);
            end
        end
        % If a block is finished (10 trials; trialCoutner = 11), move to
        % next block, reset trialCounter, reset forcedCounter
        if trialCounter == 11
            blockCounter = blockCounter + 1;
            trialCounter = 1;
            forcedCounter = 0;
        end
    else
        % Reset skipTone to false
        skipTone = 0;
        % Set outcome to 'doubleTone'
        trial(ii).outcome = 'doubleTone';
    end
end
% Determine block statistics - percent delay choices per block
% Convert structure to cell
trialCell = struct2cell(trial);
% Convert blockN cells to array
blocks = cell2mat(trialCell(4,:));
% Cycle through each block
for ii = 1:10
    % Find indices of this block
    inds = blocks == ii;
    % Calculate percent delay choices
    percentDelay = sum(strcmp('delay_free',trialCell(1,inds)))/sum(inds);
    % Input percent delay into each trial of this block
    for jj = logicFind(1,inds,'==')
        trial(jj).blockPercentDelay = percentDelay;
    end
end
%% Visualize session
figure
hold on
for ii = 1:8
    plot(eventTs.t{ii},ones(1,numel(eventTs.t{ii}))*ii,'.')
end
for ii = 1:numel(eventTs.t{7})
    text(eventTs.t{7}(ii),7.5,num2str(ii),'HorizontalAlignment','center')
end
set(gca,'yticklabel',eventTs.label(1:8),'ytick',(1:8))
ylim([0.5 8.5])
xlabel('Time (s)')