%% NEW VERSION
% Events to boxes
% Box 1 = 25:32
% Box 2 = 33:40
% Box 3 = 41:48
% Box 4 = 59:56
inds = [25:32;33:40;41:48;49:56];
stimInds = reshape(9:24,4,4)';
dInds = reshape(1:32,8,4)';
chans = cell(1,32);
for ii = 1:32
    if ii < 10
        chans(ii) = {['FP0',num2str(ii)]};
    else
        chans(ii) = {['FP',num2str(ii)]};
    end
end
chans = reshape(chans,8,4)';

% Get processed file names; for checking data
procFiles = fileSearch('F:\irdmRound2\processedContinuous3\','.mat');
procParts = cell(1,2);
for ii = 1:numel(procFiles)
    this = strsplit(procFiles{ii},'_'); 
    procParts(ii,1:2) = [this(1),this(end-1)];
end
% Set defualt box number to 4, even if not all boxes were used
nBox = 4;
% Get file names
files = fileSearch('F:\irdmRound2\toSplit\','.mat');
% Set up array of AUCs
auc = [];
for ii = 1:size(files,2)
    disp([num2str(ii),' of ',num2str(numel(files))])
    disp(files{ii})
    load(files{ii})
    oldEventTs = eventTs;
    oldLFPTs = LFPTs;
    
    parts = strsplit(files{ii},'_');
    last = logicFind(1,cellfun(@(x) isequal(x,'DDT'),parts),'==');
    ddt = 1;
    if isempty(last)
        % Look for 'ase' to grab 'base' and 'depBase'
        last = logicFind(1,contains(parts,'ase'),'==');
        ddt = 0;
    end
    % Check if the number of LFP channels outnumbers how many would be
    % expected based on the number of animals in file name
    nChan = size(oldLFPTs.data,1);
    if nChan > (last-1)*8
        warning(['There are more LFP channels (',num2str(nChan),')'...
            ' than expected (',num2str((last-1)*8),')'])
    else
        if nChan < (last-1)*8
            warning(['There are fewer LFP channels (',num2str(nChan),')'...
                ' than expected (',num2str((last-1)*8),')'])
        end
    end
    % Set rat counter
    nRat = 1;
    % Split, going through each potential box
    for k = 1:nBox
        % First check if any data from this box exists and get inds
        dInds = contains(oldLFPTs.label,chans(k,:));
        if any(dInds)
            % Check if this is a DDT file with behavior, if so split events
            if ddt
                theseInds = [inds(k,:),stimInds(k,:)];
                eventTs.t = oldEventTs.t(theseInds);
                eventTs.label = oldEventTs.label(theseInds);
            else
                % Otherwise inset empty events
                eventTs = [];
            end
            % Split LFP data
            LFPTs.data = oldLFPTs.data(dInds,:);
            LFPTs.label = oldLFPTs.label((dInds));
            % Create new file name
            newName = strjoin([parts{nRat},parts(last:end)],'_');
            % Run ddtTrials
            if ddt
                % Check for empty feeder (indicates)
                trials = ddtTrials(eventTs,LFPTs,0);
                % Add auc to table
                auc = [auc;{newName},trials(1,1).auc];
            else
                trials = [];
            end
            % Look for already split file and compare LFPTs and eventTs to
            % make sure they were split correctly
            % Check ID and date
            thisInd = logicFind(1,contains(procParts(:,1),parts{nRat}) &...
                contains(procParts(:,2),parts{end}(1:end-4)),'==');
            if isempty(thisInd)
                warning(['no matching processed file for ',parts{nRat}])
                keyboard
            else
                % Re-name LFPTs and eventTs
                newLFPTs = LFPTs;
                newEventTs = eventTs;
                % Load split LFPTs
                load(['E:\processedMat\',...
                    procFiles{thisInd}(1:end-8),'.mat'],'LFPTs','eventTs');
                % Check first column and labels
                if newLFPTs.data(:,1) ~= LFPTs.data(:,1)
                    warning('LFP mismatch')
                    keyboard
                end
                if any(~strcmp(newLFPTs.label,LFPTs.label))
                    warning('label mismatch')
                    keyboard
                end
                % Check stim events, if any
                if ~isempty(newEventTs.t{9})
                    if ~isempty(eventTs.t{9}) && contains(parts{1},'stim')
                        if eventTs.t{9}(1) ~=  newEventTs.t{9}(1)
                            warning('stim event mismatch')
                            keyboard
                        end
                    end
                end
            end
            save(['F:\irdmRound2\toProcess\',newName],'eventTs','LFPTs',...
                'adfreq','trials')
            nRat = nRat +1;
        else
            warning('Either events or LFP missing')
        end
    end
    movefile(files{ii},'F:\irdmRound2\mat\')
end
%% manual skipper
trials = [];
save(['F:\irdmRound2\split\',newName],'eventTs','LFPTs',...
    'adfreq','trials')
k = k+1;
nRat = nRat+1;
%% OLD VERSION
oldLFPTs = LFPTs; oldEventTs = eventTs;
% %
dataCh = 1:8;
behavCh = 25:32;
% behavCh = [34,33,35:40];
stimCh = 9:10;
% stimCh = [];
eventCh = [behavCh,stimCh];

LFPTs.data = oldLFPTs.data(dataCh,:); 
LFPTs.label = oldLFPTs.label(dataCh);
eventTs.t = oldEventTs.t(eventCh);
eventTs.label = oldEventTs.label(eventCh);
trials = ddtTrials(eventTs,LFPTs,1);
% keyboard
save('D:\lsdIRDM\splitMat\IRDM35-stimCore_DDT_g_2019-12-19.mat','LFPTs','eventTs','pl2','adfreq','trials')
% save('G:\GreenLab\data\irdmNew\newSplit\IRDM27_DDT_g_2019-09-20.mat','LFPTs','eventTs','pl2','adfreq','trials')
% save('G:\GreenLab\data\irdmNew\split\LD50_alcohol_g_2019-11-06.mat','LFPTs','eventTs','pl2','adfreq')
%
dataCh = 9:16;
behavCh = 33:40;
% behavCh = [34,33,35:40];
% stimCh = [];
stimCh = 13:14;
eventCh = [behavCh,stimCh];

LFPTs.data = oldLFPTs.data(dataCh,:); 
LFPTs.label = oldLFPTs.label(dataCh);
eventTs.t = oldEventTs.t(eventCh);
eventTs.label = oldEventTs.label(eventCh);
trials = ddtTrials(eventTs,LFPTs,1);
% keyboard
% save('F:\lsdIRDM\splitMat\control\IRDM41-stimCore_DDT_g_2019-12-19.mat','LFPTs','eventTs','pl2','adfreq','trials')
% save('G:\GreenLab\data\irdmNew\newSplit\IRDM28_DDT_g_2019-09-20.mat','LFPTs','eventTs','pl2','adfreq','trials')
% save('G:\GreenLab\data\irdmNew\split\IRDM40_fedBase_g_2019-10-07.mat','LFPTs','eventTs','pl2','adfreq')
%
dataCh = 9:16;
behavCh = 41:48;
% behavCh = [42,41,43:48];
% stimCh = [];
stimCh = 17:18;
eventCh = [behavCh,stimCh];

LFPTs.data = oldLFPTs.data(dataCh,:); 
LFPTs.label = oldLFPTs.label(dataCh);
eventTs.t = oldEventTs.t(eventCh);
eventTs.label = oldEventTs.label(eventCh);
trials = ddtTrials(eventTs,LFPTs,1);
% keyboard
save('F:\lsdIRDM\toProcess\IRDM41-stimCore-LSD24_DDT_g_2019-12-17.mat','LFPTs','eventTs','pl2','adfreq','trials')
% % save('G:\GreenLab\data\irdmNew\newSplit\IRDM37_DDT_g_2019-10-29.mat','LFPTs','eventTs','pl2','adfreq','trials')
% % save('G:\GreenLab\data\irdmNew\split\IRDM32_DDT_g_2019-11-05.mat','LFPTs','eventTs','pl2','adfreq')
% % 
% dataCh = 25:32;
% behavCh = 49:56;
% % stimCh = [];
% stimCh = 21:22;
% eventCh = [behavCh,stimCh];
% 
% LFPTs.data = oldLFPTs.data(dataCh,:); 
% LFPTs.label = oldLFPTs.label(dataCh);
% eventTs.t = oldEventTs.t(eventCh);
% eventTs.label = oldEventTs.label(eventCh);
% trials = ddtTrials(eventTs,LFPTs,1);
% % keyboard
% save('F:\lsdIRDM\splitMat\control\IRDM41-stimCore_DDT_g_2020-01-16.mat','LFPTs','eventTs','pl2','adfreq','trials')
% save('G:\GreenLab\data\irdmNew\newSplit\IRDM30_DDT_g_2019-09-20.mat','LFPTs','eventTs','pl2','adfreq','trials')
% save('G:\GreenLab\data\irdmNew\split\LD38_alcohol_g_2019-11-06.mat','LFPTs','eventTs','pl2','adfreq')
%
% movefile IRDM35-stimCore_IRDM41-stimCore_IRDM32_stimCore_IRDM38_DDT_g_2019-12-19.mat D:\lsdIRDM\mat\
clear
%%
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\irdm\split\ddt\','.mat');
for ii = 1:size(files,2)
    load(files{ii})
    adfreq = 2000;
    save(files{ii},'LFPTs','eventTs','pl2','adfreq')
end