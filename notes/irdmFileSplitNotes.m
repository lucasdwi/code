%% NEW VERSION
% Uses LFP channel name to determine box; FP01 = 1, FP09 = 2, FP017 = 3,
% FP25 = 4

% Box 1 = 25:32
% Box 2 = 33:40
% Box 3 = 41:48
% Box 4 = 59:56
inds = [25:32;33:40;41:48;49:56];
stimInds = reshape(9:24,4,4)';
dInds = reshape(1:32,8,4)';
chan = {'FP01','FP09','FP17','FP25'};
files = fileSearch('F:\irdmRound2\toSplit','.mat');
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
        last = logicFind(1,cellfun(@(x) isequal(x,'Base'),parts),'==');
        ddt = 0;
    end
    nRat = last-1;
    % Set box to 1
%     nBox = 1;
    for k = 1:nRat
        nBox = logicFind(oldLFPTs.label((k-1)*8+1),chan,'==');
        % First check if this is a DDT file with behavior
        if ddt
            % Check if behavior data exists in logical box (nth box for nth
            % rat), if so skip to next box with behavior
%             while isempty(oldEventTs.t{inds(nBox,1)})
%                 nBox = nBox+1;
%             end
            theseInds = [inds(nBox,:),stimInds(k,:)];
            eventTs.t = oldEventTs.t(theseInds);
            eventTs.label = oldEventTs.label(theseInds);
        else
            eventTs = [];
        end
        % LFP data inds have no gaps, so index based off nRat
        dInds = (k-1)*8+1:(k-1)*8+8;
        LFPTs.data = oldLFPTs.data(dInds,:);
        LFPTs.label = oldLFPTs.label((dInds));
        
        
        newName = strjoin([parts{k},parts(last:end)],'_');
        % Run ddtTrials
        if ddt
            trials = ddtTrials(eventTs,LFPTs,0);
            % Add auc to table
            auc = [auc;{newName},trials(1,1).auc];
        else
            trials = [];
        end
        save(['F:\irdmRound2\split\',newName],'eventTs','LFPTs',...
            'adfreq','trials')
%         nBox = nBox+1;            
    end
    movefile(files{ii},'F:\irdmRound2\mat\')
end
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