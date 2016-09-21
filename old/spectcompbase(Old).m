function [LFPTs,trls,clnTrls,clnEvents,relPower,psdTrls,TFRs,fds,avgCoh,relCoh,stdPower,stdCoh] = spectcompbase(sdir,file1,filter,dsf,thresh,onset,offset,minInt,foi,bands,cycles,ftimwin,eventInfo,overlap,comp)
%function [LFPTs,nNaN,indSkp,trls,clnTrls,clnEvents,relPower,psdTrls,TFRs,fds,avgCoh,relCoh,stdPower,stdCoh] = spectcompbase(sdir,file1,filter,dsf,thresh,onset,offset,minInt,NaNcutoff,foi,cycles,eventInfo,comp)
%% Used to compute and plot spectrogram of normalized EEG data 
% Normalization to baseline (event2) within the same animal


%% NOTE: COMMENTED OUT CHKNAN FROM THRESHFILT
%%
% Inputs: 
% sdir = source directory path; format = string
% file1 = data file with behavior of interest; format = 'filename'
%   N.B.: files need to be converted from .pl2 to .mat data structure using
%   ConvertPl2All_Files.m (by Jeffery Stott)
% filter = whether or not to apply filter; format = 'y' or 'n' with
%   apostrophes
% dsf = downsample factor; format = integer
%   N.B.: dsf, onset, and offset need to be chosen such that onset and
%   offset are multiples of dsf
% thresh = threshold above/below which data is NaNed; format = mV 
% onset = amount of data before threshold cross to NaN; format = indices 
%   (e.g. 1 = 0.5 ms; 2000 = 1 sec) 
% offset = amount of data after threshold cross to NaN; format = indices
%   (e.g. 1 = 0.5 ms; 2000 = 1 sec)  
% minInt = minimum interval length of data to keep; format = seconds 
% NaNcutoff = number of standard deviations from the average number of NaNs
%   per channel at which to get rid of channel; format = integer (e.g. 1.5) 
% foi = frequencies of interest; format = [first step last] in Hz 
%   (e.g. [1 2 150] 
% bands = definitions of frequency band; format = {'name', [start,stop];...}
% cycles = number of cycles for each wavelet
% ftimwin = size of time window for analysis
% eventInfo = structure of information about events; format: row = event;
%   column 1 = event tag (Approach = 1; Binge = 2; Rest = 3); column 2 =
%   vector [start:end]
% overlap = percent overlap for windows; format = decimal form (.5 for 50%)
% comp = events to analyze; format = integer or integer-pair of event tags
%   in [] (Approach = 1; Binge = 2; Rest = 3)

% Data structure used throughout obtained from ConvertPl2All_Files.m and
% Pl2tomvdmGenFile.m (Jeffery Stott)

% Code called:
% eventInd.m
% dwnSample.m
% threshFilt.m
%   chkNaN.m
% trialExtract.m
% powerComp.m
% cohComp.m

% Fieldtrip code:
% ft_prepare_layout.m
% ft_connectivityanalysis.m
% ft_freqanalysis.m
% ft_multiplot.m

% Output:
% LFPTsNaN = NaNed data structure
% nNaN = number of NaNs per channel and standard deviations from mean
% indSkp = indices of channels to be skipped over
% trls = structure of trialized data
% clnTrls = cell array of clean trials for the three behaviors
% clnEvents = collapsed clean trials
% relPower = matrix of relative power levels per frequency band
%   if one event: percent of total power from theta to h gamma per band
%   it two events: change in event1 power per band from event2
% psdTrls = outputs from pwelch in dB
% TFRs = outputs from ft_freqanalysis
% fds = outputs from ft_connenctivity analysis
% avgCoh = average coherence in each channel combination per band
% relCoh = matrix of relative coherence levels per band
%   if one event: percent of total coherence from theta to h gamma per band
%   if two events: change in event1 coherence per band from event2
% stdPower = the variance in power (only if two events)
% stdCoh = the variance in coherence (only if two events)

% Example:
% [...] = spectcompbase('C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed\','H10BaseSep27','y',5,2.5,5,17000,3,1.5,[1 2 150],0.5,{1,[0 0.005 3];2,[0 0.005 3];3,[0 0.005 3]},[3])
% 
% sdir='C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed\';file1='H10BaseSep27';dsf=5;thresh=2.5;onset=5;offset=17000;minInt=5;foi=[1 2 150];eventInfo={1,[0 3];2,[0 3];3,[0 3]};comp=[3];filter='y';overlap=0.5;cycles=3;bands = {'theta',[4,7];'alpha',[8,13];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};
%% Initialize varargout
stdPower = []; stdCoh = [];
%% Checks
% Check eventInfo - should be a n x 2 array with n = number of events of interest
if size(eventInfo,2) ~= 2
    disp('It looks like you are missing either an event tag or toi in eventInfo, press Ctrl+C to quit or any other key to continue.');
    pause
end
% Check last toi of eventInfo = minInt
% for e = 1:size(eventInfo,1)
%     if eventInfo{e,2}(3) ~= minInt
%         disp(strcat('Event ',num2str(e),' toi has a different length than minInt, press Ctrl+C to quit or any other key to continue.'));
%         pause
%     end
% end
%% Load file
tic
cd(sdir);
disp('Loading file and extracting event indices with eventInd.m')
% Load file
load(file1);
toc
%% Determine indices for events and create 'Rest Middle'
tic
% Counts number of channels in data
chans = size(LFPTs.data,1); 
[eventInds,eventTs,markers] = eventInd(eventTs);
% Initialize event labels
behaviors = {'Approach','Binge','Rest'};
for ii = 1:size(eventInfo,1)
    eventLabel{ii} = behaviors{eventInfo{ii,1}};
end
toc
%% Filter
if filter == 'y'
    tic
    disp('Applying 60 Hz filter...')
    % Create first order Chebychev stopband filter from 59 to 61 Hz
    [b,a] = cheby1(4,0.5,[59 61]*2/adfreq,'stop');
    % Plot filter
    %fvtool(b2,a2,'Fs',adfreq);
    %figure; plot(LFPTs.data(1,1:20000));
    % Setup filtered LFPT structure
    LFPTsFilt = LFPTs;
    % Apply filter to all chans
    for i = 1:chans
        LFPTsFilt.data(i,:) = filtfilt(b,a,LFPTsFilt.data(i,:));
    end
    % Overwrite LFPTs with LFPTsFilt
    LFPTs = LFPTsFilt;
else
    disp('Skipping filter...')
    toc
end
%% Downsample data
tic
if dsf > 1
    disp('Downsampling data with dwnSample.m')
    [LFPTs,adfreq] = dwnSample(LFPTs,dsf,adfreq,chans);
else
    disp('Skipping downsampling...')
end
toc
%% Threshold data
tic
disp('Thresholding data with threshFilt...')
% Got rid of nNaN and indSkp portion
[LFPTs] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq,dsf,chans);
%[LFPTs,nNaN,indSkp] = threshFilt(LFPTs,thresh,onset,offset,minInt,NaNcutoff,adfreq,dsf,chans);
toc

%% Trialize around events
disp('Trializing data with trialExtract.m')
[clnTrls,clnEvents,trls] = trialExtract(eventInfo,eventInds,eventTs,LFPTs,adfreq,minInt,dsf,chans);
%% Calculate power spectra and plot 
tic
disp('Calculating power spectra and plotting average total power with powerComp.m')
if length(comp) == 1
    [psdTrls,relPower,powerPlots] = powerComp(trls,adfreq,eventLabel,chans,comp,bands);
end
if length(comp) == 2
    [psdTrls,relPower,powerPlots,varargout] = powerComp(trls,adfreq,eventLabel,chans,comp,bands);
    stdPower = varargout; clear varargout;
end
toc
%% Use Fieldtrip for Fourier Analysis
% Get channel combinations
cmb = nchoosek(1:chans,2);
for c = 1:size(cmb,1)
    channelCmb(c,:) = LFPTs.label(cmb(c,:));
end
% Event 1
tic
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = foi(1):foi(2):foi(3); % frequencies to use
cfg.t_ftimwin    = cycles./cfg.foi; 
%cfg.t_ftimwin    = ones(size(cfg.foi)).*ftimwin; 
cfg.keeptrials   = 'yes';
cfg.channel      = LFPTs.label;
cfg.channelcmb   = channelCmb;
minTWin = min(cfg.t_ftimwin)*overlap;
cfg.toi          = eventInfo{comp(1),2}(1):minTWin:eventInfo{comp(1),2}(2);
%cfg.toi          = eventInfo{comp(1),2}(1):ftimwin*overlap:eventInfo{comp(1),2}(2);


TFRs{1} = ft_freqanalysis(cfg,trls{comp(1)});
toc

if length(comp) == 2
    % Event2
    tic
    cfg = []; % Create empy cfg
    cfg.output = 'powandcsd';
    cfg.channel = LFPTs.label;
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.foi = foi(1):foi(2):foi(3); % Frequencies of interest
    cfg.keeptrials = 'yes'; % For stastical comparison
    %cfg.t_ftimwin = cycles./cfg.foi; 
    cfg.t_ftimwin = ones(size(cfg.foi)).*ftimwin; 
    %minTWin = min(cfg.t_ftimwin)*overlap;
    %cfg.toi = eventInfo{comp(1),2}(1):minTWin:eventInfo{comp(1),2}(2);
    cfg.toi = eventInfo{comp(1),2}(1):ftimwin*overlap:eventInfo{comp(1),2}(2);
    cfg.channelcmb = channelCmb;

    TFRs{2} = ft_freqanalysis(cfg,trls{comp(2)});
    toc
end

%% Compute power correlations - Michael Connerney
% tic
% disp('Running powerCorr...')
% [STDCorr,MeanCorr,powerCorr] = powerCorr(TFRs,foi);
% toc
%% Plot spectrograms - Skipped for basic analysis
% tic
% [spectroPlots] = spectroComp(trl1,trl2,TFR_event1,TFR_event2,eventLabel);
% toc
%% Plot coherence
tic
if length(comp) == 1
    [avgCoh,relCoh,cohPlots,fds] = cohComp(TFRs,eventLabel,comp,bands);
end
if length(comp) == 2
    [avgCoh,relCoh,cohPlots,fds,varargout] = cohComp(TFRs,eventLabel,comp,bands);
    stdCoh = varargout;
    %varargout{2} = stdCoh;
end
toc
%% Cross Frequency Coupling - Alexander Nakhnikian
% tic
% % Concatenate clean data
% disp('Concatenating data...')
% dataCat = [];
% for t = 1:length(trls{1,comp}.trial)
%     dataCat = horzcat(dataCat,trls{1,comp}.trial{1,t});
% end
% toc
% %% Run gmwMI
% tic
% disp('Computing MIs...')
% [MIs,fSpace,Hvals,probVects,statsData] = gmwMI(dataCat,adfreq,foi(1),foi(3),18,length(trls{1,comp}.trial),minInt*adfreq+1,3,5,4);
% toc
% %% Run gmwCFCpermFDR
% tic
% disp('Running stats on MIs...')
% [sigMIs,h,pvals,critP,adjP] = gmwCFCpermFDR(MIs,statsData,18,20,0.05);
% toc
% %% Plot
% tic
% disp('Plotting MIs...')
% for p = 1:chans
%     for a = 1:chans
%         plotMIs(squeeze(MIs(:,:,p,a)),adfreq,fSpace,4,0)
%         title(strcat('Phase Channel',num2str(p),' vs Amp Channel',num2str(a)));
%     end
% end
% toc
%% Normalize
% [statData] = baselineSpectro(TFR_event1,TFR_event2,normType);
%% Save plots that exist
tic
disp('Saving plots...')
% Get name; file1 without .ext
[~,name,~] = fileparts(strcat(sdir,file1));
% Set up directory and cd to it
cond = {'approach','binge','rest'};
if length(comp) == 1
    mkdir(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\test2',name,'_',cond{comp}));
    cd(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\',name,'_',cond{comp}));
else if length(comp) == 2
        mkdir(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\',name,'_',cond{comp(1)},'_vs_',cond{comp(2)}));
        cd(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\',name,'_',cond{comp(1)},'_vs_',cond{comp(2)}));
    end
end
% Save power plots
if ~isempty(powerPlots)
    powerPlotNames = {'PSDs','TotalPower'};
    for i = 1:size(powerPlots,2)
        savefig(powerPlots{i},powerPlotNames{i});
    end
end
% Save coherence plots
if ~isempty(cohPlots)
    cohPlotNames = {'Coh'};
    for i = 1:size(cohPlots,2)
        savefig(cohPlots{i},cohPlotNames{i});
    end
end
toc
%% Save spectrograms
% spectroPlotNames = {'Event1','Event2'};
% for i = 1:size(spectroPlots,2)
%     savefig(spectroPlots{i},spectroPlotNames{i});
% end

%% Save variables
% Just save processed PSD and coherence data; skip actual TFR data
% to save space.
tic
disp('Saving data...')
if length(comp) == 1
    % Check if folder to save in exists, if not make
    if exist(strcat('C:\Users\Lucas\Desktop\GreenLab\data\',cond{comp}),'file') ~= 7
        mkdir(strcat('C:\Users\Lucas\Desktop\GreenLab\data\',cond{comp}));
    end
    save(strcat('C:\Users\Lucas\Desktop\GreenLab\data\',cond{comp},'\',name,'_',cond{comp},'test2','.mat'),'psdTrls','relPower','avgCoh','relCoh','fds','bingeSize');
end
if length(comp) == 2
    % Check if folder to save in exists, if not make
    if exist(strcat('C:\Users\Lucas\Desktop\GreenLab\data\',cond{comp(1)},'_vs_',cond{comp(2)}),'file') ~= 7
        mkdir(strcat('C:\Users\Lucas\Desktop\GreenLab\data\',cond{comp(1)},'_vs_',cond{comp(2)}));
    end
    save(strcat('C:\Users\Lucas\Desktop\GreenLab\data\',name,'_',cond{comp(1)},'_vs_',cond{comp(2)},'.mat'),'psdTrls','relPower','avgCoh','relCoh','stdCoh','fds','bingeSize');
end
toc