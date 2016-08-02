function [LFPTsNaN,nNaN,indSkp,trls,clnTrls,clnEvents,relPower,psdTrls,TFRs,fds,avgCoh,relCoh,stdPower,stdCoh] = spectcompbase(file1,filter,dsf,thresh,onset,offset,minInt,NaNcutoff,foi,ftimwin,eventInfo,comp)
%% Used to compute and plot spectrogram of normalized EEG data 
% Normalization to baseline (event2) within the same animal

% Inputs: 
% file1 = data file with behavior of interest; format = 'filename' without
%   .filetype designation (in this case .mat)
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
% ftimwin = time window steps used for both spectrograms; format = seconds
% eventInfo = structure of information about events; format: row = event;
%   column 1 = event tag (Approach = 1; Binge = 2; Rest = 3); column 2 = time
%   range of interest [first step last] in seconds; N.B. last toi will
%   usually = minInt 
% comp = events to analyze; format = integer or integer-pair of event tags
%   (Approach = 1; Binge = 2; Rest = 3) 

% Data structure used throughout obtained from ConvertPl2All_Files.m and
% Pl2tomvdmGenFile.m (Jeffery Stott)

% Code called:
% eventInd.m
% threshFilt.m
%   chkNaN.m
% trialExtract.m

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
%% Initialize varargout
stdPower = []; stdCoh = [];
%% Checks
% Check eventInfo - should be a n x 2 array with n = number of events of interest
if size(eventInfo,2) ~= 2
    disp('It looks like you are missing either an event tag or toi in eventInfo, press Ctrl+C to quit or any other key to continue.');
    pause
end
% Check last toi of eventInfo = minInt
for e = 1:size(eventInfo,1)
    if eventInfo{e,2}(3) ~= minInt
        disp(strcat('Event ',num2str(e),' toi has a different length than minInt, press Ctrl+C to quit or any other key to continue.'));
        pause
    end
end
%% Determine indices for events and create 'Rest Middle'
tic
cd('C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed\');
disp('Loading file and extracting event indices with eventInd.m')
% Adds .mat to file name for loading
load(strcat(file1,'.mat'));
% Counts number of channels in data
chans = size(LFPTs.data,1); 
[eventInds,eventTs] = eventInd(eventTs);
% Initialize event labels
behaviors = {'Approach','Binge','Rest'};
for ii = 1:size(eventInfo,1)
    eventLabel{ii} = behaviors{eventInfo{ii,1}};
end
toc
%% Filter
if filter == 'y'
    tic
    disp('Applying 60 Hz filter')
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
    toc
end
%% Downsample data
tic
if dsf > 1
    disp('Downsampling data with dwnSample.m')
    [LFPTs,adfreq] = dwnSample(LFPTs,dsf,adfreq,chans);
else
    disp('Skipping downsampling')
end
toc
%% Threshold data
tic
disp('Thresholding data with threshFilt')
[LFPTs,nNaN,indSkp] = threshFilt(LFPTs,thresh,onset,offset,minInt,NaNcutoff,adfreq,dsf,chans);
toc
%% Check channel quality
% tic
%
% toc
%% Trialize around events
disp('Trializing data with trialExtract.m')
[clnTrls,clnEvents,trls] = trialExtract(eventInfo,eventInds,eventTs,LFPTs,adfreq,minInt,dsf,chans);
%% Calculate power spectra and plot 
tic
disp('Calculating power spectra and plotting average total power with powerComp.m')
if length(comp) == 1
    [psdTrls,relPower,powerPlots] = powerComp(trls,adfreq,eventLabel,chans,comp);
end
if length(comp) == 2
    [psdTrls,relPower,powerPlots,varargout] = powerComp(trls,adfreq,eventLabel,chans,comp);
    stdPower = varargout; clear varargout;
end
toc
%% Use Fieldtrip for Fourier Analysis
% Get channel combinations
cmb = nchoosek(1:chans,2);
channelCmb = zeros(size(cmb,1));
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
cfg.t_ftimwin    = ones(size(cfg.foi)).*ftimwin; 
cfg.keeptrials   = 'yes';
cfg.channel      = LFPTs.label;
cfg.channelcmb   = channelCmb;
cfg.toi          = eventInfo{comp(1),2}(1):eventInfo{comp(1),2}(2):eventInfo{comp(1),2}(3);

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
    cfg.t_ftimwin = ones(size(cfg.foi)).*ftimwin;
    cfg.toi = eventInfo{comp(2),2}(1):eventInfo{comp(2),2}(2):eventInfo{comp(2),2}(3); % Times of interest
    cfg.channelcmb   = channelCmb;

    TFRs{2} = ft_freqanalysis(cfg,trls{comp(2)});
    toc
end

%% Plot spectrograms - Skipped for basic analysis
% tic
% [spectroPlots] = spectroComp(trl1,trl2,TFR_event1,TFR_event2,eventLabel);
% toc
%% Plot coherence
tic
if length(comp) == 1
    [avgCoh,relCoh,cohPlots,fds] = cohComp(TFRs,eventLabel,comp);
end
if length(comp) == 2
    [avgCoh,relCoh,cohPlots,fds,varargout] = cohComp(TFRs,eventLabel,comp);
    stdCoh = varargout;
    %varargout{2} = stdCoh;
end
toc
%% Normalize
% [statData] = baselineSpectro(TFR_event1,TFR_event2,normType);
%% Save plots that exist
% Set up directory and cd to it
cond = {'_approach','_binge','_rest','_vs'};
if length(comp) == 1
    mkdir(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\',file1,cond{comp}));
    cd(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\',file1,cond{comp}));
else if length(comp) == 2
        mkdir(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\',file1,cond{comp(1)},'_vs',cond{comp(2)}));
        cd(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\',file1,cond{comp(1)},'_vs',cond{comp(2)}));
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
%% Save spectrograms
% spectroPlotNames = {'Event1','Event2'};
% for i = 1:size(spectroPlots,2)
%     savefig(spectroPlots{i},spectroPlotNames{i});
% end

%% Save variabless
% Just save processed PSD and coherence data; skip actual TFR and FD data
% to save space.
if length(comp) == 1
   save(strcat('C:\Users\Lucas\Desktop\GreenLab\data\processed\',file1,cond{comp},'.mat'),'psdTrls','relPower','avgCoh','relCoh','fds'); 
end
if length(comp) == 2
    save(strcat('C:\Users\Lucas\Desktop\GreenLab\data\processed\',file1,cond{comp(1)},'_vs',cond{comp(2)},'.mat'),'psdTrls','relPower','avgCoh','relCoh','stdCoh','fds');
end