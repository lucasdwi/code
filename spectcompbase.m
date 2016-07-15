function [LFPTsNaN,nNaN,indSkp,trls,clnTrls,clnEvents,relPower,stdPower,psdTrls,TFRs,fd1,fd2,avgCoh,relCoh,stdCoh] = spectcompbase(file1,dsf,thresh,onset,offset,minInt,NaNcutoff,foi,ftimwin,eventInfo,test)
%% Used to compute and plot spectrogram of normalized EEG data 
% Normalization to baseline (event2) within the same animal

%,TFR_event1,TFR_event2,statData
% Inputs: 
% file1 = data file with behavior of interest; format = 'filename' without
%   .filetype designation (in this case .mat)
%   N.B.: files need to be converted from .pl2 to .mat data structure using
%   ConvertPl2All_Files.m (by Jeffery Stott)
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
% event1 = first event of interest; format: Approach = 1; Binge = 2; Rest =
%   3. 
% event2 = second event of interest, baseline; format: Approach = 1; Binge
%   = 2; Rest = 3. 
% foi = frequencies of interest; format = [first step last] in Hz 
%   (e.g. [1 2 150] 
% ftimwin = time window steps used for both spectrograms; format = seconds 
% toi1 = time range of interest for first event; format = [first step last]
%   in seconds; N.B. last toi will usually = minInt  
% toi2 = time range of interest for first event; format = [first step last] 
%   in seconds; N.B. last toi will usually = minInt
% normType = type of normalization to conduct; format = 'absolute',
%   'relative', 'relchange' 
% 
%Code called:
% eventInd.m
% threshFilt.m
%   chkNaN.m
% trialExtract.m

% Fieldtrip code:
% ft_prepare_layout.m
% ft_freqanalysis.m
% ft_multiplot.m

% Output:
% LFPTsNaN = NaNed data structure
% nNaN = number of NaNs per channel and standard deviations from mean
% indSkp = indices of channels to be skipped over
% trl1 = trialized data for event 1
% trl2 = trialized data for event 2
% clnTrls = cell array of clean trials for the three behaviors
% TFR_event1 = output of FFT for event 1
% TFR_event2 = output of FFT for event 2
%% Check varargin - should be a n x 2 array with n = number of events of interest
if size(eventInfo,2) ~= 2
    disp('It looks like your event information is incomplete, press Ctrl+C to quit or any other key to continue');
    pause
end
%% Determine indices for events and create Rest Middle
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
tic
disp('Applying 60 Hz filter')
[b,a] = cheby1(4,0.5,[59 61]*2/adfreq,'stop');
%fvtool(b2,a2,'Fs',adfreq);
%figure; plot(LFPTs.data(1,1:20000));
LFPTsFilt = LFPTs;
for i = 1:chans
    LFPTsFilt.data(i,:) = filtfilt(b,a,LFPTsFilt.data(i,:));
end
toc
%% Downsample data
tic
if dsf > 1
    disp('Downsampling data with dwnSample.m')
    [LFPTsdwn,adfreq] = dwnSample(LFPTsFilt,dsf,adfreq,chans);
else
    disp('Skipping downsampling')
end
toc
%% Threshold data
tic
disp('Thresholding data with threshFilt')
[LFPTsNaN,nNaN,indSkp] = threshFilt(LFPTsdwn,thresh,onset,offset,minInt,NaNcutoff,adfreq,dsf,chans);
toc
%% Trialize around events
disp('Trializing data with trialExtract.m')
[clnTrls,clnEvents,trls] = trialExtract(eventInfo,eventInds,eventTs,LFPTsNaN,adfreq,minInt,dsf,chans);
%% Calculate power spectra and plot 
tic
disp('Calculating power spectra and plotting average total power with powerComp.m')
[psdTrls,relPower,stdPower,powerPlots] = powerComp(trls,adfreq,eventLabel,chans,test);
toc
%% Use Fieldtrip for Fourier Analysis
channelCombos = {'NASL', 'NASR'; 'NASL', 'NACL'; 'NASL', 'NACR'; 'NASR','NACL';'NASR','NACR';'NACL','NACR'};
% Event 1
tic
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = foi(1):foi(2):foi(3); % frequencies to use
cfg.t_ftimwin    = ones(size(cfg.foi)).*ftimwin; 
cfg.keeptrials   = 'yes';
cfg.channel      = LFPTsNaN.label;
cfg.channelcmb   = channelCombos;
cfg.toi          = eventInfo{test(1),2}(1):eventInfo{test(1),2}(2):eventInfo{test(1),2}(3);

TFRs{1} = ft_freqanalysis(cfg,trls{test(1)});
toc

if length(test) == 2
    % Event2
    tic
    cfg = []; % Create empy cfg
    cfg.output = 'powandcsd';
    cfg.channel = LFPTsNaN.label;
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.foi = foi(1):foi(2):foi(3); % Frequencies of interest
    cfg.keeptrials = 'yes'; % For stastical comparison
    cfg.t_ftimwin = ones(size(cfg.foi)).*ftimwin;
    cfg.toi = eventInfo{test(2),2}(1):eventInfo{test(2),2}(2):eventInfo{test(2),2}(3); % Times of interest
    cfg.channelcmb   = channelCombos;

    TFRs{2} = ft_freqanalysis(cfg,trls{test(2)});
    toc
end

%% Plot spectrograms - Skipped for basic analysis
% tic
% [spectroPlots] = spectroComp(trl1,trl2,TFR_event1,TFR_event2,eventLabel);
% toc
%% Plot coherence
tic
[fd1,fd2,avgCoh,relCoh,stdCoh,cohPlots] = cohComp(TFRs,eventLabel);
toc
%% Normalize
% [statData] = baselineSpectro(TFR_event1,TFR_event2,normType);
%% Save plots that exist
% Set up directory and cd to it
% mkdir(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\',file1));
% cd(strcat('C:\Users\Lucas\Desktop\GreenLab\Plots\',file1));
% % Save power plots
% if ~isempty(powerPlots)
%     powerPlotNames = {'PSDs','TotalPower'};
%     for i = 1:size(powerPlots,2)
%         savefig(powerPlots{i},powerPlotNames{i});
%     end
% end
% % Save coherence plots
% if ~isempty(cohPlots)
%     cohPlotNames = {'Coh'};
%     for i = 1:size(cohPlots,2)
%         savefig(cohPlots{i},cohPlotNames{i});
%     end
% end
% Save spectrograms
% spectroPlotNames = {'Event1','Event2'};
% for i = 1:size(spectroPlots,2)
%     savefig(spectroPlots{i},spectroPlotNames{i});
% end

%% Save variables
% Just save processed PSD and coherence data; skip actual TFR and FD data
% to save space.
save(strcat('C:\Users\Lucas\Desktop\GreenLab\data\syn\',file1,'_processed.mat'),'psdTrls','avgCoh','relCoh','stdCoh');
end