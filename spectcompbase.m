function [LFPTsNaN,nNaN,indSkp,trls,clnTrls,clnEvents,relPower,psdTrls,TFRs,fds,avgCoh,relCoh,stdPower,stdCoh] = spectcompbase(file1,dsf,thresh,onset,offset,minInt,NaNcutoff,foi,ftimwin,eventInfo,comp)
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
% foi = frequencies of interest; format = [first step last] in Hz 
%   (e.g. [1 2 150] 
% ftimwin = time window steps used for both spectrograms; format = seconds
% eventInfo = structure of information about events; format: row = event;
%   column 1 = event tag (Approach = 1; Binge = 2; Rest = 3); column 2 = time
%   range of interest [first step last] in seconds; N.B. last toi will
%   usually = minInt 
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
%% Initialize varargout
stdPower = []; stdCoh = [];
%% Check eventInfo - should be a n x 2 array with n = number of events of interest
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
if length(comp) == 1
    [psdTrls,relPower,powerPlots] = powerComp(trls,adfreq,eventLabel,chans,comp);
    %varargout = {};
end
if length(comp) == 2
    [psdTrls,relPower,powerPlots,varargout] = powerComp(trls,adfreq,eventLabel,chans,comp);
    stdPower = varargout; clear varargout;
    %varargout{1} = stdPower;
end
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
cfg.toi          = eventInfo{comp(1),2}(1):eventInfo{comp(1),2}(2):eventInfo{comp(1),2}(3);

TFRs{1} = ft_freqanalysis(cfg,trls{comp(1)});
toc

if length(comp) == 2
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
    cfg.toi = eventInfo{comp(2),2}(1):eventInfo{comp(2),2}(2):eventInfo{comp(2),2}(3); % Times of interest
    cfg.channelcmb   = channelCombos;

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