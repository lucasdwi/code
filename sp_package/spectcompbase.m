function [LFPTs,trls,clnTrls,clnEvents,relPower,psdTrls,coh,stdPower,stdCoh,hist] = spectcompbase(sdir,file,filter,dsf,thresh,onset,offset,minInt,foi,bands,cycles,ftimwin,overlap,cohMethod,eoi,saveParent)
% ('C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\','N7_PreStim_9_8_16','y',5,2.5,5,17000,3,[1 2 150],{'theta',[4,7];'alpha',[8,13];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]},{'full',[0 3]})
% INPUTS:
% sdir = source directory of data; format = string
% file = file name w/o extension; format = string
% filter = wether or not to use a filter; format = 'y' or 'n'
% dsf = factor with which to downfactor; format = whole integer
% thresh = threshold for detecting noise artifacts; format = mV (or other
%   y-axis scale for time series)
% onset = number of samples to NaN before noise event; format = number
%   in discrete samples, e.g. samples = sec * adfreq
% offset = number of samples to NaN after noise event; format = number
%   in discrete samples, e.g. samples = sec * adfreq
% minInt = minimum interval of clean data to use, N.B.: consider the lowest
%   frequency of interest and coherence method; format = seconds 
% foi = frequencies of interest; format = [lower step upper] in Hz
% bands = structure with bands of interest starting with lowest; format =
%   {'band1',[lower upper];'band2',[lower upper];...}
% cycles = number of cycles to use in creating frequency dependent windows
%   in ft_freqanalysis; format = number 
% ftimwin = size of window to use in ft_freqanalysis; format = number in
%   seconds; N.B.: during PowerCorr.m ftimwin will be used to compute the
%   number of cycles at the lowest frequency band of interest and will warn
%   if less than 3
% overlap = amount of overlap to use with sliding windows; format = percent
%   in decimal form
% eoi = events of interest, if one event then normalizes within that event,
%   if two then compares events; format = structure {'tag1',[0 3];'tag2',[0 3]} 
%   N.B.: if all the data use the tag 'full', otherwise use tags
%   corresponding to event markers: 'app','binge','rest','sleepOut'/'sleepIn' 
% saveParent = parent directory path to save plots and files; format =
%   string N.B.: if directory/file exists will warn about overwritting

% Examples
% DBS Parameter testing
% sdir='C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\edited\';file='N7_PreStim_9_8_16';filter='y';dsf=5;thresh=2.5;onset=5;offset=17000;minInt=5;foi=[1 2 150];bands={'theta',[4,7];'alpha',[8,13];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};cycles=3;overlap=0.5;eoi={'rest',[0 3]};
% Behavior-Ephys testing
% sdir='C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed';file='H10BaseSep27';filter='y';dsf=5;thresh=2.5;onset=5;offset=17000;minInt=3;foi=[1 2 150];bands={'theta',[4,7];'alpha',[8,13];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};cycles=3;ftimwin=[];overlap=0.5;eoi={'binge',[0 3];'rest',[0 3]};
 
%spectcompbase(sdir,file1,filter,dsf,thresh,onset,offset,minInt,foi,bands,cycles,ftimwin,eventInfo,overlap,comp)
%% Checks/Initialize
% Check eventInfo - should be a n x 2 array with n = number of events of interest
if size(eoi,2) ~= 2
    disp('It looks like you are missing either an event tag or toi in eoi, press Ctrl+C to quit or any other key to continue.');
    pause
end
% Checks for either ftimwin or cycles
if isempty(cycles) && isempty(ftimwin)
    disp('You are missing either cycles or ftimwin, press Ctrl+C to quit and re-enter variables')
    pause
end
% Checks the number of cycles at the lowest frequency band of interest
% using ftimwin; especially needed for PowerCorr.mat
if ~isempty(ftimwin)
    cycFtimwin = ftimwin*bands{1,2}(1);
    if cycFtimwin < 3
       disp('Warning: With this ftimwin, your lowest frequency band will be computed with < 3 cycles; press Ctrl+c to quit and redefine or any other key to continue...')
    pause
    end
end
% Check for cohMethod
if isempty(cohMethod) || (~strcmpi(cohMethod,'ft') && ~strcmpi(cohMethod,'mat')) 
    disp('Either cohMethod is empty or incorrectly entered, press Ctrl+c to quit and re-define')
    pause
end
% Initialize varargout with placeholders 
stdPower = []; stdCoh = [];
%% Load file
tic
cd(sdir);
disp('Loading file...')
% Load file
load(file);
chans = size(LFPTs.data,1);
toc
%% Filter
if filter == 'y'
    tic
    disp('Applying 60 Hz filter...')
    % Create first order Chebychev stopband filter from 59 to 61 Hz
    [b,a] = cheby1(2,.5,[59 61]*2/adfreq,'stop');
    % Plot filter
    %fvtool(b,a,'Fs',adfreq);
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
else
    disp('Skipping filter...')
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
% Removed nNaN and indSkp portion
[LFPTs,chk_nan] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq,dsf,chans);
%[LFPTs,nNaN,indSkp] = threshFilt(LFPTs,thresh,onset,offset,minInt,NaNcutoff,adfreq,dsf,chans);
toc

%% NaN Sleep Intervals
if exist('sleep')
    for sind = 1:size(sleep.t{1,1},1)
        LFPTs.data(:,nearest_idx2(sleep.t{1,1}(sind),LFPTs.tvec):nearest_idx2(sleep.t{1,2}(sind),LFPTs.tvec)) = NaN;
    end
end
%% Trialize around events
tic
disp('Trializing data with trialize.m')
[eventTs,eventLabel,clnTrls,clnEvents,trls] = trialize(eoi,eventTs,LFPTs,adfreq,minInt,chans);
toc
%% Calculate power spectra and plot 
tic
disp('Calculating power spectra and plotting average total power with powerComp.m')
if size(eoi,1) == 1
    [psdTrls,relPower,powerPlots] = powerComp(trls,adfreq,chans,bands,filter,foi,eoi);
end
if size(eoi,1) == 2
    [psdTrls,relPower,powerPlots,powerEventComp,stdPower] = powerComp(trls,adfreq,chans,bands,filter,foi,eoi);
end
toc
%% Create n windows (per band) and run freqanalysis for each --> powerCorr
cmb = nchoosek(1:chans,2);
for c = 1:size(cmb,1)
    channelCmb(c,:) = LFPTs.label(cmb(c,:));
end
for e = 1:size(eoi,1)
    for b = 1:size(bands,1)
        tic
        disp(['Computing CSD with ft_freqanalysis.mat for ',bands{b,1},' band...'])
        cfg              = [];
        cfg.output       = 'powandcsd';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = foi(1):foi(2):foi(3); % frequencies to use
        % Use frequency dependent windows (n cycles per window, computed at
        % start, 'cycFtimwin') with (x%) overlap
        if ~isempty(cycles)
            cfg.t_ftimwin    = ones(size(cfg.foi)).*(cycles/bands{b,2}(1));
            minTWin          = min(cfg.t_ftimwin)*overlap;
            cfg.toi          = eoi{e,2}(1):minTWin:eoi{e,2}(2);
            % Or use a constant size for windows with (x%) overlap to compute
            % cycles and apply forward
        else
            cfg.t_ftimwin    = ones(size(cfg.foi)).*(cycFtimwin/bands{b,2}(1));
            cfg.toi          = eoi{e,2}(1):ftimwin*overlap:eoi{e,2}(2);
        end
        cfg.keeptrials   = 'yes';
        cfg.channel      = LFPTs.label;
        cfg.channelcmb   = channelCmb;
        
        powCorrTFR{b} = ft_freqanalysis(cfg,trls{1,e});
        toc
    end
    % Run powerCorr
    disp('Running powerCorr.m...')
    tic
    cfg.trialwindows = 'yes';
    [STDCorr{e},MeanCorr{e},TWCorr{e},powerCorrSort{e},~,~,~,~] = powerCorr(powCorrTFR,bands,cfg);
    toc
end
%STDCorr = []; MeanCorr = []; TWCorr = []; powerCorrSort = [];
%% Use inhouse/Matlab code for coherence
if strcmpi(cohMethod,'mat')
    tic
    disp('Using cohCompMat.m for coherence...')
    % Calculate window size
    if ~isempty(ftimwin)
        winSize = (ftimwin*adfreq);
    else if ~isempty(cycles)
            % Uses lowest frequency (first in bands) to compute winSize
            winSize = cycles/bands{1,2}(1)*adfreq;
        end
    end
    [coh,cohPlots] = cohCompMat(LFPTs,chans,trls,foi,winSize,overlap,adfreq,bands,eoi);
    toc
end
%% Use Fieldtrip code for coherence
if strcmpi(cohMethod,'ft')
    tic
    disp('Using cohCompFT.mat for coherence...')
    if size(eoi,1) == 1
        [coh,cohPlots,~] = cohCompFT(LFPTs,trls,bands,chans,cycles,ftimwin,overlap,foi,eoi);
    end
    if size(eoi,1) == 2
        [coh,cohPlots,stdCoh] = cohCompFT(LFPTs,trls,bands,chans,cycles,ftimwin,overlap,foi,eoi);
        %stdCoh = varargout;
        %varargout{2} = stdCoh;
    end
    toc
end
%% Create history structure - Stores all inputs used and a few variables
hist.sdir = sdir; hist.file = file; hist.filter = filter; hist.dsf = dsf;
hist.thresh = thresh; hist.onset = onset; hist.offset = offset; 
hist.minInt = minInt; hist.foi = foi; hist.bands = bands; 
hist.cycles = cycles; hist.ftimwin = ftimwin; hist.overlap = overlap;
hist.cohMethod = cohMethod; hist.eoi = eoi; hist.saveParent = saveParent;
hist.adfreq = adfreq; hist.chk_nan = chk_nan;

if exist('bingeSize')
   hist.bingeSize = bingeSize; 
end
%% Save plots that exist
tic
disp('Saving plots...')
% Get name; file without .ext
[~,name,~] = fileparts(strcat(sdir,file));
% Set up directory and cd to it
if size(eoi,1) == 1
    mkdir(strcat(saveParent,'Plots\',name,'_',eoi{1,1}));
    cd(strcat(saveParent,'Plots\',name,'_',eoi{1,1}));
else if size(eoi,1) == 2
        mkdir(strcat(saveParent,'Plots\',name,'_',eoi{1,1},'_vs_',eoi{2,1}));
        cd(strcat(saveParent,'Plots\',name,'_',eoi{1,1},'_vs_',eoi{2,1}));
    end
end
% Save power plots
if ~isempty(powerPlots)
    powerPlotNames = {strcat(name,'PSDs'),strcat(name,'TotalPower')};
    for i = 1:size(powerPlots,2)
        savefig(powerPlots{i},powerPlotNames{i});
    end
end
% Save powerCorr plots
if ~isempty(STDCorr) && ~isempty(MeanCorr) && ~isempty(TWCorr)
    for e = 1:size(eoi,1)
        savefig(STDCorr{e},[name,'STDCorr_',eoi{e,1}]);
        savefig(MeanCorr{e},[name,'MeanCorr_',eoi{e,1}]);
        savefig(TWCorr{e},[name,'TWCorr_',eoi{e,1}]);
    end
end
% Save coherence plots
if ~isempty(cohPlots)
    cohPlotNames = {strcat(name,'Coh')};
    for i = 1:size(cohPlots,2)
        savefig(cohPlots{i},cohPlotNames{i});
    end
end
toc
%% Save variables and plots
% Just save processed PSD and coherence data; skip actual TFR data
% to save space.
tic
disp('Saving data...')
if size(eoi,1) == 1
    % Check if folder to save in exists, if not make
    if exist(saveParent,'dir') == 0
        mkdir(saveParent);
    end
    % Go to directory
    cd(saveParent)
    save(strcat(name,'_',eoi{1,1},'.mat'),'psdTrls','relPower','coh','hist','trls','LFPTs','powerCorrSort');
    % Save input variables
end
if size(eoi,1) == 2
    % Check if folder to save in exists, if not make
    if exist(saveParent,'dir') == 0
        mkdir(saveParent);
    end
    cd(saveParent)
    save(strcat(name,'_',eoi{1,1},'_vs_',eoi{2,1},'.mat'),'psdTrls','relPower','powerEventComp','stdPower','coh','hist','trls','LFPTs','powerCorrSort');
    % Save input variables
end
toc