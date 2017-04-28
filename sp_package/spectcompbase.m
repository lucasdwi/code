function [LFPTs,trls,clnTrls,clnEvents,psdTrls,coh,stdPower,stdCoh,hist] = spectcompbase(sdir,file,filter,dsf,thresh,onset,offset,foi,bands,cycles,ftimwin,overlap,cohMethod,eoi,saveParent)
%% Preproccesses data and calculates the following metrics: power, 
%  coherence, and power coupling. Read input documentation for more detail
%  on the options for these analyses and necessary information to run the
%  function. N.B. scbParamsSingle.m and scbParamsMulti.m can be used to
%  help generate the inputs need to fun this function; scbParamsMulti.m is
%  especially useful for batch processing files.
%__________________________________________________________________________
% INPUTS:
% sdir = source directory of data; format: string
% file = file name w/o extension; format: string
% filter = wether or not to use a filter; format: 'y' or 'n'
% dsf = factor with which to downfactor; format: whole integer
% thresh = threshold for detecting noise artifacts; format: mV (or other
%   y-axis scale for time series)
% onset = number of samples to NaN before noise event; format: seconds
% offset = number of samples to NaN after noise event; format: seconds
% foi = frequencies of interest; format = [lower step upper] in Hz
% bands = structure with bands of interest starting with lowest; format:
%   {'band1',[lower upper];'band2',[lower upper];...}
% cycles = number of cycles to use in creating frequency dependent windows
%   in ft_freqanalysis; format: integer 
% ftimwin = size of window to use in ft_freqanalysis; format: number in
%   seconds; N.B.: during PowerCorr.m ftimwin will be used to compute the
%   number of cycles at the lowest frequency band of interest and will warn
%   if less than 3
% overlap = amount of overlap to use with sliding windows; format: percent
%   in decimal form (1-percent; e.g. 90% overlap = 0.1)
% eoi = events of interest alongside window around that event to be
%   consider a single epoch; format: cell {'tag1',[0 3];'tag2',[-1.5 1.5]}
%   indicates to use a 3 second window around tag1 (interval) and a 3
%   second window centered at the scalar tag2.
%   N.B.: if all the data use the tag 'full', otherwise use tags
%   corresponding to event markers.
% saveParent = parent directory path to save plots and files; format:
%   string N.B.: if directory/file exists will warn about overwritting
%__________________________________________________________________________
% OUTPUTS: 
% LFPTs = local field potential structure containing the following
% information
%   type = kind of data structure; if created through Pl2tomvdmGenFile.m
%       then this will be 'tsd', i.e. "time stamped data" (See
%       Pl2tomvdmGenFile.m and tsd.m) 
%   tvec = vector of time stamps
%   data = local field potential data; chan X timestamp
%   label = channel labels
%   cfg = history structure of what had been done to data when converted
%       and processed
%   oneChan = 
% trls = 
% clnTrls = 
% clnEvents = 
% relPower = 
% psdTrls = 
% coh = 
% stdPower =
% stdCoh =
% hist = 
%__________________________________________________________________________
% USE:
% DBS Parameter testing
% sdir='C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\edited\';file='N7_PreStim_9_8_16';filter='y';dsf=5;thresh=2.5;onset=5;offset=17000;minInt=5;foi=[1 2 150];bands={'theta',[4,7];'alpha',[8,13];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};cycles=3;overlap=0.5;eoi={'rest',[0 3]};
% Behavior-Ephys testing
% sdir='C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed';file='H10BaseSep27';filter='y';dsf=5;thresh=2.5;onset=5;offset=17000;minInt=3;foi=[1 2 150];bands={'theta',[4,7];'alpha',[8,13];'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};cycles=3;ftimwin=[];overlap=0.5;eoi={'binge',[0 3];'rest',[0 3]};
 
%spectcompbase(sdir,file1,filter,dsf,thresh,onset,offset,minInt,foi,bands,cycles,ftimwin,eventInfo,overlap,comp)
%% Checks/Initialize
% Check eventInfo - should be a n x 2 array with n = number of events of
% interest
if size(eoi,2) ~= 2
    disp(['It looks like you are missing either an event tag or toi in '...
        'eoi, press Ctrl+C to quit or any other key to continue.'])
    pause
end
% Checks for either ftimwin or cycles
if isempty(cycles) && isempty(ftimwin)
    disp(['You are missing either cycles or ftimwin, press Ctrl+C to '...
        'quit and re-enter variables'])
    pause
end
% Checks the number of cycles at the lowest frequency band of interest
% using ftimwin; especially needed for PowerCorr.mat
if ~isempty(ftimwin)
    cycFtimwin = ftimwin*bands{1,2}(1);
    if cycFtimwin < 3
       disp(['Warning: With this ftimwin, your lowest frequency band '...
           'will be computed with < 3 cycles; press Ctrl+c to quit and '...
           'redefine or any other key to continue...'])
    pause
    end
end
% Check for cohMethod
if isempty(cohMethod) || (~strcmpi(cohMethod,'ft') && ~strcmpi(cohMethod,'mat')) 
    disp(['Either cohMethod is empty or incorrectly entered, press '...
        'Ctrl+c to quit and re-define'])
    pause
end
% Check bands for potential overlap
[bInd] = bandIndices(bands,foi(1):foi(2):foi(3));
if numel(bInd) ~= numel(unique(bInd))
    disp(['Warning: Band indices overlap in frequency vector! This may '...
        'lead to over-representing data.'])
    disp('Press any key to continue or Ctrl+C to quit.')
    pause
end  
% Initialize varargout with placeholders 
stdPower = []; stdCoh = [];
% Convert overlap
overlap = 1-overlap;
%% Load file
cd(sdir);
disp('Loading file...')
% Load file
load(file);
[LFPTs,chk_nan,zeroedChannel,clnTrls,clnEvents,trls,adfreq] = preProcess(LFPTs,adfreq,dsf,thresh,onset,offset,eoi,eventTs); %#ok<NODEF>
%% Calculate power spectra and plot 
tic
disp(['Calculating power spectra and plotting average total power with '...
    'powerComp.m'])
chans = size(LFPTs.data,1);
if size(eoi,1) == 1
    [psdTrls,powerPlots] = powerComp(trls,adfreq,chans,bands,filter,foi,eoi);
end
if size(eoi,1) == 2
    [psdTrls,powerPlots] = powerComp(trls,adfreq,chans,bands,filter,foi,eoi);
end
toc
%% Calculate power correlations
tic
disp('Calculating power corrleations using powerCorr.m...')
[r,rVect] = powerCorr(psdTrls); %#ok<ASGLU>
toc
%% Create n windows (per band) and run freqanalysis for each --> powerCorr
% cmb = nchoosek(1:chans,2);
% for c = 1:size(cmb,1)
%     channelCmb(c,:) = LFPTs.label(cmb(c,:));
% end
% for e = 1:size(eoi,1)
%     for b = 1:size(bands,1)
%         tic
%         disp(['Computing CSD with ft_freqanalysis.mat for ',bands{b,1},' band...'])
%         cfg              = [];
%         cfg.output       = 'powandcsd';
%         cfg.method       = 'mtmconvol';
%         cfg.taper        = 'hanning';
%         cfg.foi          = foi(1):foi(2):foi(3); % frequencies to use
%         % Use frequency dependent windows (n cycles per window, computed at
%         % start, 'cycFtimwin') with (x%) overlap
%         if ~isempty(cycles)
%             cfg.t_ftimwin    = ones(size(cfg.foi)).*(cycles/bands{b,2}(1));
%             minTWin          = min(cfg.t_ftimwin)*overlap;
%             cfg.toi          = eoi{e,2}(1):minTWin:eoi{e,2}(2);
%             % Or use a constant size for windows with (x%) overlap to compute
%             % cycles and apply forward
%         else
%             cfg.t_ftimwin    = ones(size(cfg.foi)).*(cycFtimwin/bands{b,2}(1));
%             cfg.toi          = eoi{e,2}(1):ftimwin*overlap:eoi{e,2}(2);
%         end
%         cfg.keeptrials   = 'yes';
%         cfg.channel      = LFPTs.label;
%         cfg.channelcmb   = channelCmb;
%         
%         powCorrTFR{e,b} = ft_freqanalysis(cfg,trls{1,e});
%         toc
%     end
%     % Run powerCorr
%     disp('Running powerCorr.m...')
%     tic
%     cfg.trialwindows = 'yes';
%     [STDCorr{e},MeanCorr{e},TWCorr{e},powerCorrSort{e},~,~,~,~] = powerCorr(powCorrTFR(e,:),bands,cfg);
%     toc
% end
% STDCorr = []; MeanCorr = []; TWCorr = []; powerCorrSort = [];
%% Use inhouse/Matlab code for coherence
if strcmpi(cohMethod,'mat')
    tic
    disp('Using cohCompMat.m for coherence...')
    chans = size(LFPTs.data,1);
    % Calculate window size
    if ~isempty(ftimwin)
        [~,winSize] = nearestPow2(ftimwin*adfreq);
    elseif ~isempty(cycles)
        % Uses lowest frequency (first in bands) to compute winSize
        [~,winSize] = nearestPow2((cycles/bands{1,2}(1))*adfreq);
    end
    [coh,cohPlots] = cohCompMat(LFPTs,chans,trls,foi,winSize,overlap,adfreq,bands,zeroedChannel,eoi);
    toc
end
%% Use Fieldtrip code for coherence
% if strcmpi(cohMethod,'ft')
%     tic
%     disp('Using cohCompFT.mat for coherence...')
%     if size(eoi,1) == 1
%         [coh,cohPlots,~] = cohCompFT(LFPTs,trls,bands,chans,cycles,ftimwin,overlap,foi,eoi);
%     end
%     if size(eoi,1) == 2
%         [coh,cohPlots,stdCoh] = cohCompFT(LFPTs,trls,bands,chans,cycles,ftimwin,overlap,foi,eoi);
%         %stdCoh = varargout;
%         %varargout{2} = stdCoh;
%     end
%     toc
% end
%% Create history structure - Stores all inputs used and a few variables
hist.sdir = sdir; hist.file = file; hist.filter = filter; hist.dsf = dsf;
hist.thresh = thresh; hist.onset = onset; hist.offset = offset; 
hist.foi = foi; hist.bands = bands; hist.cycles = cycles; 
hist.ftimwin = ftimwin; hist.overlap = 1-overlap; 
hist.cohMethod = cohMethod; hist.eoi = eoi; hist.saveParent = saveParent;
hist.adfreq = adfreq; hist.chk_nan = chk_nan;

if exist('bingeSize','var')
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
elseif size(eoi,1) == 2
    mkdir(strcat(saveParent,'Plots\',name,'_',eoi{1,1},'_vs_',eoi{2,1}));
    cd(strcat(saveParent,'Plots\',name,'_',eoi{1,1},'_vs_',eoi{2,1}));
end
% Save power plots
if ~isempty(powerPlots)
    powerPlotNames = {strcat(name,'PSDs'),strcat(name,'TotalPower')};
    for i = 1:size(powerPlots,2)
        savefig(powerPlots{i},powerPlotNames{i});
    end
end
% Save powerCorr plots
% if ~isempty(STDCorr) && ~isempty(MeanCorr) && ~isempty(TWCorr)
%     for e = 1:size(eoi,1)
%         savefig(STDCorr{e},[name,'STDCorr_',eoi{e,1}]);
%         savefig(MeanCorr{e},[name,'MeanCorr_',eoi{e,1}]);
%         savefig(TWCorr{e},[name,'TWCorr_',eoi{e,1}]);
%     end
% end
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
    save(strcat(name,'_',eoi{1,1},'.mat'),'psdTrls','coh','hist','trls','LFPTs','r','rVect');
    % Save input variables
end
if size(eoi,1) == 2
    % Check if folder to save in exists, if not make
    if exist(saveParent,'dir') == 0
        mkdir(saveParent);
    end
    cd(saveParent)
    save(strcat(name,'_',eoi{1,1},'_vs_',eoi{2,1},'.mat'),'psdTrls','stdPower','coh','hist','trls','LFPTs','r','rVect');
    % Save input variables
end
toc