function [LFPTs,trls,clnTrls,psdTrls,coh,stdPower,stdCoh,hist] = spectcompbase(cfg)
%% Preproccesses data and calculates the following metrics: power, 
%  coherence, and power coupling. Read input documentation for more detail
%  on the options for these analyses and necessary information to run the
%  function. N.B. scbParamsSingle.m and scbParamsMulti.m can be used to
%  help generate the inputs need to fun this function; scbParamsMulti.m is
%  especially useful for batch processing files.
%__________________________________________________________________________
% INPUTS:
% cfg = config structure with the following fields:
%   sdir = source directory of data; format: string
%   file = file name w/o extension; format: string
%   nFilt = wether or not to use a filter; format: 'y' or 'n'
%   dsf = factor with which to downfactor; format: whole integer
%   thresh = threshold for detecting noise artifacts; format: mV (or other
%       y-axis scale for time series)
%   onset = number of samples to NaN before noise event; format: seconds
%   offset = number of samples to NaN after noise event; format: seconds
%   foi = frequencies of interest; format = [lower step upper] in Hz
%   bands = structure with bands of interest starting with lowest; format:
%       {'band1',[lower upper];'band2',[lower upper];...}
%   overlap = amount of overlap to use with sliding windows; format:
%       percent in decimal form (1-percent; e.g. 90% overlap = 0.1)
%   eoi = events of interest alongside window around that event to be
%       consider a single epoch; format: cell {'tag1',[0 3];'tag2',[-1.5 
%       1.5]} indicates to use a 3 second window startint at tag1
%       (interval) and a 3 second window centered at the scalar tag2.
%       N.B.: if all the data use the tag 'all', otherwise use tags
%       corresponding to event markers.
%   vis = whether or not to plot power and coherence; format = 'y' or 'n'
%   saveParent = parent directory path to save plots and files; format:
%       string N.B.: if directory/file exists will warn about overwritting
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
% sdir='C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\edited\';
% file='N7_PreStim_9_8_16';filter='y';dsf=5;thresh=2.5;onset=5;
% offset=17000;foi=[1 2 150];bands={'theta',[4,7];'alpha',[8,13];
% 'beta',[15,30];'lgam',[45,65];'hgam',[70,90]};overlap=0.5;
% eoi={'rest',[0 3]};

% Behavior-Ephys testing
% sdir='C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed';
% file='H10BaseSep27';filter='y';dsf=5;thresh=2.5;onset=5;offset=17000;
% foi=[1 2 150];bands={'theta',[4,7];'alpha',[8,13];'beta',[15,30];
% 'lgam',[45,65];'hgam',[70,90]};overlap=0.5;
% eoi={'binge',[0 3];'rest',[0 3]};
%% Unpack cfg
sdir = cfg.sdir;
file = cfg.file;
nFilt = cfg.nFilt;
dsf = cfg.dsf;
thresh = cfg.thresh;
onset = cfg.onset;
offset = cfg.offset;
foi = cfg.foi;
bands = cfg.bands;
overlap = cfg.overlap;
cohMethod = cfg.cohMethod;
skip = cfg.skip;
eoi = cfg.eoi;
vis = cfg.vis;
saveParent = cfg.saveParent;
%% Checks/Initialize
% Check eventInfo - should be a n x 2 array with n = number of events of
% interest
if size(eoi,2) ~= 2
    disp(['It looks like you are missing either an event tag or toi in '...
        'eoi, press Ctrl+C to quit or any other key to continue.'])
    pause
end
% Check for cohMethod
if isempty(cohMethod) || (~strcmpi(cohMethod,'mtm') && ~strcmpi(cohMethod,'mat')) 
    error(['Either cohMethod is empty or incorrectly entered, press '...
        'Ctrl+C to quit and re-define'])
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
% Load file - just those variables to be used
load(file,'LFPTs','eventTs','adfreq');
% Skip channels
skipInd = ~ismember(1:size(LFPTs.data,1),skip);
LFPTs.data = LFPTs.data(skipInd,:);
LFPTs.label = LFPTs.label(skipInd);
[LFPTs,chk_nan,zeroedChannel,clnTrls,trls,adfreq] = preProcess(LFPTs,...
    adfreq,dsf,thresh,onset,offset,eoi,eventTs); %#ok<NODEF>
%% Calculate power spectra and plot 
tic
disp(['Calculating power spectra and plotting average total power with '...
    'powerComp.m'])
[psdTrls,powerPlots] = powerComp(trls,adfreq,bands,nFilt,foi,eoi,vis);
toc
%% Calculate power correlations - requires at least 2 trials otherwise 
% gives NaNs
tic
disp('Calculating power corrleations using powerCorr.m...')
[r,rVect] = powerCorr(psdTrls); %#ok<ASGLU>
toc
%% Calculate coherence
if strcmpi(cohMethod,'mat')
    tic
    disp('Using mscohere.m to calculate coherence...')
    [coh,cohPlots] = cohComp(trls,adfreq,eoi,bands,zeroedChannel,foi,...
        nFilt,vis,cohMethod,'overlap',overlap);
    toc
elseif strcmpi(cohMethod,'mtm')
    tic
    disp('Using cmtm.m to calculate coherence...')
    [coh,cohPlots] = cohComp(trls,adfreq,eoi,bands,zeroedChannel,foi,...
        nFilt,vis,cohMethod,'NW',8);
    toc
end
%% Create history structure - Stores all inputs used and a few variables
hist.sdir = sdir; hist.file = file; hist.filter = nFilt; hist.dsf = dsf;
hist.thresh = thresh; hist.onset = onset; hist.offset = offset; 
hist.foi = foi; hist.bands = bands; hist.overlap = 1-overlap; 
hist.cohMethod = cohMethod; hist.eoi = eoi; hist.saveParent = saveParent; 
hist.adfreq = adfreq; hist.chk_nan = chk_nan; hist.eventTs = eventTs;

if exist('bingeSize','var')
   hist.bingeSize = bingeSize; 
end
%% Save plots that exist
disp('Saving plots...')
% Get name; file without .ext
[~,name,~] = fileparts(strcat(sdir,file));
% Set up directory and cd to it
if size(eoi,1) == 1 || size(eoi,1) > 2
    mkdir(strcat(saveParent,'Plots\',name,'_',eoi{1,1}));
    cd(strcat(saveParent,'Plots\',name,'_',eoi{1,1}));
elseif size(eoi,1) == 2
    mkdir(strcat(saveParent,'Plots\',name,'_',eoi{1,1},'_vs_',eoi{2,1}));
    cd(strcat(saveParent,'Plots\',name,'_',eoi{1,1},'_vs_',eoi{2,1}));
end
% Save power plots
if ~isempty(powerPlots)
    powerPlotNames = {strcat(name,'PSDs'),strcat(name,'TotalPower')};
    for ii = 1:size(powerPlots,2)
        savefig(powerPlots{ii},powerPlotNames{ii});
    end
end
% Save coherence plots
if ~isempty(cohPlots)
    cohPlotNames = {strcat(name,'Coh'),strcat(name,'NormCoh')};
    for ii = 1:size(cohPlots,2)
        savefig(cohPlots{ii},cohPlotNames{ii});
    end
end
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
    save(strcat(name,'_',eoi{1,1},'.mat'),'psdTrls','coh','hist','trls',...
        'LFPTs','r','rVect');
    % Save input variables
end
if size(eoi,1) == 2
    % Check if folder to save in exists, if not make
    if exist(saveParent,'dir') == 0
        mkdir(saveParent);
    end
    cd(saveParent)
    save(strcat(name,'_',eoi{1,1},'_vs_',eoi{2,1},'.mat'),'psdTrls',...
        'stdPower','coh','hist','trls','LFPTs','r','rVect');
    % Save input variables
end
if size(eoi,1) > 2
    % Check if folder to save in exists, if not make
    if exist(saveParent,'dir') == 0
        mkdir(saveParent);
    end
    cd(saveParent)
    save(strcat(name,'_',eoi{1,1},'.mat'),'psdTrls','stdPower','coh',...
        'hist','trls','LFPTs','r','rVect');
end
toc