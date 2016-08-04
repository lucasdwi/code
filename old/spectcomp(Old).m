%% Program to compute and compare spectrograms based around separate events or the same event but differnet animal or condition
function spectcomp(file1,file2,event1,event2,twin1,twin2,foi,ftimwin,toi1,toi2)
% file1 = first data file to be used; format: 'filename'
% file2 = second data file to be used; format: 'filename'
    % can be the same or different from file1
% event1 = first event of interest; format: (orientation = 1,rest start = 2, rest middle = 8, rest end = 3,approach start = 4,approach end = 5,binge start = 6,binge end = 7)
% event2 = second event of interest; format: (orientation = 1,rest start = 2, rest middle = 8, rest end = 3,approach start = 4,approach end = 5,binge start = 6,binge end = 7)
% twin1 = time window around event 1; format: [before after] in seconds
% twin2 = time window around event 2; format: [before after] in seconds
    % make sure that the abosulte amount of time for both time windows is equal if not identical before and after indices
% foi = frequencies of interest; format: [first step last] in Hz 
    % the same range will be used for both events
% ftimwin = time window steps used for both spectrograms; format: seconds
% toi1 = time range of interest for first event; format: [first step last] in seconds
% toi2 = time range of interest for second event; format: [first step last] in seconds


% E.G. spectcomp('conditionxanimal1','conditionyanimal1',6,6,[-8 2],[-8 2],[1 2 150],0.5,[-8 0.05 2],[-8 0.05 2])
% This would compare the 8 seconds before and 2 seconds after ([-8 2]) the start of
% a binge (9) between condition x and y of animal 1. The spectrograms would
% be constructed with 1 to 150 Hz in steps of 2 Hz [1 2 150] using 0.5 second
% windows, between 8 seconds before the start of the both binges and 2 seconds
% after at 0.05 second intervals [-8 0.05 2].

% spectcomp('I2FoodDep24Dec16_pl2done_plx.pl2','I2FoodDep24Dec16_pl2done_plx.pl2',6,8,[-8 2],[-5 5],[1 2 150],0.5,[-8 0.05 2],[-5 0.05 5])
% The above entry will reproduce the figures from project2.m
%% First check if .pl2 files need to be converted into .m usable data structures
chk_data = input('Do the data files need to be converted? [Y/N]','s')
if strcmpi(chk_data,'y')
    [fd] = ConvertPl2All_Files % Calls Pl2tomvdm.m to convert all .pl2 data files
end
%% Load file1
load(file1)
%% Detemine Indices for Events
[ind] = eventInd(eventTs);

%% Determine Indices for Events
% Necessary because not all data files have the same indices of events in
% eventTs
indO = find(not(cellfun('isempty',strfind(eventTs.label,'Orientation'))));
indRs = find(not(cellfun('isempty',strfind(eventTs.label,'Rest (Start)'))));
indRe = find(not(cellfun('isempty',strfind(eventTs.label,'Rest (End)'))));
indAs = find(not(cellfun('isempty',strfind(eventTs.label,'Approach (Start)'))));
indAe = find(not(cellfun('isempty',strfind(eventTs.label,'Approach (End)'))));
indBs = find(not(cellfun('isempty',strfind(eventTs.label,'Binge (Start)'))));
indBe = find(not(cellfun('isempty',strfind(eventTs.label,'Binge (End)'))));
indRm = 11; 
% Set up index matrix 
ind = [];
ind(1:8,1) = [indO,indRs,indRe,indAs,indAe,indBs,indBe,indRm];
%% Create New evnetTs, Rest Middle
eventTs.label{11} = 'Rest (Middle)';
eventTs.t{11} = eventTs.t{1,indRs} + (eventTs.t{1,indRe} - eventTs.t{1,indRs})/2;
%% Trialize Data Around Event 1
% Create Trials
cfg = [];
cfg.t = eventTs.t{1,ind(event1)}; % Timestamps corresponding to event 1
cfg.twin = twin1; % time window around event 1

clear trl; % Reset trl
trl(:,3) = nearest_idx3(cfg.t,LFPTs.tvec);
trl(:,2) = nearest_idx3(cfg.t+cfg.twin(2),LFPTs.tvec);
trl(:,1) = nearest_idx3(cfg.t+cfg.twin(1),LFPTs.tvec);

trl(:,3) = trl(:,1) - trl(:,3);
%% Extract Event1 Data
data_trl1 = []; %set up empty trialized data structure
data_trl1.label = LFPTs.label;
data_trl1.fsample = adfreq;
data_trl1.trial = {}; % Set up empty trial structure
% Fill data_trl.trial structure with trialized data series
i = 1;
for i=1:size(trl,1)
    %data_trl1.trial{1,i} = LFPTs.data(1,trl(i,1):trl(i,2)) % Just first
    %channel
    data_trl1.trial{1,i} = LFPTs.data(1:4,trl(i,1):trl(i,2)) % All four channels, or the first four channels in the case of a data file with >4 channels
    end
%% Create and fill data_trl.time structure
data_trl1.time = {};
twin = [cfg.twin(1):0.0005:cfg.twin(2)];

i = 1;
for i=1:size(trl,1)
    data_trl1.time{1,i} = twin; 
end
% Create and fill data_trl.sampleinfo
data_trl1.sampleinfo = [];
i = 1;
for i=1:size(trl,1)
    data_trl1.sampleinfo(i,1) = trl(i,1);
    data_trl1.sampleinfo(i,2) = trl(i,2);
end

% Define layout for later plotting
cfg = [];
cfg.layout = 'ordered'; cfg.channel = LFPTs.label; cfg.showlabels = 'yes';
layout = ft_prepare_layout(cfg,data_trl1)
%% Use Fieldtrip to construct spectograms
cfg = []; % Create empy cfg
cfg.output = 'pow';
cfg.channel = LFPTs.label;
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = foi(1):foi(2):foi(3); % Frequencies of interest
cfg.keeptrials = 'yes'; % For stastical comparison
cfg.t_ftimwin = ones(size(cfg.foi)).*ftimwin;
cfg.toi = toi1(1):toi1(2):toi1(3); % Times of interest

TFR_event1 = ft_freqanalysis(cfg,data_trl1);

% Plot
figure
cfg = []; cfg.channel = LFPTs.label; cfg.layout = layout;
ft_multiplotTFR(cfg,TFR_event1)

%% Trialize Data Around Event 2
% Do the same type of analysis as before, but for the second event

%% First run Matt and Wilder's code, PL2tomvdm.m, to convert pl2 data into
% usable data structures
% check if file 1 and file 2 are the same or different
tf = strcmp(file1,file2);
if tf == 0
    load(file2)
% Determine Indices for Events
% Necessary because not all data files have the same indices of events in
% eventTs
    indO = find(not(cellfun('isempty',strfind(eventTs.label,'Orientation'))));
    indRs = find(not(cellfun('isempty',strfind(eventTs.label,'Rest (Start)'))));
    indRe = find(not(cellfun('isempty',strfind(eventTs.label,'Rest (End)'))));
    indAs = find(not(cellfun('isempty',strfind(eventTs.label,'Approach (Start)'))));
    indAe = find(not(cellfun('isempty',strfind(eventTs.label,'Approach (End)'))));
    indBs = find(not(cellfun('isempty',strfind(eventTs.label,'Binge (Start)'))));
    indBe = find(not(cellfun('isempty',strfind(eventTs.label,'Binge (End)'))));
    indRm = 11; 
% Set up index matrix 
    ind = [];
    ind(1:8,1) = [indO,indRs,indRe,indAs,indAe,indBs,indBe,indRm];
% Create New evnetTs, Rest Middle
    eventTs.label{11} = 'Rest (Middle)';
    eventTs.t{11} = eventTs.t{1,indRs} + (eventTs.t{1,indRe} - eventTs.t{1,indRs})/2;
end
%% Trialize Data Around Event 2
% Create Trials
cfg = [];
cfg.t = eventTs.t{1,ind(event2)}; % Timestamps corresponding to event 2
cfg.twin = twin1; % time window around event 2

clear trl; % Reset trl
trl(:,3) = nearest_idx3(cfg.t,LFPTs.tvec);
trl(:,2) = nearest_idx3(cfg.t+cfg.twin(2),LFPTs.tvec);
trl(:,1) = nearest_idx3(cfg.t+cfg.twin(1),LFPTs.tvec);

trl(:,3) = trl(:,1) - trl(:,3);
%% Extract Event 2 Data
data_trl2 = []; %set up empty trialized data structure
data_trl2.label = LFPTs.label;
data_trl2.fsample = adfreq;
data_trl2.trial = {}; % Set up empty trial structure
% Fill data_trl.trial structure with trialized data series
i = 1;
for i=1:size(trl,1)
    % data_trlr.trial{1,i} = LFPTs.data(1,trl(i,1):trl(i,2)) % Just looking at first channel
    data_trl2.trial{1,i} = LFPTs.data(1:4,trl(i,1):trl(i,2)) % All four channels
end
% Create and fill data_trl.time structure
data_trl2.time = {};
twin = [cfg.twin(1):0.0005:cfg.twin(2)];
i = 1;
for i=1:size(trl,1)
    data_trl2.time{1,i} = twin;
end
% Create and fill data_trl.sampleinfo
data_trl2.sampleinfo = [];
i = 1;
for i=1:size(trl,1)
    data_trl2.sampleinfo(i,1) = trl(i,1);
    data_trl2.sampleinfo(i,2) = trl(i,2);
end
% Define layout for later plotting
cfg = [];
cfg.layout = 'ordered'; cfg.channel = LFPTs.label;
layout = ft_prepare_layout(cfg,data_trl2)
%% Use Fieldtrip to construct spectograms
cfg = []; % Create empy cfg
cfg.output = 'pow';
cfg.channel = LFPTs.label;
cfg.method = 'mtmconvol'
cfg.taper = 'hanning';
cfg.foi = foi(1):foi(2):foi(3); % Frequencies of interest
cfg.keeptrials = 'yes'; % For stastical comparison
cfg.t_ftimwin = ones(size(cfg.foi)).*ftimwin;
cfg.toi = toi2(1):toi2(2):toi2(3); % Times of interest

TFR_event2 = ft_freqanalysis(cfg,data_trl2);

% Plot
figure
cfg = []; cfg.channel = fname; cfg.layout = layout;
ft_multiplotTFR(cfg,TFR_event2)

%% Statistical Comparison
TFR_event2.time = TFR_event1.time; % standardize time axis to binge
% T-Test
cfg = [];
cfg.channel     = LFPTs.label;
cfg.latency     = 'all';
cfg.trials      = 'all';
cfg.frequency   = 'all';
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
cfg.parameter   = 'powspctrm';
cfg.method      = 'stats';
cfg.statistic   = 'ttest2';
cfg.alpha       = 0.05;

nTrials1 = size(TFR_event2.powspctrm,1); nTrials2 = size(TFR_event1.powspctrm,1);
cfg.design = cat(2,ones(1,nTrials1),2*ones(1,nTrials2)); % two conditions
cfg.ivar = 1; % dimension of design var which contains the independent variable (group)
 
stat = ft_freqstatistics(cfg,TFR_event1,TFR_event2);

cfg.parameter = 'stat'; cfg.layout = layout;
figure
ft_multiplotTFR(cfg,stat); % plot the t-statistic