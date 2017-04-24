%% First run Matt and Wilder's code PL2tomvdm.m to convert pl2 data into
% usable data structures
[fname,pl2,n,freqs,j,lfpchan,ad,adfreq,fn,i,LFPTs,temp,TimeSampEr,ts,WBchan,eventTs]=Pl2tomvdm('I2FoodDep24Dec16_pl2done_plx.pl2');
%% Trialize Data Around Binge Events
% Create Trials
cfg = [];
cfg.t = eventTs.t{1,9}; % Timestamps corresponding to start of a binge
cfg.twin = [-8 2];% 10 second window around binge start with temporal 
% emphasis on before initiation since we are interested in the dynmaics 
% leading up to binging
% N.B.: due to the way the data is scored this time window will invariably 
% overlap with the 'approach' behvior; although a similar analysis on this 
% behavior specifically could also provide interesting information

clear trl; % Reset trl
trl(:,3) = nearest_idx3(cfg.t,LFPTs.tvec);
trl(:,2) = nearest_idx3(cfg.t+cfg.twin(2),LFPTs.tvec);
trl(:,1) = nearest_idx3(cfg.t+cfg.twin(1),LFPTs.tvec);

trl(:,3) = trl(:,1) - trl(:,3);
% Extract Binge Data
data_trlb = []; %set up empty trialized data structure
data_trlb.label = LFPTs.label;
data_trlb.fsample = adfreq;
data_trlb.trial = {}; % Set up empty trial structure
% Fill data_trl.trial structure with trialized data series
i = 1;
for i=1:size(trl,1)
    %data_trlb.trial{1,i} = LFPTs.data(1,trl(i,1):trl(i,2)) % Just first
    %channel
    data_trlb.trial{1,i} = LFPTs.data(1:4,trl(i,1):trl(i,2)) % All four channels
    end
% Create and fill data_trl.time structure
data_trlb.time = {};
twin = [cfg.twin(1):0.0005:cfg.twin(2)];

i = 1;
for i=1:size(trl,1)
    data_trlb.time{1,i} = twin; 
end
% Create and fill data_trl.sampleinfo
data_trlb.sampleinfo = [];
i = 1;
for i=1:size(trl,1)
    data_trlb.sampleinfo(i,1) = trl(i,1);
    data_trlb.sampleinfo(i,2) = trl(i,2);
end

% Define layout for later plotting
cfg = [];
cfg.layout = 'ordered'; cfg.channel = LFPTs.label; cfg.showlabels = 'yes';
layout = ft_prepare_layout(cfg,data_trlb)
%% Use Fieldtrip to construct spectograms
cfg = []; % Create empy cfg
cfg.output = 'pow';
cfg.channel = LFPTs.label;
cfg.method = 'mtmconvol'
cfg.taper = 'hanning';
cfg.foi = 1:2:150; % Frequencies of interest 1 through 150 in 2 Hz steps
cfg.keeptrials = 'yes'; % For stastical comparison
cfg.t_ftimwin = ones(size(cfg.foi)).*0.5; %0.5 second time windows
cfg.toi = -8:0.05:2; % Times of interest

TFR_binge = ft_freqanalysis(cfg,data_trlb);

% Plot
figure
cfg = []; cfg.channel = LFPTs.label; cfg.layout = layout;
ft_multiplotTFR(cfg,TFR_binge)

%% Trialize Data Around Rest Events
% For statistical purposes, next do the same type of analysis as above, but
% for periods of rest.

% Determine midpoints for each rest epoch
restS = eventTs.t{1,5};
restE = eventTs.t{1,6};
restM = restS(:) + (restE(:) - restS(:))/2;
cfg = [];
cfg.t = restM; % Timestamps corresponding to middle of rest
cfg.twin = [-5 5];% 10 second window within rest periods

% Trialize
clear trl; % Reset trl
trl(:,3) = nearest_idx3(cfg.t,LFPTs.tvec);
trl(:,2) = nearest_idx3(cfg.t+cfg.twin(2),LFPTs.tvec);
trl(:,1) = nearest_idx3(cfg.t+cfg.twin(1),LFPTs.tvec);

trl(:,3) = trl(:,1) - trl(:,3);
%% Extract Rest Data
data_trlr = []; %set up empty trialized data structure
data_trlr.label = LFPTs.label;
data_trlr.fsample = adfreq;
data_trlr.trial = {}; % Set up empty trial structure
% Fill data_trl.trial structure with trialized data series
i = 1;
for i=1:size(trl,1)
    % data_trlr.trial{1,i} = LFPTs.data(1,trl(i,1):trl(i,2)) % Just looking at first channel
    data_trlr.trial{1,i} = LFPTs.data(1:4,trl(i,1):trl(i,2)) % All four channels
end
% Create and fill data_trl.time structure
data_trlr.time = {};
twin = [cfg.twin(1):0.0005:cfg.twin(2)];
i = 1;
for i=1:size(trl,1)
    data_trlr.time{1,i} = twin;
end
% Create and fill data_trl.sampleinfo
data_trlr.sampleinfo = [];
i = 1;
for i=1:size(trl,1)
    data_trlr.sampleinfo(i,1) = trl(i,1);
    data_trlr.sampleinfo(i,2) = trl(i,2);
end
% Define layout for later plotting
cfg = [];
cfg.layout = 'ordered'; cfg.channel = LFPTs.label;
layout = ft_prepare_layout(cfg,data_trlr)
%% Use Fieldtrip to construct spectograms
cfg = []; % Create empy cfg
cfg.output = 'pow';
cfg.channel = LFPTs.label;
cfg.method = 'mtmconvol'
cfg.taper = 'hanning';
cfg.foi = 1:2:150; % Frequencies of interest 1 through 150 in 2 Hz steps
cfg.keeptrials = 'yes'; % For stastical comparison
cfg.t_ftimwin = ones(size(cfg.foi)).*0.5; %0.5 second time windows
cfg.toi = -5:0.05:5; % Times of interest

TFR_rest = ft_freqanalysis(cfg,data_trlr);

% Plot
figure
cfg = []; cfg.channel = fname; cfg.layout = layout;
ft_multiplotTFR(cfg,TFR_rest)

%% Statistical Comparison
TFR_rest.time = TFR_binge.time; % standardize time axis to binge
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

nTrials1 = size(TFR_rest.powspctrm,1); nTrials2 = size(TFR_binge.powspctrm,1);
cfg.design = cat(2,ones(1,nTrials1),2*ones(1,nTrials2)); % two conditions
cfg.ivar = 1; % dimension of design var which contains the independent variable (group)
 
stat = ft_freqstatistics(cfg,TFR_binge,TFR_rest);

cfg.parameter = 'stat'; cfg.layout = layout;
figure
ft_multiplotTFR(cfg,stat); % plot the t-statistic
%% Test Cell
% data_2.trial = {};
% i = 1;
% for i=1:size(trl,1)
%     data_2.trial{1,i} = LFPTs.data(1,trl(i,1):trl(i,2)) % Just looking at one channel
% end
% data_2.time = {};
% twin = [cfg.twin(1):0.0005:cfg.twin(2)];
% i = 1;
% for i=1:size(trl,1)
%     data_2.time{1,i} = twin;
% end
% data_2.sampleinfo = [];
% i = 1;
% for i=1:size(trl,1)
%     data_2.sampleinfo(i,1) = trl(i,1);
%     data_2.sampleinfo(i,2) = trl(i,2);
% end
