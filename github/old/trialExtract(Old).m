function [data_trl1,data_trl2,baseTimeVect] = trialExtract(event1,event2,eventInds,twin1,twin2,eventTs,LFPTs,adfreq)
%% Trialize data around given event and extract that data
% event = event of interest to trialize around

%% Create trials
% Setup cfg for trialization
cfg = [];
cfg.t1 = eventTs.t{1,eventInds(event1)}; % Timestamps of event1
cfg.t2 = eventTs.t{1,eventInds(event2)}; % Timestamps of event2
cfg.twin1 = twin1; % time window around event 1
cfg.twin2 = twin2; % time window around event 2

% Creat event1 trials
trl1(:,3) = nearest_idx3(cfg.t1,LFPTs.tvec);
trl1(:,2) = nearest_idx3(cfg.t1+cfg.twin1(2),LFPTs.tvec);
trl1(:,1) = nearest_idx3(cfg.t1+cfg.twin1(1),LFPTs.tvec);

trl1(:,3) = trl1(:,1) - trl1(:,3);

% Create event2 trials

trl2(:,3) = nearest_idx3(cfg.t2,LFPTs.tvec);
trl2(:,2) = nearest_idx3(cfg.t2+cfg.twin2(2),LFPTs.tvec);
trl2(:,1) = nearest_idx3(cfg.t2+cfg.twin2(1),LFPTs.tvec);

trl2(:,3) = trl2(:,1) - trl2(:,3);
%% Set up data_trl structures for both events
data_trl1 = []; data_trl2 = []; %Set up empty trialized data structures
data_trl1.label = LFPTs.label; %Transfer LFPTs labels from data
data_trl1.fsample = adfreq; %Transfer adfreq from data
data_trl1.trial = {}; %Setup empty trial structure
data_trl1.time = {}; %Set up empty time structure
data_trl1.sampleinfo = []; %Setup empty sampleinfo
data_trl2 = data_trl1; %Copy data_trl1 as data_trl2, they will differentiate later

%% Fill data_trl structure with extracted trialized data series, time series, and sampleinfo 
trls = {trl1,trl2};
i = 1;
j = 1;
for i=1:size(trls,2)
    for j=1:size(trls{i},1)
        if i == 1
            data_trl1.trial{1,j} = LFPTs.data(1:4,trl1(j,1):trl1(j,2));
            data_trl1.time{1,j} = [cfg.twin1(1):0.0005:cfg.twin1(2)];
            data_trl1.sampleinfo(j,1) = trl1(j,1);
            data_trl1.sampleinfo(j,2) = trl1(j,2);
        elseif i == 2 
            data_trl2.trial{1,j} = LFPTs.data(1:4,trl2(j,1):trl2(j,2));
            data_trl2.time{1,j} = [cfg.twin2(1):0.0005:cfg.twin2(2)];
            data_trl2.sampleinfo(j,1) = trl2(j,1);
            data_trl2.sampleinfo(j,2) = trl2(j,2);
        end
    end
end
%% Create vector of timestamps for baseline
baseTimeVect = [];
for i=1:size(data_trl2.sampleinfo,1)
    thisRow = data_trl2.sampleinfo(i,1):data_trl2.sampleinfo(i,2);
    baseTimeVect = horzcat(baseTimeVect,thisRow);
end
