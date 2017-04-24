function [trl1,trl2] = trialExtract(event1,event2,eventInds,eventTs,LFPTsNaN,adfreq,minInt)
%% Trialize LFPTsNaN.data based on event markers
% Requires output of threshFilt.m and eventInd.m

%% Event Marker Cheat Sheet
% 3 = Orientation
% 8 = Rest Start; 9 = Rest Stop; 13 = Rest Interval
% 4 = Approach Start; 5 = Approach Stop; 14 = Approach Interval
% 6 = Binge Start; 7 = Binge Stop; 13 = Binge Interval
%% Setup
chans = size(LFPTsNaN.label,2); %Count channels
discEvents = [3,8,9,4,5,6,7]; %Indices for discrete behavior markers except Rest Middle (11)

%% Creates cell array of continuous behavior epoch timestamps
% intTime{i,j}(k,l)
% i = behavior (approach, binge, rest)
% j = epoch number
% k = channel
% l = timestamp

markers = [4,6,8]; %[Approach, Binge, Rest]

% Preallocate intTime array
numMark = length(markers);
for i = 1:length(markers)
    numMark(i) = numel(eventTs.t{markers(i)}); 
end
largest = max(numMark);
intTime = cell(length(markers),largest);

% Fill intTime array with .tvec indices
for i = 1:length(markers)
    tic
    for j = 1:numel(eventTs.t{markers(i)})
        intTime{i,j} = nearest_idx(eventTs.t{1,markers(i)}(j),LFPTsNaN.tvec):nearest_idx(eventTs.t{1,(markers(i)+1)}(j),LFPTsNaN.tvec);
        for k = 1:chans
            intTime{i,j}(k,:) = intTime{i,j}(1,:); %Copy channel one indices to all other channels
        end
    end
    toc
end
%% NaN timestamps corresponding to NaNed LFPTs data
for i = 1:length(markers)
    for j = 1:numel(eventTs.t{markers(i)})
        for k = 1:chans
            intTime{i,j}(k,isnan(LFPTsNaN.data(i,intTime{i,j}(k,:)))) = NaN;
        end
    end
end
%% NaN contiguous data intervals less than minInt
tic
clnTrls = cell(length(markers),largest); %Preallocate clnTrls
for i = 1:length(markers)
    for j = 1:numel(eventTs.t{markers(i)})
        for k = 1:chans
            A = intTime{i,j}(k,:); A(~isnan(A)) = 1; A(isnan(A)) = 0;
            dataStart = find(diff(A)==1)+1;
            dataStop = find(diff(A)==-1);
            if ~isnan(intTime{i,j}(k,1))
                dataStart = [1, dataStart];
            end
            if ~isnan(intTime{i,j}(k,end))
                dataStop = horzcat(dataStop,length(intTime{i,j}));
            end
            
            if ~isnan(sum(intTime{i,j}(k,:))) && LFPTsNaN.tvec(intTime{i,j}(k,end))-LFPTsNaN.tvec(intTime{i,j}(k,1)) >= (minInt + 0.0005) %No NaNed data and longer than minInt + 1 for indexing
                numTrls = floor((LFPTsNaN.tvec(intTime{i,j}(k,end))-LFPTsNaN.tvec(intTime{i,j}(k,1)))/(minInt + 0.0005));
                thisTrls = intTime{i,j}(k,1:((minInt*2000+1)*numTrls));
                clnTrls{i,j}(k,:) = thisTrls;
            else for intInd = 1:length(dataStart) %Run through data intervals
                    intLen = LFPTsNaN.tvec(dataStop(intInd)) - LFPTsNaN.tvec(dataStart(intInd));
                    thisTrls = [];
                    if intLen >= (minInt + 0.0005) %Keep if big enough
                        numTrls = floor(intLen/(minInt + 0.0005));
                        thisTrls = horzcat(thisTrls,intTime{i,j}(k,dataStart(intInd):(dataStart(intInd)+((minInt*2000+1)*numTrls))-1));
                    end
                end
                if ~isempty(thisTrls)
                    clnTrls{i,j}(k,:) = thisTrls;
                end
            end
        end
    end
end
toc      
%% Collapse data across behavior
clnEvents = cell(length(markers),1);
for i = 1:length(markers)
    for j = 1:largest
        clnEvents{i} = horzcat(clnEvents{i},clnTrls{i,j});
    end
end
%% Create trial structures
trl1 = [];
trl1.label = LFPTsNaN.label;
trl1.fsample = adfreq;
trl1.trial = {};
trl1.time = {};
trl1.sampleinfo = [];
trl2 = trl1; %Copy trl1 structure for trl2

%% Fill trl1 structure with extracted clean data
ind = find(markers == event1);
nTrls = length(clnEvents{ind,1})/((minInt*2000)+1);
for i = 1:nTrls
    trl1.time{1,i} = (0:0.0005:minInt);
    trl1.sampleinfo(i,1) = clnEvents{ind,1}(1,1+minInt*2000*(i-1));
    trl1.sampleinfo(i,2) = clnEvents{ind,1}(1,minInt*2000*i);
    trl1.trial{1,i} = LFPTsNaN.data(1:4,clnEvents{ind,1}(1,(1+minInt*2000*(i-1)):1+minInt*2000*i));
end
%% Fill trl2 structure with extracted clean data
ind = find(markers == event2);
nTrls = length(clnEvents{ind,1})/((minInt*2000)+1);
for i = 1:nTrls
    trl2.time{1,i} = (0:0.0005:minInt);
    trl2.sampleinfo(i,1) = clnEvents{ind,1}(1,1+minInt*2000*(i-1));
    trl2.sampleinfo(i,2) = clnEvents{ind,1}(1,minInt*2000*i);
    trl2.trial{1,i} = LFPTsNaN.data(1:4,clnEvents{ind,1}(1,(1+minInt*2000*(i-1)):1+minInt*2000*i));
end
% %% Use Fieldtrip for spectrogram
% cfg = []; % Create empy cfg
% cfg.output = 'pow';
% cfg.channel = LFPTs.label;
% cfg.method = 'mtmconvol';
% cfg.taper = 'hanning';
% cfg.foi = foi(1):foi(2):foi(3); % Frequencies of interest
% cfg.keeptrials = 'yes'; % For stastical comparison
% cfg.t_ftimwin = ones(size(cfg.foi)).*ftimwin;
% cfg.toi = toi(1):toi(2):toi(3); % Times of interest
% TFR_event1 = ft_freqanalysis(cfg,trl1);
% %% Use Fieldtrip for spectrogram
% cfg = []; % Create empy cfg
% cfg.output = 'pow';
% cfg.channel = LFPTs.label;
% cfg.method = 'mtmconvol';
% cfg.taper = 'hanning';
% cfg.foi = foi(1):foi(2):foi(3); % Frequencies of interest
% cfg.keeptrials = 'yes'; % For stastical comparison
% cfg.t_ftimwin = ones(size(cfg.foi)).*ftimwin;
% cfg.toi = toi(1):toi(2):toi(3); % Times of interest
% TFR_event2 = ft_freqanalysis(cfg,trl2);
