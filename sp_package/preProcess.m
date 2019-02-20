 function [LFPTs,chk_nan,zeroedChannel,clnTrls,trls,adfreq] = preProcess(LFPTs,adfreq,dsf,thresh,onset,offset,eoi,eventTs)
%% Applies preproccesing steps: filtering, downsampling, thresholding, and 
% trializing
% INPUTS:
% LFPTs = local field potential data structure; format: the following 
%   fields from Pl2tomvdmGenFile.m are necessary: 
%       .data = channel x time matrix
%       .tvec = vector of times of equal length as .data
% adfreq = sampling rate of data; format: whole number
% dsf = downsampling factor; format: integer
% thresh = y-axis threshold for removing artifacts; format: number
% onset = number of samples before threshold event to include in artifact;
%   format: integer
% offset = number of samples after threshold event to include in artifact;
%   format: integer
% eoi = events of interest; format = cell array of strings corresponding to
%   scored behaviors (e.g. 'binge', 'rest'; or special cases 'all',
%   'notbinge') and of numbers corresponding to window around that behavior
%   per trial (i.e. the minimum interval per behavioral trial)
% eventTs = event structure from ConvertPl2All_Files; fields .t and .label
%   required for script to run
%__________________________________________________________________________
% OUTPUTS:
% LFPTs = new LFPTs data; filtered, downsampled, and thresholded; 
% chk_nan = number of data points NaNed by threshFilt.m
% zeroedChannel = index of channels that have been completely zeroed due to
%   too many noise events
% clnTrls = cell array of timestamps of all clean trials; format = event x
%   trial > channel x timestamp
% clnEvents = collapsed clnTrls with all timestamps concatenated; format =
%   event > channel x timestamp
% trls = data structure of all usable trials, each event is in its own
%   structure; includes the following fields:
%       label = cell array of channel labels
%       fsample = sampling rate of data; samples per second
%       trial = cell array of clean data; format = channel x data
%       time = cell array of timestamps corresponding to clean data in 
%           .trial; format = row vector of timestamps
%       sampleinfo = start and stop times of each trial; format = trial x 
%           [start stop]
% adfreq = new sampling rate after downsampling
%__________________________________________________________________________
% USE:
% [LFPTs,chk_nan,zeroedChannel,clnTrls,clnEvents,trls,adfreq] = 
% preProcess(LFPTs,2000,5,2.5,5,17000,5,{'binge',[0 5];'rest',[0 5]},...
% eventTs)
% Will filter line noise out using filter60.m. Then downsample data and
% tvec by a factor of 5. NaN any datapoints with absolute values greater
% than 2.5 as well as 5 samples before and 17000 samples after. Lastly,
% data is trialized into two groups (binge and rest) using timestamps found
% in eventTswith each trial being 5 seconds long.
%__________________________________________________________________________
% DEPENDENCIES:
% filter60.m
% dwnSample.m
% threshFilt.m
% trialize.m
%   eventInd.m
%__________________________________________________________________________
% LLD 2016-17
%% Filter out 60 cycle line noise
disp('Applying 60 Hz filter with filter60.m...')
[LFPTs.data] = filter60(LFPTs,adfreq,0);
%% Downsample signal
disp('Downsampling signal with dwnSample.m...')
[LFPTs,adfreq] = dwnSample(LFPTs,dsf,adfreq);
%% Threshold to remove noise artifacts
disp('Removing noise artifacts with threshFilt.m...')
% Uses largest eoi minimum interval as minInt
minInt = max(diff(cell2mat(eoi(:,2)),1,2));
[LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,thresh,onset,offset,...
    minInt,adfreq);
%% NaN Sleep Intervals
% if exist('sleep')
%     for sind = 1:size(sleep.t{1,1},1)
%         LFPTs.data(:,nearest_idx2(sleep.t{1,1}(sind),LFPTs.tvec):nearest_idx2(sleep.t{1,2}(sind),LFPTs.tvec)) = NaN;
%     end
% end
%% Trialize data into equal length 'trials'
disp('Trializing data with trialize.m...')
[clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq);
