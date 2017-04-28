function [clnTrls,clnEvents,trls] = trialize(eoi,eventTs,LFPTs,adfreq)
%% Uses behavior markers to trialize data
%__________________________________________________________________________
% INPUTS:
% eoi = events of interest; format = cell array of strings corresponding to
%   scored behaviors (e.g. 'binge', 'rest'; or special cases 'all',
%   'notbinge') and of numbers corresponding to window around that behavior
%   per trial (i.e. the minimum interval per behavioral trial)
% eventTs = event structure from ConvertPl2All_Files; fields .t and .label
%   required for script to run
% LFPTs = data structure from ConvertPl2All_Files; fields .data and .tvec
%   required for script to run
% adfreq = sampling rate of data; format: integer
%__________________________________________________________________________
% OUTPUTS:
% clnTrls = cell array of timestamps of all clean trials; format = event x
%   trial > channel x timestamp
% clnEvents = collapsed clnTrls with all timestamps concatenated; format =
%   event > channel x timestamp
% trls = data structure of all usable trials, each event is in its own
%   structure; includes the following fields:
%       label = cell array of channel labels
%       fsample = sampling rate of data; samples per second
%       trial = 3D matrix of clean data; format = channel x data x trial
%       time = cell array of timestamps corresponding to clean data in 
%           .trial; format = row vector of timestamps
%       sampleinfo = start and stop times of each trial; format = trial x 
%           [start stop]
%__________________________________________________________________________
% USE: 
% [clnTrls,clnEvents,trls] = trialize({'binge',[0 3];'orient',[-1.5
% 1.5]},eventTs,LFPTs,adfreq)
% Will trialize behavior 'binge' with interval timestamps (start and stop)
% and 'orient' with scalar timestamp (each instance) using 3 second trials
% that either start at the beginning of binges or are center around
% orientations. 
%__________________________________________________________________________
% DEPENDENCIES:
% eventInd.m
% logicFind.m
% nearest_idx3.m
% nearest_idx3.mexw64
%__________________________________________________________________________
% LLD 2016-17
%% Initialize
% Count number of behaviors to trialize
nBehavior = size(eoi,1);
% Clean up eventTs, i.e. remove empty labels of empty timestamps
eventTs.label(cellfun(@isempty,eventTs.t)) = [];
eventTs.t(cellfun(@isempty,eventTs.t)) = [];
%% Prepares and runs eventInd.m to find indices of events to use.
% If behavior of interest is entire file, 'all', create new events/times
if strcmpi(eoi(:,1),'all')
    % Add new labels corresponding to start and end of data
    eventTs.label(end+1:end+2) = {'Start','End'};
    % Add new timestamps corresponding to start and end of data
    eventTs.t(end+1:end+2) = {LFPTs.tvec(1),LFPTs.tvec(end)};
    % Create eventInds 
    eventInds = [size(eventTs.label,2)-1,size(eventTs.label,2)];
    % Otherwise run eventInd.m to find correct indices for events
else%if sum(strcmpi(eoi(:,1),'app')) || sum(strcmpi(eoi(:,1),'binge')) || sum(strcmpi(eoi(:,1),'rest'))
    % If looking at 'notbinge' then remove it for eventInd to run
    if size(eoi,1) > 1 && strcmpi(eoi(2,1),'notbinge')
        [eventInds] = eventInd(eventTs,eoi(1,:));
    else
        [eventInds] = eventInd(eventTs,eoi);
    end
end
%% Get timestamps for all trials of each behavior
% Preallocate intTime array using the size of the data if looking
% at 'all' data
if strcmpi(eoi(:,1),'all')
    intTime{1,1} = 1:length(LFPTs.data);
    % Otherwise, preallocate using number of behaviors and largest number
    % of trials per behavior.
else
    nTrial = zeros(1,nBehavior);
    for iB = 1:nBehavior
        if ~strcmpi(eoi(iB,1),'notbinge')
            nTrial(iB) = numel(eventTs.t{eventInds(iB,1)});
        end
    end
    largest = max(nTrial);
    intTime = cell(nBehavior,largest);
    % Fill intTime array with .tvec indices
    for iB = 1:nBehavior
        if strcmpi(eoi(iB,1),'notbinge')
            intTime{2,1} = (1:(intTime{1,1}(1,1)-1));
            for iT = 2:size(intTime(1,:),2)+1
                if iT == size(intTime(1,:),2)+1
                    intTime{2,iT} = intTime{1,iT-1}(1,end)+1:size(LFPTs.tvec,2);
                else
                    intTime{2,iT} = intTime{1,iT-1}(1,end)+1:intTime{1,iT}(1,1)-1;
                end
            end
        else
            for iT = 1:numel(eventTs.t{eventInds(iB,1)})
                % Use eventInds corresponding to start and stop if an interval
                % behavior
                if ~isequal(eventInds(iB,1),eventInds(iB,2))
                    intTime{iB,iT} = nearest_idx3(eventTs.t{1,eventInds(iB,1)}(iT),LFPTs.tvec):nearest_idx3(eventTs.t{1,(eventInds(iB,2))}(iT),LFPTs.tvec);
                    % Otherwise it is a scalar behavior, so look around those
                    % timestamps in a window defined by the corresponding entry in
                    % the second column of 'eoi'
                else
                    intTime{iB,iT} = nearest_idx3(eventTs.t{1,eventInds(iB,1)}(iT)+eoi{iB,2}(1),LFPTs.tvec):nearest_idx3(eventTs.t{1,(eventInds(iB,2))}(iT)+eoi{iB,2}(2),LFPTs.tvec);
                end
            end
        end
    end
end
%% If looking at 'notbinge' then uses time stamps for the end of one binge
% and the beginning of the next as start and stops.
% N.B.: ASSUMES THAT THE RECORDING DOES NOT START OR STOP DURING A BINGE
% if nBehavior >1 && strcmpi(eoi(2,1),'notbinge')
%    intTime{2,1} = repmat((1:(intTime{1,1}(1,1)-1)),4,1);
%    for iT = 2:size(intTime(1,:),2)+1
%       if iT == size(intTime(1,:),2)+1
%           intTime{2,iT} = intTime{1,iT-1}(1,end)+1:size(LFPTs.tvec,2);
%       else
%           intTime{2,iT} = intTime{1,iT-1}(1,end)+1:intTime{1,iT}(1,1)-1;
%       end
%    end
% end
%% NaN timestamps corresponding to NaNed LFPTs data
for iB = 1:nBehavior
    for iT = 1:sum(~cellfun(@isempty,intTime(iB,:)))
        intTime{iB,iT}(isnan(LFPTs.data(iB,intTime{iB,iT}))) = NaN;
    end
end
%% NaN contiguous data intervals less than minInt
clnTrls = cell(size(intTime,1),sum(~cellfun(@isempty,intTime(iB,:))));
for iB = 1:nBehavior
    % Get this behaviors minimum interval
    minInt = diff(cell2mat(eoi(iB,2)),1,2);
    for iT = 1:sum(~cellfun(@isempty,intTime(iB,:)))
        dummy = intTime{iB,iT}; 
        dummy(~isnan(dummy)) = 1; 
        dummy(isnan(dummy)) = 0;
        dataStart = find(diff(dummy)==1)+1;
        dataStop = find(diff(dummy)==-1);
        % Check for clean data at first index
        if ~isnan(intTime{iB,iT}(1,1))
            dataStart = [1, dataStart]; %#ok<AGROW>
        end
        % Check for clean data at last index
        if ~isnan(intTime{iB,iT}(end))
            dataStop = horzcat(dataStop,length(intTime{iB,iT})); %#ok<AGROW>
        end
        % Only keep data that is not NaNed and in continuous intervals
        % longer than minInt (+ 1 for indexing)
        if ~isnan(sum(intTime{iB,iT})) && (size(intTime{iB,iT},2) >= (minInt*adfreq + 1/adfreq)) 
            numTrls = floor(size(intTime{iB,iT},2)/(minInt*adfreq + 1/adfreq));
            thisTrls = intTime{iB,iT}(1:((minInt*adfreq)*numTrls));
            clnTrls{iB,iT} = thisTrls;
        else
            thisTrls = [];
            % Run through data intervals
            for intInd = 1:length(dataStart) 
                intLen = dataStop(intInd) - dataStart(intInd);
                % Double check that each interval is long enough
                if intLen >= (minInt*adfreq + 1/adfreq) 
                    numTrls = floor(intLen/(minInt*adfreq + 1/adfreq));
                    thisTrls = horzcat(thisTrls,intTime{iB,iT}(dataStart(intInd):(dataStart(intInd)+((minInt*adfreq)*numTrls))-1)); %#ok<AGROW>
                end
            end
            if ~isempty(thisTrls)
                clnTrls{iB,iT} = thisTrls;
            end
        end
    end
end
%% Collapse data across behavior
clnEvents = cell(size(clnTrls,1),1);
for iB = 1:nBehavior
    for iT = 1:size(clnTrls,2)
        clnEvents{iB} = horzcat(clnEvents{iB},clnTrls{iB,iT});
    end
end
%% Create empty trial structure
thisTrl = [];
thisTrl.label = LFPTs.label;
thisTrl.fsample = adfreq;
thisTrl.trial = [];
thisTrl.time = {};
thisTrl.sampleinfo = [];
% Create trl structures for each event
trls = cell(1,size(eoi,1));
for iB = 1:size(eoi,1)
    if ~isempty(clnEvents{iB,1})
        nTrls = length(clnEvents{iB,1})/(minInt*adfreq);
        for iT = 1:nTrls
            % To account for the transition of zero indexed time, start
            % 1 sample after 0
            thisTrl.time{1,iT} = (1/adfreq:1/adfreq:minInt);
            thisTrl.sampleinfo(iT,1) = clnEvents{iB,1}(1,1+minInt*adfreq*(iT-1));
            thisTrl.sampleinfo(iT,2) = clnEvents{iB,1}(1,minInt*adfreq*iT);
            thisTrl.trial(:,:,iT) = LFPTs.data(:,clnEvents{iB,1}(1,(1+minInt*adfreq*(iT-1)):minInt*adfreq*iT));
        end
        trls{iB} = thisTrl;
    else
        disp('Warning: This event has 0 clean trials!')
    end
end