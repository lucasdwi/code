function [clnTrls,trls] = trialize(eoi,eventTs,LFPTs,adfreq,fixed,discrete)
% [clnTrls,clnEvents,trls] = trialize(eoi,eventTs,LFPTs,adfreq)
%% Uses behavior markers to trialize data
%__________________________________________________________________________
% INPUTS:
% eoi = events of interest; format = cell array of strings corresponding to
%   scored behaviors (e.g. 'binge', 'rest'; or special cases 'all',
%   'notbinge') and of numbers corresponding to window around that behavior
%   per trial; if a single positive number, then defines the minimum trial
%   length; if a set of numbers, then defines the 'window' around the given
%   marker
% eventTs = event structure from ConvertPl2All_Files; fields .t and .label
%   required for script to run
% LFPTs = data structure from ConvertPl2All_Files; fields .data and .tvec
%   required for script to run
% adfreq = sampling rate of data; format: integer 
% fixed = whether or not to use fixed trializing; if 1 (fixed), then will
%   use contiguous trials that disregard NaNs; if 0 (not fixed) then will
%   maximize the number of trials without NaNs
% discrete = whether or not analysis is discrete; 0 for continuous, 1 for
%   discrete
%__________________________________________________________________________
% OUTPUTS:
% clnTrls = cell array of timestamps of all clean trials; format = event X
%   trial > channel X timestamp
% clnEvents = collapsed clnTrls with all timestamps concatenated; format =
%   event > channel X timestamp
% trls = data structure of all usable trials, each event is in its own
%   structure; includes the following fields:
%       label = cell array of channel labels
%       fsample = sampling rate of data; samples per second
%       trial = 3D matrix of clean data; format = channel X data X trial
%       time = cell array of timestamps corresponding to clean data in
%           .trial; format = row vector of timestamps
%       sampleinfo = start and stop times of each trial; format = trial X
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
% Count events
nEvent = size(eoi,1);
%% Check for negated, '~', behaviors. If they exist, remove negation for
% eventInd.m
neg = false(1,nEvent);
for ii = 1:nEvent
    neg(ii) = contains(eoi{ii,1},'~');
end
% Remove ~ from eoi so eventInd.m runs; if the positive version of event
% also exists, there will be a repeat which is fine
if any(neg)
    eoi{neg,1} = eoi{neg,1}(2:end);
end
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
else
    [eventInds] = eventInd(eventTs,eoi);
end
% Re-negate events to avoid confusion
if any(neg)
    eoi{neg,1} = ['~',eoi{neg,1}];
end
%% Get timestamps for all trials of each behavior
% Preallocate intTime array using the size of the data if looking
% at 'all' data
if strcmpi(eoi(:,1),'all')
    intTime{1,1} = 1:length(LFPTs.data);
end
% Fill intTime array with .tvec indices
for iE = 1:nEvent
    % Check that given event has any timestamps
    if ~isempty(eventTs.t{eventInds(iE,1)})
        % Check if event is negated
        if neg(iE)
            % Check that first positive event does not occur at exact start of
            % file
            if eventTs.t{eventInds(iE,1)}(1)~=LFPTs.tvec(1)
                % If not, use time 1 as beginning of first interval and
                % first Start time as end of interval
                intTime{iE,1} = 1:nearest_idx3(eventTs.t{eventInds(iE,1)}(1)...
                    ,LFPTs.tvec);
                % Set starting index at 2 to skip first event
                c = 2;
            end
        end
        % Run through all events
        for iT = 1:numel(eventTs.t{eventInds(iE,1)})
            % Use eventInds corresponding to start and stop if an interval
            % behavior
            if ~isequal(eventInds(iE,1),eventInds(iE,2))
                % Check if event is negated
                if neg(iE)
                    % Check if last interval
                    if iT == numel(eventTs.t{eventInds(iE,1)})
                        % If so first check if last interval goes to end of
                        % file in which case skip, otherwise make interval that
                        % goes from end of last interval to end of recording
                        
                        % Use counter, c, as storage index to account for
                        % potential first interval created outside of loop
                        if eventTs.t{eventInds(iE,2)}(end) ~= LFPTs.tvec(end)
                            intTime{iE,c} = nearest_idx3(...
                                eventTs.t{1,eventInds(iE,2)}(iT),LFPTs.tvec)...
                                :numel(LFPTs.tvec);
                        end
                    else
                        % If not last then use the start of an event as the end
                        % of the previous interval, and the end as the start of
                        % the current
                        intTime{iE,c} = nearest_idx3(...
                            eventTs.t{1,eventInds(iE,2)}(iT),LFPTs.tvec):...
                            nearest_idx3(eventTs.t{1,eventInds(iE,1)}(iT+1),...
                            LFPTs.tvec);
                    end
                    % Add one to counter
                    c = c+1;
                else
                    intTime{iE,iT} = nearest_idx3(eventTs.t{1,...
                        eventInds(iE,1)}(iT),LFPTs.tvec):nearest_idx3(...
                        eventTs.t{1,eventInds(iE,2)}(iT),LFPTs.tvec);
                end
                % Otherwise it is a scalar behavior, so look around those
                % timestamps in a window defined by the corresponding entry in the
                % second column of 'eoi'
            else
                intTime{iE,iT} = nearest_idx3(...
                    eventTs.t{1,eventInds(iE,1)}(iT)+eoi{iE,2}(1),...
                    LFPTs.tvec):nearest_idx3(...
                    eventTs.t{1,(eventInds(iE,2))}(iT)+eoi{iE,2}(2),...
                    LFPTs.tvec);
            end
        end
    else
        if neg(iE)
            intTime{iE,1} = 1:length(LFPTs.data);
        else
            disp(['Warning: ',eoi{iE,1},' has no associated timestamps'])
        end
    end
end
%% NaN timestamps corresponding to NaNed LFPTs data
% Only do this if doing discrete analysis; for continuous analysis the NaNs
% will be taken care of during processing
if discrete
    for iE = 1:nEvent
        % Skip if event has no timestamps
        if ~isempty(intTime{iE})
            for iT = 1:sum(~cellfun(@isempty,intTime(iE,:)))
                % Sums across channels of data to propogate missing data
                intTime{iE,iT}(isnan(sum(LFPTs.data(:,intTime{iE,iT}),1))) =...
                    NaN;
            end
        end
    end
end
%% NaN contiguous data intervals less than the minimum trial length, which
% in the case of two time values in the eoi will be calculated as the time
% between them
clnTrls = cell(size(intTime));
for iE = logicFind(0,cellfun(@isempty,intTime(:,1)),'==')
    % Get this behaviors minimum interval
    if size(eoi{iE,2},2) == 2
        minInt = diff(cell2mat(eoi(iE,2)),1,2);
    elseif size(eoi{iE,2},2) == 1
        minInt = eoi{iE,2};
    end
    for iT = 1:sum(~cellfun(@isempty,intTime(iE,:)))
        % If not doing discrete analysis, then put all intTime into clnTrls
        if ~discrete
            clnTrls{iE,iT} = intTime{iE,iT};
        % If using fixed windows, then use maximal amount of data that
        % will evenly split into the most windows
        elseif fixed
            numTrls = floor(size(intTime{iE,iT},2)/(minInt*adfreq));
            clnTrls{iE,iT} = intTime{iE,iT}(1:numTrls*(minInt*adfreq));
        % Otherwise, figure out where clean data windows exist that are
        % long enough
        else
            dummy = intTime{iE,iT};
            dummy(~isnan(dummy)) = 1;
            dummy(isnan(dummy)) = 0;
            dataStart = find(diff(dummy)==1)+1;
            dataStop = find(diff(dummy)==-1);
            % Check for clean data at first index
            if ~isnan(intTime{iE,iT}(1,1))
                dataStart = [1, dataStart]; %#ok<AGROW>
            end
            % Check for clean data at last index
            if ~isnan(intTime{iE,iT}(end))
                dataStop = horzcat(dataStop,length(intTime{iE,iT}));%#ok<AGROW>
            end
            % Only keep data that is not NaNed and in continuous intervals
            % longer than minInt
            if ~isnan(sum(intTime{iE,iT})) && (size(intTime{iE,iT},2) >= ...
                    (minInt*adfreq))
                numTrls = floor(size(intTime{iE,iT},2)/(minInt*adfreq));
                thisTrls = intTime{iE,iT}(1:((minInt*adfreq)*numTrls));
                clnTrls{iE,iT} = thisTrls;
            else
                thisTrls = [];
                % Run through data intervals
                for intInd = 1:length(dataStart)
                    % Add one for indexing
                    intLen = length(dataStart(intInd):dataStop(intInd));
                    % Double check that each interval is long enough
                    if intLen >= (minInt*adfreq)
                        numTrls = floor(intLen/(minInt*adfreq));
                        thisTrls = horzcat(thisTrls,intTime{iE,iT}(...
                            dataStart(intInd):(dataStart(intInd)+...
                            ((minInt*adfreq)*numTrls))-1)); %#ok<AGROW>
                    end
                end
                if ~isempty(thisTrls)
                    clnTrls{iE,iT} = thisTrls;
                end
            end
        end
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
for iE = 1:size(eoi,1)
    % Check if behavior has any data
    if any(~cell2mat(cellfun(@isempty,clnTrls(iE,:),'UniformOutput',0)))
        % Start counter
        c = 1;
        for iT = 1:size(clnTrls,2)
            thisTrl = [];
            thisTrl.label = LFPTs.label;
            thisTrl.fsample = adfreq;
            thisTrl.trial = [];
            thisTrl.time = {};
            thisTrl.sampleinfo = [];
            % If doing discrete analysis, figure out maximal windowing
            if discrete
                % Check if trial exists
                if ~isempty(clnTrls{iE,iT})
                    % Split into as many trials as possible
                    for iN = 1:size(clnTrls{iE,iT},2)/(minInt*adfreq)
                        thisTrl.time{1,c} = (1/adfreq:1/adfreq:minInt);
                        thisTrl.sampleinfo(c,1) = clnTrls{iE,iT}...
                            (1+(iN-1)*minInt*adfreq);
                        thisTrl.sampleinfo(c,2) = clnTrls{iE,iT}...
                            (iN*minInt*adfreq);
                        % If using fixed trials, then check if any NaNs exist
                        % in which case make a whole trial of NaNs as
                        % placeholder
                        if fixed && sum(isnan(clnTrls{iE,iT}(1+(iN-1)*minInt...
                                *adfreq:iN*minInt*adfreq)))
                            thisTrl.trial(:,:,c) = NaN([size(LFPTs.data,1),...
                                minInt*adfreq]);
                        else
                            thisTrl.trial(:,:,c) = LFPTs.data(:,...
                                clnTrls{iE,iT}(1+(iN-1)*minInt*adfreq:...
                                iN*minInt*adfreq));
                        end
                        c = c+1;
                    end
                end
            % Otherwise, if doing continuous analysis, keep all data
            else
                thisTrl.trial = LFPTs.data(:,clnTrls{iE,iT});
                thisTrl.sampleinfo = clnTrls{iE,iT};
            end
        end
        trls{iE} = thisTrl;
        thisTrl = [];
    else
        disp('Warning: This event has 0 clean trials!')
    end
end