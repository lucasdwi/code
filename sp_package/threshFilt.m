function [LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq)
%% Find indices of |data| greater than or equal to thresh and NaNs. Useful 
% for filtering out noise artifacts from data.
%__________________________________________________________________________
% INPUTS:
% LFPTs = LFPTs data structure from Pl2tomvdmGenFile.m
% thresh = maximum y-axis |value| at which to NaN data; format: number
%   (e.g. 2.5 will NaN any data above/equal to 2.5 or below/equal to 2.5
% onset = amount of data before threshold cross to NaN; format: seconds
%   N.B. onset*adfreq will be rounded up if not a whole number
% offset = amount of data after threshold cross to NaN; format: seconds
%   N.B. offset*adfreq will be rounded up if not a whole number
% minInt = minimum interval length of data to keep; format: seconds 
% adfreq = sampling frequency; format: Hz
%__________________________________________________________________________
% OUTPUTS:
% LFPTs = NaNed LFPT data structure
% nNaN = number of NaNed samples per channel
% indSkp = indices of channels skipped due to too many NaNs
%__________________________________________________________________________
% USE:
% [LFPTs,chk_nan,zeroedChannel] = threshFilt(LFPTs,2.5,0.0125,40,5,400)
% Will find all instances where |LFPTs.data| > 2.5 and then NaN windows
% around each of those instances that extends 0.0125 seconds before the
% crossing and 40 seconds after the crossing (5 samples before and 16000
% samples after [e.g. 0.0125*400=5]). Lastly, all intervals of clean data
% that are less than 5 seconds long (1200 samples) will also be NaNed.
%__________________________________________________________________________
% DEPENDENCIES:
% logicFind.m
%__________________________________________________________________________
% LLD 2016-17
%% Create threshInd structure; 1 = above threshhold, 0 = not
% Count number of channels
chans = size(LFPTs.data,1);
% Set up threshInd to have one cell per channel to keep track of all
% threshold crossing indices
threshInd = cell(1,chans); % Columns = channels; {1,:} = indices above thresh; {2,:} = onset indices; {3,:} = offset indices
% Set up zeroedChannel to keep track of channels that were too noisy (>30%
% noise)
zeroedChannel = [];
% Set chk_nan to zero
chk_nan = 0;
for iI = 1:chans
    threshInd{iI} = logicFind(thresh,abs(LFPTs.data(iI,:)),'>');
    % Checks if amount of data beyond threhold exceeds 50%, if so set whole
    % channel to 0
    if length(threshInd{iI}) >= 0.5*length(LFPTs.data(iI,:))
       LFPTs.data(iI,:) = 0;
       % Keep trach of zeroed channels
       zeroedChannel = [zeroedChannel,iI]; %#ok<AGROW>
       % Reset threshInd to empty
       threshInd{iI} = [];
    end
    % Checks if none of the channels exceed threshold
    chk_nan = chk_nan + numel(threshInd{iI});
end
%%
if chk_nan ~= 0
    %% Find onset and offset cuttoff
    % Convert onse and offset to samples (seconds*samples/second = samples)
    % and round up if needed to avoid fractional indexing
    onset = ceil(onset*adfreq);
    offset = ceil(offset*adfreq);
    for iI=1:chans
        % Skip to the first usable index to avoid negative indexing
        firstIndex = logicFind(onset,threshInd{1,iI},'>');
        threshInd{2,iI} = threshInd{1,iI}(firstIndex)-onset;
        threshInd{3,iI} = threshInd{1,iI}(firstIndex)+offset;
    end

    %% Create overlapping intervals
    % Initialize overlapping interval array
    overInt = cell(2,chans);
    for iI = 1:chans
        if ~isempty(threshInd{1,iI})
            startI = threshInd{2,iI}(1);
            endI = threshInd{3,iI}(1);
            k = 0;
            for j=1:length(threshInd{2,iI})
                if threshInd{2,iI}(j)> endI
                    k = k+1;
                    overInt{1,iI}(k) = startI;
                    overInt{2,iI}(k) = endI;
                    startI = threshInd{2,iI}(j);
                    endI = threshInd{3,iI}(j);
                else
                    endI = max(endI,threshInd{3,iI}(j));
                end
            end
            %Stack last interval
            k = k+1;
            overInt{1,iI}(k) = startI;
            overInt{2,iI}(k) = endI;
            % Check if last interval goes beyond length of data, if so cap it
            if overInt{2,iI}(end) > size(LFPTs.data,2)
                overInt{2,iI}(end) = size(LFPTs.data,2);
            end
        end
    end
    % Idea borrowed from MergeBrackets.m Author: Bruno Luong
    % <brunoluong@yahoo.com> Original: 25-May-2009
    %% NaN overlapping intervals
    % Copy data to be NaNed
    LFPTsNaN = LFPTs;
    for iI = 1:chans
        if ~isempty(overInt{1,iI})
            % NaN all data across channels corresponding to indices in
            % overInt
            for j=1:length(overInt{1,iI})
                LFPTsNaN.data(iI,overInt{1,iI}(j):overInt{2,iI}(j)) = NaN;
            end
        end
    end
    %% Find chunks of existing data
    % Sum each column across rows to get all channels to have the same number
    % of NaNed indices per channel; N.B. a sum that includes >= 1 NaN will be
    % set to NaN
    oneChan = sum(LFPTsNaN.data,1);
    % Spread oneChan across all channels
    LFPTsNaN.data(:,isnan(oneChan)) = NaN; 
    % Set oneChan to 1 if good data and to 0 if NaN
    oneChan(~isnan(oneChan)) = 1; 
    oneChan(isnan(oneChan)) = 0;
    % Use diff() to find where good data starts and stops
    dataStart = logicFind(1,diff(oneChan),'==')+1;
    dataStop = logicFind(-1,diff(oneChan),'==');
    % Above misses first and last indices; check these for clean data and
    % concatenate with above vectors
    if oneChan(1) == 1
        dataStart = [1, dataStart];
    end
    if oneChan(end) == 1
        dataStop = horzcat(dataStop,length(oneChan));
    end
    %% Check size of intervals and if < minInt, then NaN
    for iI = 1:length(dataStart)
        if dataStop(iI) - dataStart(iI) <= (minInt*adfreq+1)
            LFPTsNaN.data(:,dataStart(iI):dataStop(iI)) = NaN;
        end
    end
    %% Overwrite LFPTs with LFPTsNaN
    LFPTs = LFPTsNaN;
    % Record that LFPTs was processed with threshFilt.m
%     LFPTs.cfg.history.mfun{end+1} = 'threshFilt';
end