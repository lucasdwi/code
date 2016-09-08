function [LFPTs] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq,dsf,chans)
%function [LFPTs,nNaN,indSkp] = threshFilt(LFPTs,thresh,onset,offset,minInt,NaNcutoff,adfreq,dsf,chans)
%% Find indices of values > |thres| and NaN

% Inputs:
% LFPTs = LFPTs data structure from ConvertPl2All_Files.m (Stott)
% thresh = maximum value at which to NaN data; format: mV (e.g. 2.5 will
%   NaN any data above/equal to 2.5 or below/equal to 2.5
% onset = amount of data before threshold cross to NaN; format: samples 
%   (e.g. 1 = 0.5 ms; 2000 = 1 sec) 
% offset = amount of data after threshold cross to NaN; format: samples
%   (e.g. 1 = 0.5 ms; 2000 = 1 sec)  
% minInt = minimum interval length of data to keep; format: seconds 
% NaNcutoff = number of standard deviations from the average number of NaNs
%   per channel at which to get rid of channel; format = positive number
%   (e.g. 1.5)
% adfreq = sampling frequency; format = Hz
% dsf = downsampling factor; format = integer
% chans = number of channels; format = integer

% Outputs:
% LFPTs = NaNed LFPT data structure
% nNaN = number of NaNed samples per channel
% indSkp = indices of channels skipped due to too many NaNs

% Code:
% chkNaN.m

%% Create threshInd structure; 1 = above threshhold, 0 = not
threshInd = cell(1,chans); % Columns = channels; {1,:} = indices above thresh; {2,:} = onset indices; {3,:} = offset indices
chk_empty = 0;
for i=1:chans
    threshInd{i} = find(abs(LFPTs.data(i,:))>thresh);
    % Checks if none of the channels exceed threshold
    chk_empty = chk_empty + numel(threshInd{i});
end

%%
if chk_empty ~= 0
    %% Find onset and offset cuttoff
    % If downsampled, reset offset relative to downsampling; onset stays that
    % same due to being so low
    if dsf > 1
        offset = offset/dsf;
    end
    for i=1:chans
        % Skip to the first usable index to avoid negative indexing
        firstIndex = find(threshInd{1,i} > onset);
        threshInd{2,i} = threshInd{1,i}(firstIndex)-onset;
        threshInd{3,i} = threshInd{1,i}(firstIndex)+offset;
    end

    %% Create overlapping intervals
    overInt = cell(2,chans); % Creates interval structure
    for i = 1:chans
        if ~isempty(threshInd{1,i})
            s = threshInd{2,i}(1);
            e = threshInd{3,i}(1);
            k = 0;
            for j=1:length(threshInd{2,i})
                if threshInd{2,i}(j)> e
                    k = k+1;
                    overInt{1,i}(k) = s;
                    overInt{2,i}(k) = e;
                    s = threshInd{2,i}(j);
                    e = threshInd{3,i}(j);
                else
                    e = max(e,threshInd{3,i}(j));
                end
            end
            %Stack last interval
            k = k+1;
            overInt{1,i}(k) = s;
            overInt{2,i}(k) = e;
        %     overInt{1,i}(k+1:end) = [];
        %     overInt{2,i}(k+1:end) = [];
        end
    end
    % Idea borrowed from MergeBrackets.m Author: Bruno Luong <brunoluong@yahoo.com> Original: 25-May-2009
    %% NaN overlapping intervals
    LFPTsNaN = LFPTs;% Copy data to be NaNed
    for i=1:chans
        if ~isempty(overInt{1,i})
            for j=1:length(overInt{1,i})
                LFPTsNaN.data(i,overInt{1,i}(j):overInt{2,i}(j)) = NaN;
            end
        end
    end
    %% Check NaN distribution across channels
    % Creates nNaN matrix of number of NaNs per channel and standard deviations
    % from the mean, and vector of channels above cutoff
    %[nNaN,indSkp] = chkNaN(LFPTsNaN,chans,NaNcutoff);
    %% Find chunks of existing data
    % Sum each column across rows to get all channels to have the same number
    % of NaNed indices per channel
    LFPTsNaN.oneChan = sum(LFPTsNaN.data,1);
%     if isempty(indSkp)
%         LFPTsNaN.oneChan = sum(LFPTsNaN.data,1); 
%     end
%     % If channel was skipped, then don't include in sum to avoid NaNing all
%     % data
%     if ~isempty(indSkp)
%         c = (1:chans);
%         c = c(c~=indSkp);
%         LFPTsNaN.oneChan = sum(LFPTsNaN.data(c,:),1);
%     end
    % Spread LFPTsNaN.oneChan across all channels
    LFPTsNaN.data(:,isnan(LFPTsNaN.oneChan)) = NaN; 
    % Find where data starts and stops
    A = LFPTsNaN.oneChan; A(~isnan(A)) = 1; A(isnan(A)) = 0;
    dataStart = find(diff(A)==1)+1;
    dataStop = find(diff(A)==-1);
    % Above misses first and last indices; check these for clean data and
    % concatenate with above vectors
    if ~isnan(LFPTsNaN.oneChan(1))
        dataStart = [1, dataStart];
    end
    if ~isnan(LFPTsNaN.oneChan(end))
        dataStop = horzcat(dataStop,length(LFPTsNaN.oneChan));
    end
    %% Check size of intervals and if < minInt, NaN
    for i = 1:length(dataStart)
        if dataStop(i) - dataStart(i) <= (minInt*adfreq+1)
            LFPTsNaN.data(:,dataStart(i):dataStop(i)) = NaN;
        end
    end
    % Overwrite LFPTs with LFPTsNaN
    LFPTs = LFPTsNaN; 
% else
%     nNaN = [];
%     indSkp = [];
end
