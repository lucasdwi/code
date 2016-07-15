function [LFPTsNaN] = threshFilt(LFPTs,thresh,onset,offset,minInt)
% Find values > |thres| and indices
% Create int onset->offset
% thresh = 2.5 mV
% onset = 5 milliseconds (0.005 seconds)
% offset = 17000 milliseconds (17 seconds)
% minInt = 5 seconds
chans = size(LFPTs.data,1);
%% Create threshInd structure; 1 = above threshhold, 0 = not
threshInd = {}; % Columns = channels; {1,:} = indices above thresh; {2,:} = onset indices; {3,:} = offset indices
for i=1:chans
    threshInd{i} = find(abs(LFPTs.data(i,:))>thresh);
end
%% Find onset and offset cuttoff
for i=1:chans
    threshInd{2,i} = threshInd{1,i}-onset;
    threshInd{3,i} = threshInd{1,i}+offset;
end
%% Create overlapping intervals

overInt = {}; % Creates interval structure
for i = 1:chans
    s = threshInd{2,i}(1);
    e = threshInd{3,i}(1);
    k = 0;
    for j=1:length(threshInd{1,i})
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
end
% k = k+1;
% overInt{1,i}(k) = s;
% overInt{2,i}(k) = e;
% 
% overInt{1,i}(k+1:end) = [];
% overInt{2,i}(k+1:end) = [];
% Idea borrowed from MergeBrackets.m Author: Bruno Luong <brunoluong@yahoo.com> Original: 25-May-2009

%% NaN overlapping intervals
LFPTsNaN = LFPTs;% Copy data to be NaNed
for i=1:chans
    for j=1:length(overInt{1,i})
        LFPTsNaN.data(i,overInt{1,i}(j):overInt{2,i}(j)) = NaN;
    end
end
%% Find chunks of existing data
dataInt = {}; %Intervals of existing data
for i=1:chans
     for j=1:length(overInt{1,i})
        if ~isnan(sum(LFPTsNaN.data(i,1:(minInt*2000+1)))) %Checks for 5 secs clean data at beginning
            dataInt{1,i}(1) = 1; %Sets start to 1
            dataInt{1,i}(j+1) = overInt{2,i}(j);
        else  dataInt{1,i}(j) = overInt{2,i}(j);
        end
        dataInt{2,i}(j) = overInt{1,i}(j); %End of clean data
        dataInt{2,i}(length(overInt{1,i})+1) = overInt{1,1}(end); %Last intveral to end of data
     end
end
%% Get rid of last noise spike
% For some reason the above code consistently leaves the last noise event;
% CHECK WHEN ABLE
lastInd = {};
lastInt = zeros(2,chans);
for i = 1:chans
    lastInd{i} = find(abs(LFPTsNaN.data(i,:))>thresh);
    lastInt(1,i) = lastInd{i}(1) - onset;
    lastInt(2,i) = lastInd{i}(end) + offset;
    LFPTsNaN.data(i,lastInt(1,i):lastInt(2,i)) = NaN;
end
%% Check size of intervals and if < minInt, NaN
for i=1:chans
    for j=1:length(overInt{1,i})
        if dataInt{2,i}(j)-dataInt{1,i}(j) < (minInt*2000)
            LFPTsNaN.data(i,dataInt{1,i}(j):dataInt{2,i}(j)) = NaN;
        end
    end
end
