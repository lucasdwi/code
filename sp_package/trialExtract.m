function [clnTrls,clnEvents,trls] = trialExtract(eventInfo,eventInds,eventTs,LFPTsNaN,adfreq,minInt,dsf,chans)
%% Trialize LFPTsNaN.data based on event markers
% Inputs:
% eventInfo = structure of information about events; format: row = event;
%   column 1 = event tag (Approach = 1; Binge = 2; Rest = 3); column 2 = time
%   range of interest [first step last] in seconds; N.B. last toi will
%   usually = minInt 
% eventInds = vector of indices in eventTs of behaviors; format:
%   [integers]; e.g. [indO,indRs,indRe,indAs,indAe,indBs,indBe,indRm]; from
%   eventInd.m
% eventTs = event structure with behavior time stamps and labels; from
%   ConvertPl2All_Files.m (Stott)
% LFPTs = LFPTs data structure
% adfreq = sampling rate of data; format = Hz; from ConvertPl2All_Files.m (Stott)
% minInt = minimum interval for clean data; format: seconds
% dsf = downsample factor; format = integer
% chans = number of channels; format = integer

% Outputs:
% clnTrls = cell array of all clean trials of all behavior; row = behavior
%   (1 = Approach, 2 = Binge, 3 = Rest)
% clnEvents = collapsed clnTrls
% varargout = cell array of data structure for all events in eventInfo

% Requires output of eventInd.m
%% Creates cell array of continuous behavior epoch timestamps

% intTime{i,j}(k,l)
% i = behavior (approach, binge, rest)
% j = epoch number
% k = channel
% l = timestamp

markers = [eventInds(4), eventInds(6), eventInds(2)]; %[Approach, Binge, Rest]; gives eventTs index of behavior

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
        intTime{i,j} = nearest_idx3(eventTs.t{1,markers(i)}(j),LFPTsNaN.tvec):nearest_idx3(eventTs.t{1,(markers(i)+1)}(j),LFPTsNaN.tvec);
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
thisTrls = [];
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
            
            if ~isnan(sum(intTime{i,j}(k,:))) && LFPTsNaN.tvec(intTime{i,j}(k,end))-LFPTsNaN.tvec(intTime{i,j}(k,1)) >= (minInt + (0.0005*dsf)) %No NaNed data and longer than minInt + 1 for indexing
                numTrls = floor((LFPTsNaN.tvec(intTime{i,j}(k,end))-LFPTsNaN.tvec(intTime{i,j}(k,1)))/(minInt + (0.0005*dsf)));
                thisTrls = intTime{i,j}(k,1:((minInt*adfreq+1)*numTrls));
                clnTrls{i,j}(k,:) = thisTrls;
            else for intInd = 1:length(dataStart) %Run through data intervals
                    intLen = LFPTsNaN.tvec(dataStop(intInd)) - LFPTsNaN.tvec(dataStart(intInd));
                    thisTrls = [];
                    if intLen >= (minInt + (0.0005*dsf)) %Keep if big enough
                        numTrls = floor(intLen/(minInt + (0.0005*dsf)));
                        thisTrls = horzcat(thisTrls,intTime{i,j}(k,dataStart(intInd):(dataStart(intInd)+((minInt*adfreq+1)*numTrls))-1));
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

%% Create empty trial structure
thisTrl = [];
thisTrl.label = LFPTsNaN.label;
thisTrl.fsample = adfreq;
thisTrl.trial = {};
thisTrl.time = {};
thisTrl.sampleinfo = [];
%thisTrl.event = [];
% Create trl structures for each event
trls = cell(1,size(eventInfo,1));
for ii = 1:size(eventInfo,1)
    ind = eventInfo{ii,1};
    if ~isempty(clnEvents{ind,1})
        nTrls = length(clnEvents{ind,1})/((minInt*adfreq)+1);
        for j = 1:nTrls
            thisTrl.time{1,j} = (0:(0.0005*dsf):minInt);
            thisTrl.sampleinfo(j,1) = clnEvents{ind,1}(1,1+minInt*adfreq*(j-1));
            thisTrl.sampleinfo(j,2) = clnEvents{ind,1}(1,minInt*adfreq*j);
            thisTrl.trial{1,j} = LFPTsNaN.data(1:4,clnEvents{ind,1}(1,(1+minInt*adfreq*(j-1)):1+minInt*adfreq*j));
            %thisTrl.event = ind;
        end
        trls{ii} = thisTrl;
    else
        disp('Warning: This event has 0 clean trials!')
    end
end
%% Create trial structures
% trl1 = [];
% trl1.label = LFPTsNaN.label;
% trl1.fsample = adfreq;
% trl1.trial = {};
% trl1.time = {};
% trl1.sampleinfo = [];
% trl2 = trl1; %Copy trl1 structure for trl2
% Fill trl1 structure with extracted clean data
% ind = event1;
% if ~isempty(clnEvents{ind,1})
%     nTrls = length(clnEvents{ind,1})/((minInt*adfreq)+1);
%     for i = 1:nTrls
%         trl1.time{1,i} = (0:(0.0005*dsf):minInt);
%         trl1.sampleinfo(i,1) = clnEvents{ind,1}(1,1+minInt*adfreq*(i-1));
%         trl1.sampleinfo(i,2) = clnEvents{ind,1}(1,minInt*adfreq*i);
%         trl1.trial{1,i} = LFPTsNaN.data(1:4,clnEvents{ind,1}(1,(1+minInt*adfreq*(i-1)):1+minInt*adfreq*i));
%     end
% else
%     disp('Warning: This event has 0 clean trials!')
% end
% %% Fill trl2 structure with extracted clean data
% ind = event2;
% if ~isempty(clnEvents{ind,1})
%     nTrls = length(clnEvents{ind,1})/((minInt*adfreq)+1);
%     for i = 1:nTrls
%         trl2.time{1,i} = (0:(0.0005*dsf):minInt);
%         trl2.sampleinfo(i,1) = clnEvents{ind,1}(1,1+minInt*adfreq*(i-1));
%         trl2.sampleinfo(i,2) = clnEvents{ind,1}(1,minInt*adfreq*i);
%         trl2.trial{1,i} = LFPTsNaN.data(1:4,clnEvents{ind,1}(1,(1+minInt*adfreq*(i-1)):1+minInt*adfreq*i));
%     end
% else
%     disp('Warning: This event has 0 clean trials!')
% end
