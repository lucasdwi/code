function [eventTs,eventLabel,clnTrls,clnEvents,trls] = trialize(eoi,eventTs,LFPTs,adfreq,minInt,chans)
%% Uses behavior marker to trialize data
% INPUTS:
% behavior = string corresponding to behavior(s) of interest; possible
%   inputs: 'full','app','binge','rest'
% eventInds = 

% If behavior of interest is entire file, 'full', create new events/times
if strcmpi(eoi(:,1),'full')
    eventTs.label = []; eventTs.label = {'Start','Stop'};
    eventTs.t = []; eventTs.t = {LFPTs.tvec(1),LFPTs.tvec(end)};
    markers = [1,1];
    eventLabel = {'All Data'};
% Otherwise run eventInd to find correct indices for events
else if sum(strcmpi(eoi(:,1),'app')) || sum(strcmpi(eoi(:,1),'binge')) || sum(strcmpi(eoi(:,1),'rest'))
       [eventInds,eventTs,eventLabel,markers] = eventInd(eventTs,eoi);
    end
end
%%
% Preallocate intTime array
numMark = size(markers,1);
for i = 1:size(markers,1)
    numMark(i) = numel(eventTs.t{markers(i,2)}); 
end
largest = max(numMark);
intTime = cell(size(markers,1),largest);
% Fill intTime array with .tvec indices
for i = 1:size(markers,1)
    tic
    for j = 1:numel(eventTs.t{markers(i,2)})
        intTime{i,j} = nearest_idx3(eventTs.t{1,markers(i,2)}(j),LFPTs.tvec):nearest_idx3(eventTs.t{1,(markers(i,2)+1)}(j),LFPTs.tvec);
        for k = 1:chans
            intTime{i,j}(k,:) = intTime{i,j}(1,:); %Copy channel one indices to all other channels
        end
    end
    toc
end
%% NaN timestamps corresponding to NaNed LFPTs data
for i = 1:size(markers,1)
    for j = 1:numel(eventTs.t{markers(i,2)})
        for k = 1:chans
            intTime{i,j}(k,isnan(LFPTs.data(i,intTime{i,j}(k,:)))) = NaN;
        end
    end
end
%% NaN contiguous data intervals less than minInt
tic
clnTrls = cell(size(markers,1),largest); %Preallocate clnTrls
thisTrls = [];
for i = 1:size(markers,1)
    for j = 1:numel(eventTs.t{markers(i,2)})
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
            
            if ~isnan(sum(intTime{i,j}(k,:))) && LFPTs.tvec(intTime{i,j}(k,end))-LFPTs.tvec(intTime{i,j}(k,1)) >= (minInt + 1/adfreq) %No NaNed data and longer than minInt + 1 for indexing
                numTrls = floor((LFPTs.tvec(intTime{i,j}(k,end))-LFPTs.tvec(intTime{i,j}(k,1)))/(minInt + 1/adfreq));
                thisTrls = intTime{i,j}(k,1:((minInt*adfreq)*numTrls));
                %thisTrls = intTime{i,j}(k,1:((minInt*adfreq+1)*numTrls));
                clnTrls{i,j}(k,:) = thisTrls;
            else for intInd = 1:length(dataStart) %Run through data intervals
                    intLen = LFPTs.tvec(dataStop(intInd)) - LFPTs.tvec(dataStart(intInd));
                    thisTrls = [];
                    if intLen >= (minInt + 1/adfreq) %Keep if big enough
                        numTrls = floor(intLen/(minInt + 1/adfreq));
                        thisTrls = horzcat(thisTrls,intTime{i,j}(k,dataStart(intInd):(dataStart(intInd)+((minInt*adfreq)*numTrls))-1));
                        %thisTrls = horzcat(thisTrls,intTime{i,j}(k,dataStart(intInd):(dataStart(intInd)+((minInt*adfreq+1)*numTrls))-1));
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
clnEvents = cell(size(markers,1),1);
for i = 1:size(markers,1)
    for j = 1:largest
        clnEvents{i} = horzcat(clnEvents{i},clnTrls{i,j});
    end
end
%% Create empty trial structure
thisTrl = [];
thisTrl.label = LFPTs.label;
thisTrl.fsample = adfreq;
thisTrl.trial = {};
thisTrl.time = {};
thisTrl.sampleinfo = [];
%thisTrl.event = [];
% Check if doing comparison, then grab correct trials
if size(eoi,1)>1
    
end
% Create trl structures for each event
trls = cell(1,size(eoi,1));
for ii = 1:size(eoi,1)
    %ind = markers(ii);
    if ~isempty(clnEvents{ii,1})
        nTrls = length(clnEvents{ii,1})/(minInt*adfreq);
        %nTrls = length(clnEvents{ii,1})/((minInt*adfreq)+1);
        for j = 1:nTrls
            % To account for the transition of zero indexed time, start
            % 1 sample after 0
            thisTrl.time{1,j} = (1/adfreq:1/adfreq:minInt);
            thisTrl.sampleinfo(j,1) = clnEvents{ii,1}(1,1+minInt*adfreq*(j-1));
            thisTrl.sampleinfo(j,2) = clnEvents{ii,1}(1,minInt*adfreq*j);
            thisTrl.trial{1,j} = LFPTs.data(1:4,clnEvents{ii,1}(1,(1+minInt*adfreq*(j-1)):minInt*adfreq*j));
            %thisTrl.trial{1,j} = LFPTs.data(1:4,clnEvents{ii,1}(1,(1+minInt*adfreq*(j-1)):1+minInt*adfreq*j));
            %thisTrl.event = ind;
        end
        trls{ii} = thisTrl;
    else
        disp('Warning: This event has 0 clean trials!')
    end
end