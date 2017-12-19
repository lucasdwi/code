function [psdTrls,powerPlots,hist] = powerComp(trlData,adfreq,chans,bands,filter,foi,eoi,disp)
%% Uses pwelch to compute power and plots overall PSD per channel and 
% either a distribution of power per frequency band for each channel or a 
% relative change in each frequency band 

% Inputs:
% trlData = cell array of data structure for each event; format: {1} =
%   Approach; {2} = Binge; {3} = Rest
% adfreq = sampling rate; format = Hz
% eventLabel = cell array of behaviors; format = string
% chans = number of channels; format = integer
% comp = events to analyze; format = integer or integer-pair of event tags
%   (Approach = 1; Binge = 2; Rest = 3) 


% Outputs:
% psdTrls = power data structure for each event; contains PSD for each
%   trial of data (.Pow), average PSD (.Overall), average power within each
%   band (.Avg), standard deviation within each band (.Std), and power from
%   lowest to highest band (.totalPower)
% relPower = either percent of total power in each band per channel (1
%   value in comp) or relative change in power in each band per channel
%   between 2 events (2 values in comp)
% powerPlots = array of plot information
% varargout = if 2 values in comp, then carries error from both

%% Initialize
% If only one value in comp, then will look at just one event; if two
% values then will compare those two events
nTrials = {};
if size(eoi,1) == 1
    nTrials{1} = size(trlData{1,1}.trial,3);
    trls = trlData(1);
    events = {'event1'};
end

if size(eoi,1) == 2
   nTrials{1} = size(trlData{1,1}.trial,3);
   nTrials{2} = size(trlData{1,2}.trial,3);
   trls = {trlData{1,1},trlData{1,2}};
   events = {'event1','event2'};
end

if size(eoi,1) > 2
    for ii = 1:size(eoi,1)
        nTrials{ii} = size(trlData{1,ii}.trial,3);
%         nTrials{2} = size(trlData{1,2}.trial,3);
    end
    events = 1:size(eoi,1);
    trls = trlData;
end
%%
%trl1.nTrials = length(trl1.trial); trl2.nTrials = length(trl2.trial);
%trls = {trl1,trl2};
%events = {'event1','event2'};

% behaviors = {'Approach','Binge','Rest'};
% eventLabel{1} = behaviors{event1};
% eventLabel{2} = behaviors{event2};
%%
psdTrls = cell(1,size(trls,2));
% Determine hammSize from frequency step in foi as desired freq resolution
% hammSize = adfreq/foi(2);
% Uses nearest power of 2 for speed.
[~,hammSize] = nearestPow2(adfreq/foi(2));
% Check hammSize vs. trial size
%if hammSize > length(trls)
for ii = 1:length(trls)
    for j = 1:nTrials{ii}
        for k = 1:chans
            % Uses 1024 nfft
            %[Pxx,F] = pwelch(trls{ii}.trial{1,j}(k,:),hamming(hammSize),[],1024,adfreq);
            % Uses nfft of the next power of 2 of Hamming window size
            %[Pxx,F] = pwelch(trls{ii}.trial{1,j}(k,:),hamming(hammSize),[],2^nextpow2(hammSize),adfreq);
            % Uses full trial length nfft
            %[Pxx,F] = pwelch(trls{ii}.trial{1,j}(k,:),hamming(hammSize),[],length(trls{ii}.trial{1,j}(k,:)),adfreq);
            % Returns two-sided Welch PSD estimate at fois
            [Pxx,F] = pwelch(trls{ii}.trial(k,:,j),hamming(hammSize),[],(foi(1):foi(2):foi(3)),adfreq);
            % Convert PSD into dB and store
            %             psdTrls.(events{ii}).Pow{1,j}(k,:) = (10*log10(Pxx))';
            psdTrls{ii}.Pow(k,:,j) = (10*log10(Pxx))';
            
            if strcmpi(filter,'y')
                notchInd = [nearest_idx3(57.5,F);nearest_idx3(62.5,F)];
                %             psdTrls.(events{ii}).Pow{1,j}(:,notchInd(1):notchInd(2)) = NaN;
                psdTrls{ii}.Pow(k,notchInd(1):notchInd(2),j) = NaN;
                psdTrls{ii}.Pow(k,notchInd(1):notchInd(2),j) = interp1(find(~isnan(psdTrls{ii}.Pow(k,:,j))),psdTrls{ii}.Pow(k,~isnan(psdTrls{ii}.Pow(k,:,j)),j),find(isnan(psdTrls{ii}.Pow(k,:,j))),'linear');
                %             for c = 1:chans
                %                 psdTrls.(events{ii}).Pow{1,j}(c,notchInd(1):notchInd(2)) = interp1(find(~isnan(psdTrls.(events{ii}).Pow{1,j}(c,:))),psdTrls.(events{ii}).Pow{1,j}(c,~isnan(psdTrls.(events{ii}).Pow{1,j}(c,:))),find(isnan(psdTrls.(events{ii}).Pow{1,j}(c,:))),'linear');
                %             end
            end
        end
    end
end
% Also store frequency vectors and hammSize
hist.powF = F;
hist.hammSize = hammSize;
%% Plot PSDs
% Setup notch info and interpolate over data
% notchInd = [nearest_idx3(57.5,F);nearest_idx3(62.5,F)];

% Set up figure
if strcmpi(disp,'y')
    powerPlots{1} = figure('Position',[1 1 1500 500]);
    ax = cell(1,length(events));
end
% Find average overall PSD and plot
for ii = 1:length(events)
    psdTrls{ii}.Overall = mean(psdTrls{ii}.Pow,3);
    psdTrls{ii}.OverallStd = std(psdTrls{ii}.Pow,0,3);
    % NaN and interpolate over notch filter
    %     if strcmpi(filter,'y')
    %         psdTrls.(events{ii}).Overall(:,notchInd(1):notchInd(2)) = NaN;
    %         for c = 1:chans
    %              psdTrls.(events{ii}).Overall(c,notchInd(1):notchInd(2)) = interp1(find(~isnan(psdTrls.(events{ii}).Overall(c,:))),psdTrls.(events{ii}).Overall(c,~isnan(psdTrls.(events{ii}).Overall(c,:))),find(isnan(psdTrls.(events{ii}).Overall(c,:))),'linear');
    %         end
    %     end
    if strcmpi(disp,'y')
        % Check number of events and set up subplots accordingly
        if length(events) == 1
            subplot(1,2,ii)
        elseif length(events) == 2
            subplot(1,3,ii)
        elseif length(events) > 2
            subplot(2,ceil(length(events)/2),ii)
        end
        for j = 1:chans
            hold on;
            ax{ii} = plot(F,psdTrls{ii}.Overall(j,:));
        end
        if ii == 1
            ylabel('Power (dB)');
        end
        xlim([0 foi(3)]);
        xlabel('Frequency (Hz)');
        title(eoi(ii));
        % Linkaxes of both event subplots
        if length(events) == 2 && ii == 2
            linkaxes([ax{1,1}.Parent,ax{1,2}.Parent],'y');
        end
    end
end
%% Calculate total power in frequency bands of interest and normalize by 
% total power across all bands. N.B. This includes space between bands and
% so total percent will not equal 100%.
% Get frequency indices corresponding to bands of interest
bInd = bandIndices(bands,F);
% Go through each event separately
for ii = 1:length(trls)
    % Then each trial
    for k = 1:nTrials{ii}
        % And each band
        for j = 1:size(bInd,1)
            % Get trapezoidal area (numerical integral) under the curve
            % for each band and store (band,channel,trial)
            psdTrls{ii}.bandPow(j,:,k) = trapz(psdTrls{ii}.Pow(:,bInd(j,1):bInd(j,2),k),2);
        end
        % Get total trapezoidal area from beginning of first band to end of
        % last. Replicate to match dimension of .bandPow 
        % (band,channel,trial) where columns will be identical because the
        % total power is the same (i.e. the bands don't matter here)
%         psdTrls{ii}.totPow(:,:,k) = repmat(trapz(psdTrls{ii}.Pow(:,bInd(1,1):bInd(end,2),k),2)',size(bands,1),1,1);
    end
    psdTrls{ii}.totPow = repmat(trapz(psdTrls{ii}.Overall(:,bInd(1,1):bInd(end,2),:),2)',size(bands,1),1,nTrials{ii});
    % Use element-wise division to obtain percent of total power within
    % each band
    psdTrls{ii}.relPow = psdTrls{ii}.bandPow./psdTrls{ii}.totPow;
end
%% Get average relative power
for ii = 1:length(trls)
    psdTrls{ii}.avgRelPow = mean(psdTrls{ii}.relPow,3);
    psdTrls{ii}.stdRelPow = std(psdTrls{ii}.relPow,0,3);
end
%% If two behaviors, compare relPow  
if length(trls) == 2 
   % Preallocate p matrix
   p = zeros(chans,size(bInd,1));
   for iC = 1:chans
       for iB = 1:size(bInd,1)
           [~,p(iB,iC)] = ttest2(psdTrls{1,1}.relPow(iB,iC,:),psdTrls{1,2}.relPow(iB,iC,:));
       end
   end
   [~,~,pAdj] = fdr_bh(p,0.05,'dep');
end
%% Plot average power differences
if strcmpi(disp,'y')
    % Plot average percent change from event 2 to event 1 with 2 sigma confidence bars
    if length(events) == 2
        subplot(1,3,3)
        h = bar((psdTrls{1,1}.avgRelPow-psdTrls{1,2}.avgRelPow)./psdTrls{1,2}.stdRelPow);
        set(gca,'XTick',1:size(bands,1),'XTickLabel',{'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
        l = legend(trls{1,1}.label);
        set(l,'Location','northwest');
        title(['Change in ',eoi{1,1},' from ',eoi{2,1}]);
        xlabel('Frequency Band'); ylabel(['Number of standard deviations from ',eoi{2,1}]);
    end
    tightfig(powerPlots{1});
else
    powerPlots = [];
end
