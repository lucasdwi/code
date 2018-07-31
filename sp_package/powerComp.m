function [psdtrlData,powerPlots] = powerComp(trlData,adfreq,bands,nFilt,foi,eoi,vis)
%% Uses pwelch to compute power and plots overall PSD per channel and 
% either a distribution of power per frequency band for each channel or a 
% relative change in each frequency band 
%__________________________________________________________________________
% Inputs:
% trlData = cell array of data structure for each event; format: {1} =
%   Approach; {2} = Binge; {3} = Rest
% adfreq = sampling rate; format = Hz
% bands = cell array defining frequency band names and ranges; format =
% cell with first column as strings with band names and second column as
%   frequency ranges in Hz; e.g. ['theta',[4 10];'lgam',[45 65]]
% nFilt = frequencies to interpolate over to account for notch filter; if
%   no notch filter was applied, leave empty; format = [low high]
% foi = frequencies of interest; format = integer row vector [low step
%   high]
% eoi = events of interest
% vis = whether or not to plot power spectra; format = 'y' or 'n'
%__________________________________________________________________________
% Outputs:
% psdtrlData = power data structure for each event; contains PSD for each
%   trial of data (.Pow), average PSD (.Overall), average power within each
%   band (.Avg), standard deviation within each band (.Std), and power from
%   lowest to highest band (.totalPower)
% powerPlots = array of plot information
% hist = history structure of important 
%__________________________________________________________________________
%% LLD 2016-2018
%% Initialize
% Find if any events have 0 trials
empt = cellfun(@isempty,trlData);
% Count channels using first non-empty dataset
nChan = size(trlData{1,logicFind(0,empt,'==','first')}.trial,1);
% Count number of channels per event
nTrials = zeros(1,size(eoi,1));
for ii = 1:size(eoi,1)
    % Skip if empty
    if ~empt(ii)
        nTrials(ii) = size(trlData{1,ii}.trial,3);
    end
end
% Set nTrials to zero if no events
nTrials(empt) = 0;
% Count number of events, including potnetially empty ones
nEvents = size(eoi,1);
% Count number of bands
nBand = size(bands,1);
%% Calculate power
% Preallocate psdtrlData
psdtrlData = cell(1,size(trlData,2));
% Uses nearest power of 2 for speed.
[~,hammSize] = nearestPow2(adfreq/foi(2));
% Only calculate power spectra for events with data
for ii = logicFind(0,empt,'==')
    for jj = 1:nTrials(ii)
        for k = 1:nChan
            % Set up pwelch inputs
            data = trlData{ii}.trial(k,:,jj);
            win = hamming(hammSize);
            fVect = foi(1):foi(2):foi(3);
            % Returns two-sided Welch PSD estimate at fois
            [Pxx,F] = pwelch(data,win,[],fVect,adfreq);
            % Convert PSD into dB and store
            psdtrlData{ii}.Pow(k,:,jj) = (10*log10(Pxx))';
            % Check if notch filter needs to be interpolated over
            if ~isempty(nFilt)
                notchInd = [nearest_idx3(nFilt(1),F);...
                    nearest_idx3(nFilt(2),F)];
                notch = notchInd(1):notchInd(2);
                psdtrlData{ii}.Pow(k,notch,jj) = NaN;
                % Set up interp1 inputs
                x = find(~isnan(psdtrlData{ii}.Pow(k,:,jj)));
                v = psdtrlData{ii}.Pow(k,~isnan(...
                    psdtrlData{ii}.Pow(k,:,jj)),jj);
                xq = find(isnan(psdtrlData{ii}.Pow(k,:,jj)));
                psdtrlData{ii}.Pow(k,notch,jj) = interp1(x,v,xq,'linear');
            end
        end
    end
    % Calculate mean and standard deviation across trials
    psdtrlData{ii}.Overall = mean(psdtrlData{ii}.Pow,3);
    psdtrlData{ii}.OverallStd = std(psdtrlData{ii}.Pow,0,3);
    % Also store frequency vectors and hammSize
    psdtrlData{ii}.hammSize = hammSize;
    psdtrlData{ii}.f = F;
end
%% Plot PSDs
% Set up figure
if strcmpi(vis,'y')
    powerPlots{1} = figure;
    ax = cell(1,length(nEvents));
    for ii = 1:nEvents
        subplot(1,nEvents,ii)
        % Only plot if data exists
        if ~empt(ii)
            for jj = 1:nChan
                hold on;
                ax{ii} = plot(F,psdtrlData{ii}.Overall(jj,:));
            end
        else
            text(50,0.5,'No trials exist!','HorizontalAlignment','center')
        end
        if ii == 1
            ylabel('Power (dB)');
        end
        xlim([0 foi(3)]);
        xlabel('Frequency (Hz)');
        title(eoi(ii));
    end
else
    powerPlots = [];
end
%% Calculate total power in frequency bands of interest and normalize by 
% total power across all bands. N.B. This includes space between bands and
% so total percent will not equal 100%.
% Get frequency indices corresponding to bands of interest
bInd = bandIndices(bands,F);
% Go through each event separately, if data exist for it
for ii = logicFind(0,empt,'==')
    % Then each trial
    for k = 1:nTrials(ii)
        % And each band
        for jj = 1:nBand
            % Set up data to get trapezoidal area of
            data = psdtrlData{ii}.Pow(:,bInd(jj,1):bInd(jj,2),k);
            % Get trapezoidal area (numerical integral) under the curve
            % for each band and store (band,channel,trial)
            psdtrlData{ii}.bandPow(jj,:,k) = trapz(data,2);
        end
    end
    % Get total trapezoidal area from beginning of first band to end of
    % last. Replicate to match dimension of .bandPow (band,channel,trial)
    % where columns will be identical because the total power is the same
    % (i.e. the bands don't matter here)
    data = trapz(psdtrlData{ii}.Overall(:,bInd(1,1):bInd(end,2),:),2)';
    psdtrlData{ii}.totPow = repmat(data,size(bands,1),1,nTrials(ii));
    % Use element-wise division to obtain percent of total power across
    % bands
    psdtrlData{ii}.relPow = psdtrlData{ii}.bandPow./psdtrlData{ii}.totPow;
end
%% Get average relative power
for ii = logicFind(0,empt,'==')
    psdtrlData{ii}.avgRelPow = mean(psdtrlData{ii}.relPow,3);
    psdtrlData{ii}.stdRelPow = std(psdtrlData{ii}.relPow,0,3);
end
