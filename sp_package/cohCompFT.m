function [coh,cohPlots,varargout] = cohCompFT(LFPTs,trls,bands,chans,cycles,ftimwin,overlap,foi,eoi)
%% Uses output of ft_freqanalysis.m to compute coherence from cross-spectral densities

% Inputs:
% TFRs = data structures from ft_freqanalysis.m
% eventLabel = cell array of behaviors; format = string

% Outputs:
% avgCoh = matrix of average coherence values, channel-pair X frequency
%   band
% relCoh = matrix of coherence normalized to either the average coherence
%   in each pair (1 value in comp) or the first event normalized the the
%   second (2 values in comp)
% cohPlots = array of plot information
% fds = output of ft_connectivityanalysis.m with added mean and standard deviation 
% varargout = if 2 values in comp, then the error carried through from each
%   event coherence average
%% Use Fieldtrip for Fourier Analysis
% Get channel combinations
cmb = nchoosek(1:chans,2);
for c = 1:size(cmb,1)
    channelCmb(c,:) = LFPTs.label(cmb(c,:));
end
% Event 1
tic
cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = foi(1):foi(2):foi(3); % frequencies to use
% Use frequency dependent windows (n cycles per window) with (x%) overlap 
if ~isempty(cycles)
    cfg.t_ftimwin    = cycles./cfg.foi;
    minTWin          = min(cfg.t_ftimwin)*overlap;
%     cfg.toi          = eoi{1,2}(1):minTWin:eoi{1,2}(2);
    cfg.toi = eoi{1,2}(1):eoi{1,2}(2);
% Or use a constant size for windows with (x%) overlap
else
    cfg.t_ftimwin    = ones(size(cfg.foi)).*ftimwin;
    cfg.toi          = eoi{1,2}(1):ftimwin*overlap:eoi{1,2}(2);
end
cfg.keeptrials   = 'yes';
cfg.channel      = LFPTs.label;
cfg.channelcmb   = channelCmb;

TFRs{1} = ft_freqanalysis(cfg,trls{1});
toc

if size(eoi,1) == 2
    % Event2
    tic
    cfg             = []; % Create empy cfg
    cfg.output      = 'powandcsd';
    cfg.channel     = LFPTs.label;
    cfg.method      = 'mtmconvol';
    cfg.taper       = 'hanning';
    cfg.foi         = foi(1):foi(2):foi(3); % Frequencies of interest
    cfg.keeptrials  = 'yes'; % For stastical comparison
    % Use frequency dependent windows (n cycles per window) with (x%) overlap 
    if ~isempty(cycles)
        cfg.t_ftimwin    = cycles./cfg.foi;
        minTWin          = min(cfg.t_ftimwin)*overlap;
%         cfg.toi          = eoi{2,2}(1):minTWin:eoi{2,2}(2);
        cfg.toi = eoi{1,2}(1):eoi{1,2}(2);
    % Or use a constant size for windows with (x%) overlap
    else
        cfg.t_ftimwin    = ones(size(cfg.foi)).*ftimwin;
        cfg.toi          = eoi{2,2}(1):ftimwin*overlap:eoi{2,2}(2);
    end
    cfg.channelcmb  = channelCmb;

    TFRs{2} = ft_freqanalysis(cfg,trls{1,2});
    toc
end
%% Calculate coherence
cfg = {};
cfg.method = 'coh';
 
fd1 = ft_connectivityanalysis(cfg,TFRs{1});
% Check if comparing 2 events
if size(eoi,1) == 2
    fd2 = ft_connectivityanalysis(cfg,TFRs{2});
end
% Combine fds into one structure
if size(eoi,1) == 1
    fds = {fd1};
    %varargout{1} = fd1;
else if size(eoi,1) == 2
        fds = {fd1,fd2};
        %varargout{1} = fd1; varargout{2} = fd2;
    end
end
%% Setup notch info
% N.B. fds freq axes need to be the same
notchInd = [nearest_idx3(59,fds{1}.freq);nearest_idx3(61.5,fds{1}.freq)];
% Average across trials and interpolate across notch filter
for f = size(fds)
    meanCoh = nanmean(fds{f}.cohspctrm(:,:,:),3);
    sdCoh = std(fds{f}.cohspctrm(:,:,:),0,3,'omitnan');
    % NaN notch filter
    meanCoh(:,notchInd(1):notchInd(2)) = NaN; sdCoh(:,notchInd(1):notchInd(2)) = NaN;
    % Interpolate
    for c = 1:size(cmb,1)
        meanCoh(c,notchInd(1):notchInd(2)) = interp1(find(~isnan(meanCoh(c,:))),meanCoh(c,~isnan(meanCoh(c,:))),find(isnan(meanCoh(c,:))));
        sdCoh(c,notchInd(1):notchInd(2)) = interp1(find(~isnan(sdCoh(c,:))),sdCoh(c,~isnan(sdCoh(c,:))),find(isnan(sdCoh(c,:))));
    end
end
%% Plot
tic
cohPlots {1} = figure('Position',[1 1 1500 500]);
cols = {[0 1 1];[1 0 0];[0 1 0]; [0 0 1]; [0 0 0]; [1 0 1]};
h = {};
lbl = {};
s = {};

for j = 1:size(eoi,1)
    s{j} = subplot(1,length(TFRs)+1,j);
    for iCmb = 1:size(fds{j}.labelcmb,1)
        lbl{iCmb} = cat(2,fds{j}.labelcmb{iCmb,1},'-',fds{j}.labelcmb{iCmb,2});
%         temp = nanmean(sq(fds{j}.cohspctrm(iCmb,:,:)),2);
%         % NaN and interpolate notch filter
%         temp(notchInd(1):notchInd(2)) = NaN;
%         temp(notchInd(1):notchInd(2)) = interp1(find(~isnan(temp)),temp(~isnan(temp)),find(isnan(temp)),'linear');
        % Plot
        hold on;
        h{iCmb} = plot(fds{j}.freq,meanCoh(iCmb,:),'color',cols{iCmb});
        title(eoi{j,1});
        xlabel('Frequency (Hz)');
        if j == 1
            set(s{j},'Position',[0.043 0.169 0.29 0.76]);
            ylabel('Coherence');
        else
            set(s{j},'Position',[0.361 0.169 0.29 0.76]);
        end
        xlim([0 150])        
    end
end
%hL = legend(lbl,'Orientation','horizontal'); %,'Position',[0.16 0.006 0.692 0.047]);
toc
%% Calculate average coherence per frequency band
% Set frequency band limits
%bands.thet.limit = [4,7]; bands.alph.limit = [8,13]; bands.bet.limit = [15,30]; bands.lgam.limit = [45,65]; bands.hgam.limit = [70,90];
% Find indices corresponding to frequency bands
%bandNames = fieldnames(bands);

for i = 1:length(TFRs)
    %for j = 1:length(bandNames)
    for j = 1:size(bands,1)
        % Create frequency band intervals from indices in freq
        bandInd(j,1) = find(fd1.freq>=bands{j,2}(1),1);
        bandInd(j,2) = find(fd1.freq<=bands{j,2}(2),1,'last');
        %bands.(bandNames{j}).ind = [find(fd1.freq>=bands.(bandNames{j}).limit(1),1),find(fd1.freq<=bands.(bandNames{j}).limit(2),1,'last')];
    end
    %for j = 1:length(bandNames)
    for j = 1:size(bands,1)
        for k = 1:length(fd1.labelcmb)
            tempMean = nanmean(sq(fds{i}.cohspctrm(k,:,:)),2);
            tempStd{i}(:,j,k) = std(sq(fds{i}.cohspctrm(k,:,:)),0,2,'omitnan');
            if j == 4
                %avgCoh(k,j,i) = nanmean([tempMean(bands.(bandNames{j}).ind(1):notchInd(1)-1);tempMean(notchInd(2)+1:bands.(bandNames{j}).ind(2))]);
                avgCoh(k,j,i) = nanmean([tempMean(bandInd(j,1):notchInd(1)-1);tempMean(notchInd(2)+1:bandInd(j,2))]);
            else
            %avgCoh(k,j,i) = nanmean(tempMean(bands.(bandNames{j}).ind(1):bands.(bandNames{j}).ind(2)));
            avgCoh(k,j,i) = nanmean(tempMean(bandInd(j,1):bandInd(j,2)));
            end
            % Calculate average coherence in each band in each combo, not counting the
            % notch interval; creates matrix instead of column for easier
            % division later
            %totalCoh(k,j,i) = (nanmean(tempMean(bands.thet.ind(1):notchInd(1))) + nanmean(tempMean(notchInd(2):bands.hgam.ind(2))))./2;
            totalCoh(k,j,i) = (nanmean(tempMean(bandInd(1,1):notchInd(1))) + nanmean(tempMean(notchInd(2):bandInd(size(bands,1),2))))./2;
        end
    end
end
%% Get pooled std
for i = 1:length(TFRs)
    for j = 1:length(fd1.labelcmb)
        %for k = 1:length(bandNames)
        for k = 1:size(bands,1)
            %s{2,i}(j,k) = sqrt(sum((length(fd1.time)-1).*tempStd{i}(bands.(bandNames{k}).ind(1):bands.(bandNames{k}).ind(2),k,j).^2)/(length(fd1.time)*(bands.(bandNames{k}).ind(2) - bands.(bandNames{k}).ind(1) + 1)-(bands.(bandNames{k}).ind(2) - bands.(bandNames{k}).ind(1) +1)));
            s{2,i}(j,k) = sqrt(sum((length(fd1.time)-1).*tempStd{i}(bandInd(k,1):bandInd(k,2),k,j).^2)/(length(fd1.time)*(bandInd(k,2) - bandInd(k,1) + 1)-(bandInd(k,2) - bandInd(k,1) +1)));
        end
    end
end
%% Normalize coherence values
% If one condition, normalize average coherence within each frequency band
% to the average coherence from beginning of alpha to end of high gamma

% maxCoh = max(totalCoh(:,1); for Matt

if size(eoi,1) == 1
    relCoh = avgCoh./totalCoh;
    varargout{1} = [];
    %relCoh = avgCoh./maxCoh; For Matt
end
% If two conditions, find percent change from event 2 to event 1 in each
% frequency band 
if size(eoi,1) == 2
    relCoh = (avgCoh(:,:,1)-avgCoh(:,:,2))./avgCoh(:,:,2);
    stdCoh = relCoh.*(s{2,1}./avgCoh(:,:,1) + s{2,2}./avgCoh(:,:,2));
    varargout{1} = stdCoh;
end
%% Plot
if size(eoi,1) == 1
    s = subplot(1,2,2);
    bar(relCoh'.*100); colormap(cell2mat(cols));
    hline(100);
    set(gca,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'});
    title('% of maximum average coherence');
    xlabel('Frequency Band'); ylabel('Percent of average coherence'); 
    legend(lbl);
    
end
if size(eoi,1) == 2
    s = subplot(1,3,3);
    h = barwitherr(stdCoh'.*100,relCoh'.*100);
    set(gca,'XTick',1:5,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'});
    set(s,'Position',[0.7 0.169 0.29 0.76]);
    %legend(lbl);
    title('Change in coherence from baseline');
    xlabel('Frequency Band'); ylabel('Percent change from baseline');
    % Set color scheme
    for i = 1:6
         h(i).FaceColor = cols{i};
    end
    hL = legend(lbl,'Orientation','horizontal','Position',[0.2 0.002 0.604 0.055]);
end
%% Setup output structure 'coh'
coh.avgCoh = avgCoh;
coh.relCoh = relCoh;
coh.fds = fds;
coh.TFRs = TFRs;