function [avgCoh,relCoh,cohPlots,fds,varargout] = cohComp(TFRs,eventLabel,comp,bands)
%% Uses output of ft_freqanalysis.m to compute coherence from cross-spectral densities

% Inputs:
% TFRs = data structures from ft_freqanalysis.m
% eventLabel = cell array of behaviors; format = string
% comp = events to analyze; format = integer or integer-pair of event tags
%   (Approach = 1; Binge = 2; Rest = 3) 

% Outputs:
% avgCoh = matrix of average coherence values, channel-pair X frequency
%   band
% relCoh = matrix of coherence normalized to either the average coherence
%   in each pair (1 value in comp) or the first event normalized the the
%   second (2 values in comp)
% cohPlots = array of plot information
% fds = output of ft_connectivityanalysis.m 
% varargout = if 2 values in comp, then the error carried through from each
%   event coherence average
%% Calculate coherence
cfg = {};
cfg.method = 'coh';
 
fd1 = ft_connectivityanalysis(cfg,TFRs{1});
% Check if comparing 2 events
if length(comp) == 2
    fd2 = ft_connectivityanalysis(cfg,TFRs{2});
end

%% Setup notch info
if length(comp) == 1
    fds = {fd1};
    %varargout{1} = fd1;
else if length(comp) == 2
        fds = {fd1,fd2};
        %varargout{1} = fd1; varargout{2} = fd2;
    end
end

% N.B. fds freq axes need to be the same
notchInd = [nearest_idx(59,fds{1}.freq);nearest_idx(61.5,fds{1}.freq)];
%% Plot
tic
cohPlots {1} = figure('Position',[1 1 1500 500]);
cols = {[0 1 1];[1 0 0];[0 1 0]; [0 0 1]; [0 0 0]; [1 0 1]};
h = {};
lbl = {};
s = {};

for j = 1:length(comp)
    s{j} = subplot(1,length(TFRs)+1,j);
    
    for iCmb = 1:size(fds{j}.labelcmb,1)
        lbl{iCmb} = cat(2,fds{j}.labelcmb{iCmb,1},'-',fds{j}.labelcmb{iCmb,2});

        temp = nanmean(sq(fds{j}.cohspctrm(iCmb,:,:)),2);
        temp(notchInd(1):notchInd(2)) = NaN;
        h{iCmb} = plot(fds{j}.freq,temp,'color',cols{iCmb});
        % Fill NaN data in plot
        hold on;
        thisstart = [fds{j}.freq(notchInd(1)-1),fds{j}.freq(notchInd(2)+1)];
        thisend = [temp(notchInd(1)-1),temp(notchInd(2)+1)];
        plot(thisstart,thisend,'color',cols{iCmb});
       
        title(eventLabel{comp(j)});
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

if length(comp) == 1
    relCoh = avgCoh./totalCoh;
    %relCoh = avgCoh./maxCoh; For Matt
end
% If two conditions, find percent change from event 2 to event 1 in each
% frequency band 
if length(comp) == 2
    relCoh = (avgCoh(:,:,1)-avgCoh(:,:,2))./avgCoh(:,:,2);
    stdCoh = relCoh.*(s{2,1}./avgCoh(:,:,1) + s{2,2}./avgCoh(:,:,2));
    varargout{1} = stdCoh;
end
%% Plot
if length(comp) == 1
    s = subplot(1,2,2);
    bar(relCoh'.*100); colormap(cell2mat(cols));
    hline(100);
    set(gca,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'});
    title('% of maximum average coherence');
    xlabel('Frequency Band'); ylabel('Percent of average coherence'); 
    legend(lbl);
    
end
if length(comp) == 2
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
