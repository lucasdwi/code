function [fd1,fd2,avgCoh,relCoh,stdCoh,cohPlots] = cohComp(TFRs,eventLabel)
%% Calculate coherence
cfg = {};
cfg.method = 'coh';
 
fd1 = ft_connectivityanalysis(cfg,TFRs{1});
% Check if comparing 2 events
if length(TFRs) == 2
    fd2 = ft_connectivityanalysis(cfg,TFRs{2});
end
%% Setup notch info
fds = {fd1,fd2};

% N.B. fds freq axes need to be the same
notchInd = [nearest_idx(59,fds{1}.freq);nearest_idx(61.5,fds{1}.freq)];
%% Plot
tic
cohPlots {1} = figure('Position',[1 1 1500 500]);
cols = {[0 1 1];[1 0 0];[0 1 0]; [0 0 1]; [0 0 0]; [1 0 1]};
h = {};
lbl = {};
s = {};

for j = 1:2
    s{j} = subplot(1,3,j);
    
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
       
        title(eventLabel{j});
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
bands.thet.limit = [4,7]; bands.alph.limit = [8,13]; bands.bet.limit = [15,30]; bands.lgam.limit = [45,65]; bands.hgam.limit = [70,90];
% Find indices corresponding to frequency bands
bandNames = fieldnames(bands);
for i = 1:length(fds)
    for j = 1:length(bandNames)
        % Create frequency band intervals from indices in freq
        bands.(bandNames{j}).ind = [find(fd1.freq>bands.(bandNames{j}).limit(1),1),find(fd1.freq<bands.(bandNames{j}).limit(2),1,'last')];
        for k = 1:length(fd1.labelcmb)
            tempMean = nanmean(sq(fds{i}.cohspctrm(k,:,:)),2);
            tempStd{i}(:,j,k) = std(sq(fds{i}.cohspctrm(k,:,:)),0,2,'omitnan');
            if j == 4
                avgCoh(k,j,i) = nanmean([tempMean(bands.(bandNames{j}).ind(1):notchInd(1)-1);tempMean(notchInd(2)+1:bands.(bandNames{j}).ind(2))]);
            else
            avgCoh(k,j,i) = nanmean(tempMean(bands.(bandNames{j}).ind(1):bands.(bandNames{j}).ind(2)));
            end
        end
    end
end
%% Get pooled std
for i = 1:2
    for j = 1:length(fd1.labelcmb)
        for k = 1:length(bandNames)
            s{2,i}(j,k) = sqrt(sum((length(fd1.time)-1).*tempStd{i}(bands.(bandNames{k}).ind(1):bands.(bandNames{k}).ind(2),k,j).^2)/(length(fd1.time)*(bands.(bandNames{k}).ind(2) - bands.(bandNames{k}).ind(1) + 1)-(bands.(bandNames{k}).ind(2) - bands.(bandNames{k}).ind(1) +1)));
        end
    end
end
%% Compare relative change from baseline in total power in each band
relCoh = (avgCoh(:,:,1)-avgCoh(:,:,2))./avgCoh(:,:,2);
stdCoh = relCoh.*(s{2,1}./avgCoh(:,:,1) + s{2,2}./avgCoh(:,:,2));
%% Plot
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
