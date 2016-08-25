function [psdTrls,relPower,powerPlots,varargout] = powerComp(trlData,adfreq,eventLabel,chans,comp)
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
% Colors to use for each channel in plotting
clrs = {[0.2081 0.1663 0.5292];[0.0265 0.6137 0.8135];[0.6473 0.7456 0.4188];[0.9990 0.7653 0.2164]};

% If only one value in comp, then will look at just one event; if two
% values then will compare those two events
nTrials = {};
if length(comp) == 1
    nTrials = {length(trlData{comp}.trial)};
    trls = trlData(comp);
    events = {'event1'};
end

if length(comp) == 2
   nTrials{1} = length(trlData{1,comp(1)}.trial);
   nTrials{2} = length(trlData{1,comp(2)}.trial);
   trls = {trlData{comp(1)},trlData{comp(2)}};
   events = {'event1','event2'};
end
%%
%trl1.nTrials = length(trl1.trial); trl2.nTrials = length(trl2.trial);
%trls = {trl1,trl2};
%events = {'event1','event2'};

% behaviors = {'Approach','Binge','Rest'};
% eventLabel{1} = behaviors{event1};
% eventLabel{2} = behaviors{event2};
%%
psdTrls = {};
for i = 1:length(trls)
    for j = 1:nTrials{i}
        for k = 1:chans
            [Pxx,F] = pwelch(trls{i}.trial{1,j}(k,:),hamming(length(trls{i}.trial{1,j}(k,:))),[],length(trls{i}.trial{1,j}(k,:)),adfreq);
            % Convert PSD into dB and store
            psdTrls.(events{i}).Pow{1,j}(k,:) = (10*log10(Pxx))';
        end
    end
end
% Only need one of the frequency vectors
psdTrls.F = F;
%% Setup notch info
notchInd = [nearest_idx3(59,F);nearest_idx3(61.5,F)];
%% Plot PSDs
% Find average overall PSD and plot
powerPlots{1} = figure('Position',[1 1 1500 500]);
ax = cell(1,length(events));
for i = 1:length(events)
    psdTrls.(events{i}).Overall = mean(cat(3,psdTrls.(events{i}).Pow{1,:}),3);
    % NaN notch data
    psdTrls.(events{i}).Overall(:,notchInd(1):notchInd(2)) = NaN;
    % Check number of events and set up subplots accordingly
    if length(events) == 1
        subplot(1,2,i)
    else if length(events) == 2
    subplot(1,3,i);
        end
    end
    for j = 1:chans
        hold on; ax{i} = plot(F,psdTrls.(events{i}).Overall(j,:),'color',clrs{j});
        % Fill NaN data in plot
        thisstart = [F(notchInd(1)-1),F(notchInd(2)+1)];
        thisend = [psdTrls.(events{i}).Overall(j,notchInd(1)-1),psdTrls.(events{i}).Overall(j,notchInd(2)+1)];
        plot(thisstart,thisend,'color',clrs{j});
    end 
    if i == 1
        ylabel('Power (dB)');
    end
    xlim([0 150]);
    xlabel('Frequency (Hz)');
    title(eventLabel{comp(i)});
end
% Linkaxes of both event subplots
if length(events) == 2
linkaxes([ax{1,1}.Parent,ax{1,2}.Parent],'y');
end
% axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5,1,'\bf Power Spectral Densities','HorizontalAlignment','Center','VerticalAlignment','top');
hold off;
%% Calculate total power in frequency bands of interest
% Set frequency band limits
bands.thet.limit = [4,7]; bands.alph.limit = [8,13]; bands.bet.limit = [15 30]; bands.lgam.limit = [45,65]; bands.hgam.limit = [70,90];
% Find indices corresponding to frequency bands
bandNames = fieldnames(bands);
for i = 1:length(trls)
    for j = 1:length(bandNames)
        % Create frequency band intervals from indices in F
        bands.(bandNames{j}).ind = [find(F>=bands.(bandNames{j}).limit(1),1),find(F<=bands.(bandNames{j}).limit(2),1,'last')];
        for k = 1:nTrials{i}
            for c = 1:chans
                % Integrates across frequency bands and puts in 5x4 matrix
                % in second row;[theta; alpha; beta; low gamma; high gamma] 
                psdTrls.(events{i}).Pow{2,k}(j,c) = trapz(bands.(bandNames{j}).ind(1):bands.(bandNames{j}).ind(2),psdTrls.(events{i}).Pow{1,k}(c,bands.(bandNames{j}).ind(1):bands.(bandNames{j}).ind(2)));
                % Removes values around notch filter
                if j == 4
                    notchOut = trapz(notchInd(1):notchInd(2),psdTrls.(events{i}).Pow{1,k}(c,notchInd(1):notchInd(2)));
                    psdTrls.(events{i}).Pow{2,k}(j,c) = psdTrls.(events{i}).Pow{2,k}(j,c) - notchOut;
                end
                % Sums
                % psdTrls.(events{i}).Pow{3,k}(j,c) = sum(psdTrls.(events{i}).Pow{1,k}(c,bands.(bandNames{j}).ind(1):bands.(bandNames{j}).ind(2)));
            end
        end
    end
end
%% Calculate absolute average PSD and standard deviation for each channel across frequency bands
for i = 1:length(events)
    for j = 1:length(bandNames)
        for c = 1:chans
            psdTrls.(events{i}).Avg = mean(cat(3,psdTrls.(events{i}).Pow{2,:}),3);
            psdTrls.(events{i}).Std = std(cat(3,psdTrls.(events{i}).Pow{2,:}),0,3);
        end
    end
end
%% Normalize power in each band
% If one event: compare power of each band to total power from theta to
% high gamma minus notched data
% N.B.: the sum of each band's percent of total will NOT = 100; some data
% falls outside of the bands
if length(events) == 1
    relPower = zeros(length(bandNames),chans);
        for c = 1:chans
            psdTrls.event1.totalPower(c) = trapz(bands.thet.ind(1):notchInd(1)-1,psdTrls.event1.Overall(c,bands.thet.ind(1):notchInd(1)-1))+trapz(notchInd(2)+1:bands.hgam.ind(2),psdTrls.event1.Overall(c,notchInd(2)+1:bands.hgam.ind(2)));
            relPower(:,c) = psdTrls.event1.Avg(:,c)./psdTrls.event1.totalPower(c);
        end
end

% If two events: compare relative change from baseline in total power in each band
% Use absolute power to compute percent change; represents a reduction in
% negative power as an increase
if length(events) == 2
relPower = (psdTrls.event1.Avg-psdTrls.event2.Avg)./abs(psdTrls.event2.Avg);
stdPower = relPower.*(psdTrls.event1.Std./psdTrls.event1.Avg + psdTrls.event2.Std./psdTrls.event2.Avg);
varargout{1} = stdPower;
end

%% Plot average power differences
% If one event, use stacked bargraph to visualize percent of total power
% per frequency band per channel
if length(events) == 1
    subplot(1,2,2)
    bar(relPower'.*100,'stacked');
    set(gca,'XTickLabel',trls{1,1}.label);
    title('Distribution of power across frequency bands');
    xlabel('Channel'); ylabel('Percent of total power');
end
% Plot average percent change from event 2 to event 1 with 2 sigma confidence bars
if length(events) == 2
    subplot(1,3,3)
    h = barwitherr(stdPower.*100,relPower.*100);
    set(gca,'XTick',1:5,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'});
    l = legend(trls{comp(1)}.label);
    set(l,'Location','southeast');
    title(['Percent change in ',eventLabel{comp(1)},' from ',eventLabel{comp(2)}]);
    xlabel('Frequency Band'); ylabel(['Percent change from ',eventLabel{comp(2)}]);
    % Set color scheme
    for i = 1:chans
         h(1,i).FaceColor = clrs{i};
    end
end
tightfig(powerPlots{1});

%% T-Test: Average total band powers between events
if length(events) == 2
    powerArray1 = cat(3,psdTrls.event1.Pow{2,:});
    powerArray2 = cat(3,psdTrls.event2.Pow{2,:});
    hy = cell(length(bandNames),chans); p = cell(length(bandNames),chans);
    for i = 1:chans
        for j = 1:length(bandNames)
            [hy{j,i,1},p{j,i,1}] = ttest2(powerArray1(j,i,:),powerArray2(j,i,:));
        end
    end
end
%% Plot histogram of average total power in bands of both events
if length(events) == 2
    powerPlots{2} = figure('Position',[1 1 1226 472]);
    hBar = cell(1,chans);
    for i = 1:chans
        subplot(1,chans,i)
        [hBar{i}] = barwitherr([psdTrls.event1.Std(:,i),psdTrls.event2.Std(:,i)],[psdTrls.event1.Avg(:,i),psdTrls.event2.Avg(:,i)]);
        set(gca,'XTick',1:5,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'})
        title([trls{1}.label{i}],'FontSize',9);
        hold on
        for j = 1:length(bandNames)
            if hy{j,i,1} == 1
                if max(hBar{i}(1).YData(j),hBar{i}(2).YData(j)) > 0
                    plot(j,max(hBar{i}(1).YData(j)+psdTrls.event1.Std(j,i),hBar{i}(2).YData(j)+psdTrls.event2.Std(j,i))+0.2*(hBar{i}(1).Parent.YTick(end)/length(hBar{i}(1).Parent.YTick)),'*k');
                else
                    plot(j,min(hBar{i}(1).YData(j)-psdTrls.event1.Std(j,i),hBar{i}(2).YData(j)-psdTrls.event2.Std(j,i))+0.2*(hBar{i}(1).Parent.YTick(1)/length(hBar{i}(1).Parent.YTick)),'*k');
                end
            end
            if i == 1
                ylabel('Total Power (dBs)');
            end
        end
    end
    l = legend({eventLabel{comp(1)};eventLabel{comp(2)}}); 
    set(l,'Position',[0.909 0.828 0.089 0.095]);
    axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized','clipping','off');
    text(0.5,1,'\bf Total average power','HorizontalAlignment','Center','VerticalAlignment','top');
end
%% Propagate error from baseline average
% for i = 1:size(psdTrls.event2.Std,1)
%     for j = 1:size(psdTrls.event2.Std,2)
%         propErr(i,j) = sqrt(2*(psdTrls.event2.Std(i,j)/psdTrls.event2.Avg(i,j))^2);
%     end
% end
%% Compare all possible pairs of behavior trials
% % Get number of trials per behavior
% trl1Inds = [1:length(trl1.trial)]; trl2Inds = [1:length(trl2.trial)];
% % Setup all pairs of trials
% [A,B] = meshgrid(trl1Inds,trl2Inds);
% c = cat(2,A',B');
% perms = reshape(c,[],2);
% % Setup relative power array; percent change in event1 from baseline
% % (event2); Format: (event1power - event2power) / event2power
% relPower = {};
% for p = 1:length(perms)
%     relPower{p} = (psdTrls.event1.Pow{2,perms(p,1)} - psdTrls.event2.Pow{2,perms(p,2)})./psdTrls.event2.Pow{2,perms(p,2)};
% end
% % Get average percent and standard deviation
% avgPerc = mean(cat(3,relPower{1,:}),3);
% stdPerc = std(cat(3,relPower{1,:}),0,3);
% % Plot average percent change with error bars
% figure; h = barwitherr(stdPerc,avgPerc);
% set(gca,'XTick',1:5,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'})
% legend('Channel1','Channel2','Channel3','Channel4')
% title('Change in total power from baseline')
% % Set color scheme
% for i = 1:chans
%      h(1,i).FaceColor = clrs{i};
% end
%% Permuted 
% trl1Len = length(trl1.trial); trl2Len = length(trl2.trial);
% minTrls = min(trl1Len,trl2Len); maxTrls = max(trl1Len,trl2Len);
% relPower = {}; theseMeans = {}; theseStds = {};
% perms = 10000;
% for n = 1:perms
%     r = randperm(maxTrls,minTrls);
%     if minTrls == trl1Len
%         for i = 1:minTrls
%         psdTrls.event2.Rand{1,i} = psdTrls.event2.Pow{1,r(i)};
%         psdTrls.event2.Rand{2,i} = psdTrls.event2.Pow{2,r(i)};
%         % Calculate percent change from event2 -> event1
%         relPower{i} = (psdTrls.event1.Pow{2,i} - psdTrls.event2.Rand{2,i})./ psdTrls.event2.Rand{2,i};
%         end
%     else  
%         for i = 1:minTrls
%         psdTrls.event1.Rand{1,i} = psdTrls.event1.Pow{1,r(i)};
%         psdTrls.event1.Rand{2,i} = psdTrls.event1.Pow{2,r(i)};
%         % Calculate percent change from event2 -> event1
%         relPower{i} = (psdTrls.event1.Rand{2,i} - psdTrls.event2.Pow{2,i})./ psdTrls.event2.Pow{2,i};
%         end
%     end
%     theseMeans{n} = mean(cat(3,relPower{1,:}),3);
%     theseStds{n} = std(cat(3,relPower{1,:}),0,3);
% end
% % Average and standard deviation percent change
% thisavgPerc = mean(cat(3,theseMeans{1,:}),3);
% 
% sqStd = cellfun(@(theseStds) theseStds.^2,theseStds,'UniformOutput',false);
% thisstdPerc = sqrt(sum(cat(3,sqStd{1,:}),3)./perms);