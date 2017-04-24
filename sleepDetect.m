function [smInstAmp,LFPTs] = sleepDetect(LFPTs,minInt,thresh,onset,offset,adfreq,smMethod)
%[sleepInt,perc,smInstAmp,LFPTs] = sleepDetect(LFPTs,p,minInt,thresh,onset,offset,adfreq,smMethod)
%% INPUTS
% LFPTs
% minInt
% thresh
% onset
% offset
% adfreq
% smMethod
% eventTs
%%
tic
chans = size(LFPTs.data,1);
[LFPTs,~] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq,1);
%% Get NaN indices and replace with randomized data
nanInd = find(isnan(LFPTs.data(1,:)));
goodInd = find(~isnan(LFPTs.data(1,:)));
LFPTs.data(:,nanInd) = 0;
%% Filter signal for sleep, get instantaneous amplitude via Hilbert 
% transform and determine 90 %tile for cutoff for each channel
% Create filter
[b,a] = cheby1(3,.5,[8 31]/1000);
%fvtool(b,a,'Fs',2000);
% Apply filter and get %tile from amplitude distribution
for c = 1:chans
    %LFPTs.data(c,isnan(LFPTs.data(c,:))) = nanmean(LFPTs.data(c,:));
    filtData.data(c,:) = filtfilt(b,a,LFPTs.data(c,:));%(goodInd==1)));
    instAmp(c,:) = abs(hilbert(filtData.data(c,:)));
    % Smooth instAmp
    if strcmpi(smMethod,'g')
        smInstAmp(c,:) = filtfilt(gausswin(10000),1,instAmp(c,:));
    else if strcmpi(smMethod,'sg')
            smInstAmp(c,:) = sgolayfilt(instAmp(c,:),3,1001);
        end
    end
    % Get local maxima (peaks) of amplitude
    %pks{c} = findpeaks(smInstAmp(c,goodInd));
    % Split peaks into bimodal distribution
    %idx{c} = kmeans(pks{c}',2);
    %pk1{c} = pks{c}(idx{c}==1); pk2{c} = pks{c}(idx{c}==2);
    % Get distribution stats
    %[y1{c},x1{c}] = ksdensity(pk1{c}); [y2{c},x2{c}] = ksdensity(pk2{c});
    % Use second distribution mean as threshold for detection
    %perc2(c) = prctile(pk2{c},25);
    % Percentile threshold of all data
%     [f{c},xi{c}] = ksdensity(smInstAmp(c,goodInd));
%     perc(c) = prctile(smInstAmp(c,goodInd),p);
end
% % Check if sleep events seem to be happening; cut off determined by
% % examining exemplary data files
% if ~any(smInstAmp >= 3e6) == 1
%     disp('It looks like no sleep events occur in this file, review signal to make sure; skipping file.')
%     % Set up empty sleepInt
%     sleepInt = [];
%     return
% end
% % Plot instantaneous amplitude and distribution
% figure;
% for c = 1:4
%     subplot(1,2,1); h = plot((1:size(smInstAmp,2))./adfreq,smInstAmp(c,:)); hold on;
%     col{c} = get(h,'color');
%     subplot(1,2,1); hlin = hline(perc(c),'--'); set(hlin,'color',col{c}); hold on;
%     %subplot(1,2,1); hlin2 = hline(perc2(c)); set(hlin2,'color',col{c});
%     subplot(1,2,2); plot(xi{c},f{c},'k'); hold on;
%     %subplot(1,2,2); plot(x1{c},y1{c},'r'); hold on;
%     %subplot(1,2,2); plot(x2{c},y2{c},'b'); hold on;
% end
% %% Use amplitude cutoffs to determine time stamps of sleep events
% for c = 1:4
%    sleepInds{c} = LFPTs.tvec(smInstAmp(c,:) >= perc(c));
%    %sleepInds{c} = LFPTs.tvec(Z >= mu+2*sigma);
% end
% % Uses minInt to compute onset and offset 
% onset = 15; offset = 15;
% % Turn each time stamp into interval and then merge
% for c = 1:chans
%     % Skip to the first usable index to avoid negative indexing
%     firstIndex = find(sleepInds{1,c} > onset);
%     sleepInds{2,c}(1,:) = sleepInds{1,c}(firstIndex) - onset;
%     sleepInds{2,c}(2,:) = sleepInds{1,c}(firstIndex) + offset;
% end
% %% Merge overlapping intervals
% [overInt] = intMaker(sleepInds(2,:));
% %% Combine channels
% allChan = [];
% for c = 1:chans
%    allChan = horzcat(allChan,overInt{1,c}); 
% end
% % Sort allChan and put in cell for intMaker
% allChanSort{1,1} = sortrows(allChan',1)';
% Check if intervals exist in at least one other channel
% badInts = [];
% for ii = 1:length(allChanSort{1,1})
%     lessThisIdx = allChanSort{1,1}(1,:); lessThisIdx(ii) = 0;
%     % If the closest index is more than 5 seconds away, NaN that interval
%     if abs(allChanSort{1,1}(1,ii) - lessThisIdx(1,nearest_idx2(allChanSort{1,1}(1,ii),lessThisIdx(1,:)))) > 5
%         badInts = vertcat(badInts,ii);
%     end
% end
% % Remove all columns corresponding to bad intervals
% allChanSort{1,1}(:,badInts) = NaN;
% allChanSort{1,1} = allChanSort{1,1}(:,~isnan(allChanSort{1,1}(1,:)));
% %%
% Re-merge intervals across all channels
% [overIntAll] = intMaker(allChanSort);
% %% Check if last interval goes past end of recording and truncate
% if overIntAll{1,1}(2,end) > LFPTs.tvec(1,end)
%     overIntAll{1,1}(2,end) = LFPTs.tvec(1,end);
% end
% %%
% sleepInt = overIntAll{1,1};
% 
% % Plot channels with intervals
% for c = 1%:chans
%     figure;
%     plot(LFPTs.tvec,LFPTs.data(c,:));
%     hold on; 
%     for int = 1:length(sleepInt)
%         plot([sleepInt(1,int) sleepInt(2,int)],[-.9 -.9],'-k')
%     end
% end
toc
%%
% % Get nearest idices in data from tvec
% for c = 1:chans
%     for ii = 1:length(overInt{1,c})
%         dataInds{c}(1,ii) = nearest_idx3(overInt{1,c}(1,ii),LFPTs.tvec);
%         dataInds{c}(2,ii) = nearest_idx3(overInt{1,c}(2,ii),LFPTs.tvec);
%         % Grab just sleep data
%         powTest{c}(ii) = pwelch(LFPTs.data(c,dataInds{c}(1,ii):dataInds{c}(2,ii)),hamming(2000),[],[0 100],2000);
%     end
% end
%% Checks if the number of sleep event is the same across channels, if so 
% % then concatenates 
% if ~isequal(length(overInt{1,1}),length(overInt{1,2}),length(overInt{1,3}),length(overInt{1,4}))
% error('Warning: Inequal number of intervals across channels, check signal and onset/offset parameters.')
% end
% 
% catInt = cell(2,1);
% for c = 1:chans
%     catInt{1,1} = vertcat(catInt{1,1},overInt{1,c}(1,:));
%     catInt{2,1} = vertcat(catInt{2,1},overInt{1,c}(2,:));
% end
% %%
% % For each interval, find the earliest start and latest stop across
% % channels, then trim beginning and end to 5 seconds before/after sleep
% % event
% trim = minInt - 5; 
% for ii = 1:length(catInt{1,1});
%     % Get earliest start, then trim
%     sleepInt(1,ii) = min(catInt{1,1}(:,ii));
%     % Get latest stop, then trim
%     sleepInt(2,ii) = max(catInt{2,1}(:,ii));
% end
%%
% for c = 1:chans
%    mu(c) = nanmean(LFPTs.data(c,:));
%    sigma(c) = std(LFPTs.data(c,:),'omitnan');
%    pseudoDat(c,:) = normrnd(mu(c),sigma(c),1,length(nanInd));
%    LFPTs.data(c,nanInd) = pseudoDat(c,:);
% end
