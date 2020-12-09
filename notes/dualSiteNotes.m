files = fileSearch('E:\dualSite\','.mat');
%%
fI = 1;
thisFile = files{fI};
%%
load(files{fI},'LFPTs','eventTs','adfreq')
thisFile = thisFile(1:strfind(thisFile,'.')-1);
parts = strsplit(thisFile,'_');
site = '';
%%
% Number of animals in recording
n = 4;
oldLFPTs = LFPTs;
oldEventTs = eventTs;
% Set up eventInds since box 2 and 3 are switched for stim
eventInds = [9:12;17:20;13:16;21:24];
for ii = 1:4
    LFPTs = []; eventTs = [];
    LFPTs.data = oldLFPTs.data(8*ii-8+1:8*ii,:);
    LFPTs.tvec = oldLFPTs.tvec;
    LFPTs.label = oldLFPTs.label(8*ii-8+1:8*ii);
    
    eventTs.t = oldEventTs.t(eventInds(ii,:));
    eventTs.label = oldEventTs.label(eventInds(ii,:));
    save(['E:\dualSite\splitMat\',strjoin([parts(ii),parts(n+1),site,...
        parts{n+2},parts(n+3)],'_'),'.mat'],'LFPTs','eventTs','adfreq',...
        '-v7.3')
end
%% Split single site files
files = fileSearch('E:\dualSite\','singleSite');
for fI = 1:size(files,2)
    disp(fI)
    thisFile = files{fI};
    load(files{fI},'LFPTs','eventTs','adfreq')
    thisFile = thisFile(1:strfind(thisFile,'.')-1);
    parts = strsplit(thisFile,'_');
    
    oldLFPTs = LFPTs;
    oldEventTs = eventTs;
    
    firstStart = 1; 
    firstEnd = nearest_idx3(min([oldEventTs.t{3}(1),...
        oldEventTs.t{4}(1)]),oldLFPTs.tvec)-1;
    
    secondStart = nearest_idx3(max([oldEventTs.t{1}(end),...
        oldEventTs.t{2}(end)]),oldLFPTs.tvec)+1;
    secondEnd = size(oldLFPTs.tvec,2);
    
    inds = [firstStart,firstEnd;secondStart,secondEnd];
    for ii = 1:2
        LFPTs = []; eventTs = [];
        LFPTs.tvec = oldLFPTs.tvec(inds(ii,1):inds(ii,2));
        LFPTs.data = oldLFPTs.data(:,inds(ii,1):inds(ii,2));
        LFPTs.label = oldLFPTs.label;
        
        eventTs.t = oldEventTs.t(2*ii-1:2*ii);
        eventTs.label = oldEventTs.label(2*ii-1:2*ii);
        save(['E:\dualSite\splitMat\',strjoin([parts(1),parts(2),...
            parts(ii+2),parts(5),parts(6)],'_'),'.mat'],'LFPTs',...
            'eventTs','adfreq');
    end
end
%% Update eventTs with baseline and interstim times
files = fileSearch('E:\dualSite\splitMat\','.mat');
for fI = 37:size(files,2)
    disp(num2str(fI))
    % Load file
    load(files{fI})
    % Find first non-empty eventTs.t
    eInd = logicFind(1,cellfun(@(x) ~isempty(x),eventTs.t),'==','first');
    % Use first eventTs for approximate stim times
    stimBlockOff = [logicFind(1,round(diff(eventTs.t{eInd}))>600,'=='),...
        numel(eventTs.t{eInd})];
    % Check how many stim blocks there are; 
    if numel(stimBlockOff==1)
        % If only one (singleSite) 
        stimBlockOn = [1,stimBlockOff+1];
    else
        % Otherwise (dualSite) use all but last time
        stimBlockOn = [1,stimBlockOff(1:end-1)+1];
    end
    % Convert stimBlock to LFPTs.tvec inds; jitter to ensure no stim
    % artifact is included in baselines
    stimBlockOffInds = nearest_idx3(eventTs.t{eInd}(stimBlockOff)+0.02,...
        LFPTs.tvec);
    stimBlockOnInds = nearest_idx3(eventTs.t{eInd}(stimBlockOn(1:end-1))...
        -0.03,LFPTs.tvec);
    % Convert inds to times (seconds)
    stimBlockOffTime = LFPTs.tvec(stimBlockOffInds);
    stimBlockOnTime = LFPTs.tvec(stimBlockOnInds);
    % Use stimBlock inds to get baseline inds and times
    baseOnInds = [1;stimBlockOffInds+1];
    baseOffInds = [stimBlockOnInds-1;numel(LFPTs.tvec)];

    baseOnTime = LFPTs.tvec(baseOnInds);
    baseOffTime = LFPTs.tvec(baseOffInds);
    % Use stimBlockOn and -Off to find 'clean' data between stims
    dataStartTimes = eventTs.t{eInd}(logicFind(10,round(diff(...
            eventTs.t{eInd})),'=='))+0.02;
    dataStopTimes = eventTs.t{eInd}(logicFind(10,round(diff(...
            eventTs.t{eInd})),'==')+1);
    % Split up times into the 7 blocks and add the last chunk of clean
    % data
    [starts,stops] = deal(zeros(30,numel(baseOnTime)-1));
    for ii = 1:numel(baseOnTime)-1
        starts(:,ii) = [dataStartTimes(29*ii-28:29*ii);baseOnTime(ii+1)];
        stops(:,ii) = [dataStopTimes(29*ii-28:29*ii);baseOnTime(ii+1)+10];
    end
    % Set up structure for new eventTs
    struct.label = [];
    struct.t = [];
    % Add in basline times
    for ii = 1:numel(baseOnTime)
        struct.label = [struct.label,{['Base',num2str(ii),' (Start)']},...
            {['Base',num2str(ii),' (End)']}];
        struct.t = [struct.t,{baseOnTime(ii)},{baseOffTime(ii)}];
    end
    % Add in stim times
    for ii = 1:numel(baseOnTime)-1
        struct.label = [struct.label,{['Stim',num2str(ii),' (Start)']},...
            {['Stim',num2str(ii),' (End)']}];
        struct.t = [struct.t,{starts(:,ii)},{stops(:,ii)}];
    end
    % Combine eventTs and struct
    eventTs.t = [eventTs.t,struct.t];
    eventTs.label = [eventTs.label,struct.label];
    % Save
    save(['E:\dualSite\toProcess\',files{fI}],'LFPTs','adfreq',...
        'eventTs','-v7.3')
    clearvars -except fi files
end
%% Test plot
figure
plot(LFPTs.tvec,LFPTs.data(1,:))
hold on
plot(stimBlockOnTime,zeros(1,numel(stimBlockOnTime)),'sk')
plot(stimBlockOffTime,zeros(1,numel(stimBlockOnTime)),'ok')
plot(baseOnTime,zeros(1,numel(baseOnTime)),'sr')
plot(baseOffTime,zeros(1,numel(baseOnTime)),'or')

plot(starts,ones(1,size(starts,1)).*4,'oc')
plot(stops,ones(1,size(starts,1)).*4,'sc')
%% Go through files and get timings
offsetFiles = fileSearch('E:\dualSite\processed\toUse\','.mat')';
[offset100,offset100std] = deal(zeros(size(offsetFiles,2),7));
events = cell(1,size(offsetFiles,2));
for fI = 1:size(offsetFiles,1)
    disp(num2str(fI))
    load(offsetFiles{fI})
    oC = 1;
    for ii = 6:2:18
        idx1 = nearest_idx3(hist.eventTs.t{ii},hist.eventTs.t{1});
        idx2 = nearest_idx3(hist.eventTs.t{ii},hist.eventTs.t{3});
        offset100(fI,oC) = mean(hist.eventTs.t{1}(idx1:idx1+100)-...
            hist.eventTs.t{3}(idx2:idx2+100));
        offset100std(fI,oC) = std(hist.eventTs.t{1}(idx1:idx1+100)-...
            hist.eventTs.t{3}(idx2:idx2+100));
%         offset1000(fI,c) = mean(hist.eventTs.t{1}(idx1:idx1+1000)-...
%             hist.eventTs.t{3}(idx2:idx2+1000));
        oC = oC+1;
    end
    events{fI} = hist.eventTs;
end
save('E:\dualSite\offsets.mat','offset100','offsetFiles','offset100std',...
    'events')
%% Collate dualSite data
[data,samp,files] = collateData('E:\dualSite\processed\toUse\',...
    {'dualSite'},{'pow','coh'},'trl','raw');
%% Calculate z-scores of stim data and washout data, then put into matrix
% allZ = time X feature X offset X animal X recording
%   time = all stim data (up to 60 bins max; two 5 second bins every minute
%       for 30 minutes)
%   feature = all features, 216
%   offset = all offsets, 7; ordered in ascending order
%   animal = each animal; 6 due to having both dual-site stims (NAcS-OFC & 
%       NAcS-IL)
%   recording = each recording; 3 recordings per animal per dual-site pair
%__________________________________________________________________________
% N.B. NAcS-OFC indices = [2:2:12] along fourth dimension and NAcS-IL
% indices = [1:2:12] along fourth dimension
%__________________________________________________________________________
load('E:\dualSite\dualSiteData.mat')
load('E:\dualSite\offsets.mat')
% Get sort indices for offsets in ascending order
[~,sortI] = sort(offset100,2,'ascend');
% Use the last 5 minutes of the previous baseline to z-score the stim data
% and the first 5 minutes of the next baseline (washout)
allZ = zeros(120,216,7,12,3);
% Set up recording counter
rC = 1;
% Set up animal counter
aC = 1;
for ii = 1:size(data{1},1)
    disp(num2str(ii))
    % Set up offset counter
    oC = 1;
    % Set up thisZ
    thisZ = [];
    % Find 'Base1 (End)' index in events{ii}
    eI = eventInd(events{ii},{'Base1 (End)'});
    for jj = 9:15 % Indices of the 7 stim chunks
        
        % Determine cutoff for last 5 minutes of previous baseline
        % Convert samp to time
        if ismember(ii,[4,10,23,34]) % Stupid files that don't start at t=0
            baseTime = samp{1,1}{ii,jj-8}.*0.0025+1842;
        else
            baseTime = samp{1,1}{ii,jj-8}.*0.0025;
        end
        % Find split
        baseSplit = nearest_idx3(events{ii}.t{eI(1)+2*(oC-1)}-5*60,...
            baseTime(:,1));
        
        % Calculate mean of last 5 minutes of previous (jj-8) baseline;
        % replicate to same size as stim data
        mu = mean(data{1}{ii,jj-8}(baseSplit:end,:),1);
        % Calculate standard deviation of last 5 minutes of previous (jj-8)
        % baseline; replicate to same size as stim data
        sig = std(data{1}{ii,jj-8}(baseSplit:end,:),[],1);
        % Replicate mu and sig to same size as stim data
        thisMu = repmat(mu,size(data{1}{ii,jj},1),1);
        thisSig = repmat(sig,size(data{1}{ii,jj},1),1);
        % Preallocate thisStimZ with NaNs; 60x216 (max possible)
        thisStimZ = zeros(60,216);
        thisStimZ(thisStimZ==0) = NaN;
        % Use mu and sig to calculate z-score of stim data
        thisStimZ(1:size(data{1}{ii,jj},1),:) = (data{1}{ii,jj}-mu)./sig;
        % Also build array of times corresponding to data used
        stimTime{ii,oC} = samp{1,1}{ii,jj}.*0.0025;
        
        % Determine cutoff for first 5 minutes of next baseline; washout
        if ismember(ii,[4,10,23,34]) % Stupid files that don't start at t=0
            thisWashTime = samp{1,1}{ii,jj-7}.*0.0025+1842;
        else
            thisWashTime = samp{1,1}{ii,jj-7}.*0.0025;
        end
        washSplit = nearest_idx3(events{ii}.t{eI(1)+2*oC-1}+5*60,...
            thisWashTime(:,2));
        % Replicate mu and sig to same size as washout data
        thisMu = repmat(mu,size(data{1}{ii,jj-7}(1:washSplit,:),1),1);
        thisSig = repmat(sig,size(data{1}{ii,jj-7}(1:washSplit,:),1),1);
        % Preallocate thisWashZ with NaNs; 60x216 (max possible)
        thisWashZ = zeros(60,216);
        thisWashZ(thisWashZ==0) = NaN;
        % Use mu and sig to calculate z-score of washout data
        thisWashZ(1:washSplit,:) = (data{1}{ii,jj-7}(1:washSplit,:)-...
            thisMu)./thisSig;
%         thisWashZ(1:washSplit,:) = (data{1}{ii,jj-7}(1:washSplit,:)-...
%             mu)./sig;
        % Also build array of times corresponding to data used
        washTime{ii,oC} = thisWashTime(1:washSplit,:);
        % Put thisStimZ and thisWashZ together
        thisZ = cat(3,thisZ,[thisStimZ;thisWashZ]);
        
        % Add one to counter
        oC = oC+1;
    end
    % Sort thisZ by ascending offsets using sortI
    allZ(:,:,:,aC,rC) = thisZ(:,:,sortI(ii,:));
    if rC == 3
        rC = 1;
        aC = aC+1;
    else
        rC = rC+1;
    end
end
% Split into NAcS-OFC and NAcS-IL matrices
nacs_ofc = allZ(:,:,:,2:2:12,:);
nacs_il = allZ(:,:,:,1:2:12,:);
%%
load('E:\dualSite\z.mat')
m_nacs_il = squeeze(mean(nacs_il(1:60,:,:,:,:),1,'omitnan'));
vec = reshape(m_nacs_il,1,numel(m_nacs_il));
m = mean(vec,'omitnan');
s = std(vec,'omitnan');
c = 1;
for ii = 1:0.5:6
    perc(c) = sum(abs(vec)>=(m+ii*s))/numel(m_nacs_il);
    num(c) = sum(abs(vec)>=(m+ii*s));
    c = c+1;
end
%% Find features that surpass z-score for at least 2 of the 3 recordings 
c = 1;
for ii = 1:0.5:6
    this = abs(m_nacs_il) >= ii;
    thisSum = sum(this,4)>=2;
    days(:,:,:,c) = thisSum.*sign(mean(m_nacs_il,4));
    c = c+1;
end
posThresh = squeeze(sum(m_nacs_il >= 2,4)>=2);
negThresh = squeeze(sum(m_nacs_il <= -2,4)>=2);
% Get indices of those features surviving
% features x offsets x animals x z-score threshold
[x,y,z,a] = ind2sub(size(days),logicFind(0,days,'~='));
for ii = 1:7
   figure
   imagesc(squeeze(days(181,:,4,:)))
end
%% 
for ii = 1:11
    for jj = 1:7
        for l = 1:216
            animals(l,jj,ii) = sum(days(l,jj,:,ii))>=3;
        end
    end
end
%% Try looking at distribution of z-scores
load('E:\dualSite\dualSiteData.mat')
% Just grab the NAcS-IL data
dual = data{1}{[1,2,3,7,8,9,13,14,15,19,20,21,25,26,27,31,32,33],:};
dualFiles = files([1,2,3,7,8,9,13,14,15,19,20,21,25,26,27,31,32,33]);
load('E:\dualSite\singleSiteData.mat')
single = data;
singleFiles = files;
for ii = 1:6
    
end

for ii = 1:size(dual{1},1)
    parts = strsplit(dualFiles{ii},'_');
    animal = parts{1};
    site = parts{3};
    
    baseMu = mean(dual{1,1}{ii,1},1);
    baseSig = std(dual{1,1}{ii,1},[],1);
    
    baseZ = (dual{1,1}{ii,1}-baseMu)./baseSig;
    baseP = sum(baseZ>=0,1)./size(baseZ,1);
    for jj = 1:7
        stimZ{jj}(:,:) = (dual{1,1}{ii,8+jj}-baseMu)./baseSig;
        p(ii,jj,:) = sum(stimZ{jj}>=0,1)./size(stimZ{jj},1);
    end
    singleZ = (single)
end
%%
for ii = 1:216
    for jj = 1:7
        [~,pvalue(ii,jj),chi(ii,jj),~] = prop_test(...
            [sum(stimZ{jj}(:,ii)>=0,1),sum(baseZ(:,ii)>=0,1)],...
            [size(stimZ{jj},1),size(baseZ,1)],true);
    end
end
pAdj = pvalue.*(216*7);
h = pAdj<=0.05;
for ii = 1:216
   if sum(h(ii,:))>=1
      figure
      title(['Feature ',num2str(ii)])
      bar([[baseP(ii),p(1,:,ii)]',1-[baseP(ii),p(1,:,ii)]'],'stacked')
      set(gca,'xticklabel',{'Base','-0.15','-0.1','-0.05','0','0.05','0.1','0.15'})
   end
end
%% Collate singleSite data
[data,samp,files] = collateData('E:\dualSite\processed\toUseSingle\',...
    {'IL','in';'NAcS','in';'OFC','in'},{'pow','coh'},'trl','raw');
% Go through files and grab eventTs
for ii = 1:size(files,1)
    for jj = 1:size(files{ii},2)
        load(files{ii}{jj},'hist')
        events{ii,jj} = hist.eventTs;
    end
end
allZ = [];
for ii = 1:3
    thisZ = [];
    rC = 1; aC = 1;
    for jj = 1:size(data{1},1)
        % Find split
        baseSplit = nearest_idx3(samp{ii}{jj,1}(end,1)-120000,...
            samp{ii}{jj,1}(:,1));
        % First get mean and std of baseline
        mu = mean(data{ii}{jj,1}(baseSplit:end,:),1);
        sig = std(data{ii}{jj,1}(baseSplit:end,:),[],1);
        % Preallocate thisStimZ
        thisStimZ = zeros(60,216);
        thisStimZ(thisStimZ==0) = NaN;
        % Calculate zscore of stim
        thisStimZ(1:size(data{ii}{jj,3},1),:) = (data{ii}{jj,3}-mu)./sig;
        % Grab stim times
        stimTime{ii,jj} = samp{ii}{jj,3};
        
        % Determine cutoff for first 5 minutes of next baseline; washout
        washSplit = nearest_idx3(samp{ii}{jj,2}(1,1)+120000,...
            samp{ii}{jj,2}(:,1));
        if washSplit > 60
            washSplit = 60;
        end
        washTime{ii,jj} = samp{ii}{jj,2}(1:washSplit,:).*0.0025;
        % Replicate mu and sig to same size as washout data
        thisMu = repmat(mu,size(data{ii}{jj,2}(1:washSplit,:),1),1);
        thisSig = repmat(sig,size(data{ii}{jj,2}(1:washSplit,:),1),1);
        % Preallocate thisWashZ
        thisWashZ = zeros(60,216);
        thisWashZ(thisWashZ==0) = NaN;
        % Use mu and sig to calculate z-score of washout data
        thisWashZ(1:washSplit,:) = (data{ii}{jj,2}(1:washSplit,:)-...
            thisMu)./thisSig;
        
        % Put stim and washout together
        thisZ(:,:,rC,aC) = [thisStimZ;thisWashZ];
        if rC == 3
            rC = 1;
            aC = aC+1;
        else
            rC = rC+1;
        end
    end
    allZ(:,:,:,:,ii) = thisZ;
end
%% Find features that are either consistenly above/below threshold
posThresh = squeeze(sum(mean(allZ(1:60,:,:,:,:),1,'omitnan') >= 2,3)>=2);
negThresh = squeeze(sum(mean(allZ(1:60,:,:,:,:),1,'omitnan') <= -2,3)>=2);
for ii = 1:3
   figure
   imagesc(squeeze(negThresh(:,:,ii)))
end
%% Load all files; look at change in features over time
[data,samp,files] = collateData('E:\dualSite\processed\all\',{'single';...
    'dual'},{'pow','coh'},'trl','raw');
% IRDM2
clear z day d
test = [files{1}(1:10);files{2}(1:8)];
for ii = 1:size(test,1)
    parts = strsplit(test{ii},'_');
    d(ii,1) = parts(end-1);
end
[~,ind] = sort(d);
useInds = [15,16,17,18,5,9,6,12,2,3,14,10];
% calculate mu and sigma from first baseline
mu = mean(data{1,2}{5,1},1);
sig = std(data{1,2}{5,1},[],1);
for ii = 2:size(useInds,2)
    if useInds(ii) > 10
        this = data{1,2}{useInds(ii)-10,1};
    else
        this = data{1,1}{useInds(ii),1};
    end
    z(ii,:) = mean((this-mu)./sig,1);
    day(ii) = days(datetime(d{useInds(ii)})-datetime(d{ind(1)}));
end
figure
h = imagesc(abs(z'));
colormap('viridis')
set(gca,'xtick',1:12,'xticklabel',day)
title('IRDM2')
% IRDM5
clear z day d
test = [files{1}(11:21);files{2}(9:16)];
for ii = 1:size(test,1)
    parts = strsplit(test{ii},'_');
    d(ii,1) = parts(end-1);
end
[~,ind] = sort(d);
useInds = [9,1,2,8,12,3,4,14,10,11];
mu = mean(data{1,1}{9+10,1},1);
sig = std(data{1,1}{9+10,1},[],1);
for ii = 2:size(useInds,2)
    if useInds(ii) > 10
        this = data{1,2}{useInds(ii)-10+8,1};
    else
        this = data{1,1}{useInds(ii)+10,1};
    end
    z(ii,:) = mean((this-mu)./sig,1);
    day(ii) = days(datetime(d{useInds(ii)})-datetime(d{ind(1)}));
end
figure
h = imagesc(abs(z'));
colormap('viridis')
set(gca,'xtick',1:numel(day),'xticklabel',day)
title('IRDM5')
% IRDM6
clear z day d
test = [files{1}(22:33);files{2}(17:25)];
for ii = 1:size(test,1)
    parts = strsplit(test{ii},'_');
    d(ii,1) = parts(end-1);
end
[~,ind] = sort(d);
useInds = [18,19,5,6,1,2,3,14,15,4,16,17,12];
mu = mean(data{1,2}{8+16,1},1);
sig = std(data{1,2}{8+16,1},[],1);
for ii = 2:size(useInds,2)
    if useInds(ii) > 10
        this = data{1,2}{useInds(ii)-10+16,1};
    else
        this = data{1,1}{useInds(ii)+21,1};
    end
    z(ii,:) = mean((this-mu)./sig,1);
    day(ii) = days(datetime(d{useInds(ii)})-datetime(d{ind(1)}));
end
figure
h = imagesc(abs(z'));
colormap('viridis')
set(gca,'xtick',1:numel(day),'xticklabel',day)
title('IRDM6')
% IVSA74
clear z day d
test = [files{1}(34:44);files{2}(26:34)];
for ii = 1:size(test,1)
    parts = strsplit(test{ii},'_');
    d(ii,1) = parts(end-1);
end
[~,ind] = sort(d);
useInds = [16,17,18,19,20,4,5,9,10,3,12,13,1,2,15,11];
mu = mean(data{1,2}{6+25,1},1);
sig = std(data{1,2}{6+25,1},[],1);
for ii = 2:size(useInds,2)
    if useInds(ii) > 10
        this = data{1,2}{useInds(ii)-10+25,1};
    else
        this = data{1,1}{useInds(ii)+33,1};
    end
    z(ii,:) = mean((this-mu)./sig,1);
    day(ii) = days(datetime(d{useInds(ii)})-datetime(d{ind(1)}));
end
figure
h = imagesc(abs(z'));
colormap('viridis')
set(gca,'xtick',1:numel(day),'xticklabel',day)
title('IVSA74')
% IVSA75
clear z day d
test = [files{1}(45:54);files{2}(35:40)];
for ii = 1:size(test,1)
    parts = strsplit(test{ii},'_');
    d(ii,1) = parts(end-1);
end
[~,ind] = sort(d);
useInds = [4,5,1,2,11,3,10];
mu = mean(data{1,1}{4+34,1},1);
sig = std(data{1,1}{4+34,1},[],1);
for ii = 2:size(useInds,2)
    if useInds(ii) > 10
        this = data{1,2}{useInds(ii)-10+34,1};
    else
        this = data{1,1}{useInds(ii)+44,1};
    end
    z(ii,:) = mean((this-mu)./sig,1);
    day(ii) = days(datetime(d{useInds(ii)})-datetime(d{ind(1)}));
end
figure
h = imagesc(abs(z'));
colormap('viridis')
set(gca,'xtick',1:numel(day),'xticklabel',day)
title('IVSA75')
% IVSA76
clear z day d
test = [files{1}(55:65);files{2}(41:49)];
for ii = 1:size(test,1)
    parts = strsplit(test{ii},'_');
    d(ii,1) = parts(end-1);
end
[~,ind] = sort(d);
useInds = [17,18,19,5,6,1,2,8,14,3,15,4,11];
mu = mean(data{1,2}{7+40,1},1);
sig = std(data{1,2}{7+40,1},[],1);
for ii = 2:size(useInds,2)
    if useInds(ii) > 10
        this = data{1,2}{useInds(ii)-10+40,1};
    else
        this = data{1,1}{useInds(ii)+54,1};
    end
    z(ii,:) = mean((this-mu)./sig,1);
    day(ii) = days(datetime(d{useInds(ii)})-datetime(d{ind(1)}));
end
figure
h = imagesc(abs(z'));
colormap('viridis')
set(gca,'xtick',1:numel(day),'xticklabel',day)
title('IVSA76')