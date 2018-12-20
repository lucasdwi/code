%% Get fileNames of preDrink data
fNames = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper3\preDrink\','.mat');
for ii = 1:size(fNames,2)
    names(ii,:) = strsplit(fNames{ii},'_');
end
%% Combine prebinge data with non-overlapping non-drink data
for ii = 1:size(names,1)
    disp(num2str(ii))
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\drinkNot\',names{ii,1},'_',names{ii,2},'_drink_vs_~drink.mat'])
    notCoh = coh;
    notPow = psdTrls;
    notTrls = trls;
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\mat\',names{ii,1},'_',names{ii,2},'.mat'],'eventTs')
    ind = logicFind(1,cell2mat(cellfun(@(x) strncmpi('drink (s',x,8),eventTs.label,'UniformOutput',0)),'==');
    times{ii} = [eventTs.t{ind}*hist.adfreq,eventTs.t{ind+1}*hist.adfreq];
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\preDrink\',names{ii,1},'_',names{ii,2},'_',names{ii,3}])
    preCoh = coh;
    prePow = psdTrls;
    preTrls = trls;
    clear coh psdTrls trls
    %% Check that pre trials do no overlap with drinking
    overlap = [];
    k = 1;
    while ~any(overlap) &&  k <= size(prePow,2)
%         disp(num2str(k))
        if ~isempty(preTrls{1,k})
            % Reset overlap
            overlap = [];
            for jj = 1:size(preTrls{1,k}.sampleinfo,1)
                for m = 1:size(times{ii},1)
                    thisTrl = ismember(preTrls{1,k}.sampleinfo(jj,1):preTrls{1,k}.sampleinfo(jj,2),round(times{ii}(m,1)):round(times{ii}(m,2)));
                    % Force last index to 0
                    thisTrl(end) = 0;
                    thisOverlap(jj,m) = any(thisTrl(:) == 1);
                end
                overlap(jj) = any(thisOverlap(jj,:)==1);
            end
            trls{1,k}.trial = preTrls{1,k}.trial(:,:,~overlap);
            trls{1,k}.time = preTrls{1,k}.time(~overlap);
            trls{1,k}.sampleinfo = preTrls{1,k}.sampleinfo(~overlap,:);
            
            psdTrls{1,k}.Pow = prePow{1,k}.Pow(:,:,~overlap);
            psdTrls{1,k}.f = prePow{1,k}.f;
            psdTrls{1,k}.hammSize = prePow{1,k}.hammSize;
            psdTrls{1,k}.bandPow = prePow{1,k}.bandPow(:,:,~overlap);
            psdTrls{1,k}.totPow = prePow{1,k}.totPow(:,:,~overlap);
            psdTrls{1,k}.relPow = prePow{1,k}.relPow(:,:,~overlap);
            psdTrls{1,k}.Pow = prePow{1,k}.Pow(:,:,~overlap);
            
            coh{1,k}.Cxy = preCoh{1,k}.Cxy(:,:,~overlap);
            coh{1,k}.rel = preCoh{1,k}.normBandCoh(:,:,~overlap);
            coh{1,k}.band = preCoh{1,k}.mBandCoh(:,:,~overlap);
            coh{1,k}.freq = preCoh{1,k}.f;
            coh{1,k}.mtCxy = preCoh{1,k}.mtCxy(:,:,~overlap);
        end
        k = k+1;
    end
    %% Find notTrls that do not overlap with preTrls
    samps = [];
    for jj = 1:size(trls,2)
        if~isempty(trls{1,jj})
            samps = [samps;trls{1,jj}.sampleinfo];
        end
    end
    for jj = 1:size(notTrls{1,2}.sampleinfo,1)
        thisTrl = ismember(samps,notTrls{1,2}.sampleinfo(jj,1):notTrls{1,2}.sampleinfo(jj,2));
        overlap(jj) = any(thisTrl(:) == 1);
    end
    %% Extract non-overlapping trials
    trls{1,k}.label = notTrls{1,2}.label;
    trls{1,k}.fsample = notTrls{1,2}.fsample;
    trls{1,k}.trial = notTrls{1,2}.trial(:,:,~overlap);
    trls{1,k}.time = notTrls{1,2}.time(~overlap);
    trls{1,k}.sampleinfo = notTrls{1,2}.sampleinfo(~overlap,:);
    
    psdTrls{1,k}.Pow = notPow{1,2}.Pow(:,:,~overlap);
    psdTrls{1,k}.f = notPow{1,2}.f;
    psdTrls{1,k}.hammSize = notPow{1,2}.hammSize;
    psdTrls{1,k}.bandPow = notPow{1,2}.bandPow(:,:,~overlap);
    psdTrls{1,k}.totPow = notPow{1,2}.totPow(:,:,~overlap);
    psdTrls{1,k}.relPow = notPow{1,2}.relPow(:,:,~overlap);
    psdTrls{1,k}.Pow = notPow{1,2}.Pow(:,:,~overlap);
    
    coh{1,k}.Cxy = notCoh{1,2}.Cxy(:,:,~overlap);
    coh{1,k}.rel = notCoh{1,2}.normBandCoh(:,:,~overlap);
    coh{1,k}.band = notCoh{1,2}.mBandCoh(:,:,~overlap);
    coh{1,k}.freq = notCoh{1,2}.f;
    coh{1,k}.mtCxy = notCoh{1,2}.mtCxy(:,:,~overlap);
    overlap = [];
    %%
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\preDrinkCombined2\',names{ii,1},'_',names{ii,2},'.mat'],'trls','psdTrls','coh')
end
%% Grab first and last drink data point
fNames = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper3\drinkNot\','.mat');
c = 1;
for ii = 1:size(fNames,2)
    load(fNames{ii})
    % Check for empty drink trials
    if ~isempty(trls{1,1})
        inds = logicFind(1,[1;diff(trls{1,1}.sampleinfo(:,1))~=2000],'==');
        for jj = 1:numel(inds)
            % Check if last ind or if next ind is the same as the current
            % ind +1; if so then replicate data for first and last
            if jj == numel(inds) || inds(jj)+1 == inds(jj+1)
                firstPsd(:,:,c) = psdTrls{1,1}.relPow(:,:,inds(jj));
                lastPsd(:,:,c) = psdTrls{1,1}.relPow(:,:,inds(jj));
                firstCoh(:,:,c) = coh{1,1}.normBandCoh(:,:,inds(jj));
                lastCoh(:,:,c) = coh{1,1}.normBandCoh(:,:,inds(jj));
                c = c+1;
                % Otherwise, grab first and last (the next ind-1)
            else
                firstPsd(:,:,c) = psdTrls{1,1}.relPow(:,:,inds(jj));
                lastPsd(:,:,c) = psdTrls{1,1}.relPow(:,:,inds(jj+1)-1);
                firstCoh(:,:,c) = coh{1,1}.normBandCoh(:,:,inds(jj));
                lastCoh(:,:,c) = coh{1,1}.normBandCoh(:,:,inds(jj+1)-1);
                c = c+1;
            end
        end
    end
end
save('C:\Users\Pythia\Documents\GreenLab\data\paper3\firstLastDrink.mat','firstPsd','firstCoh','lastPsd','lastCoh')
%% Check postDrink doesn't overlap with drink
fNames = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper3\postDrink\','.mat');
for ii = 1:size(fNames,2)
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\drinkNot\',fNames{ii}(1:end-13),'_drink_vs_~drink.mat'],'trls')
    drink = trls;
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\postDrink\',fNames{ii}])
    if ~isempty(drink{1,1})
        for jj = 1:size(drink{1,1}.sampleinfo,1)
            times(jj*2000-1999:jj*2000) = drink{1,1}.sampleinfo(jj,1):drink{1,1}.sampleinfo(jj,2);
        end
        for jj = 1:size(trls,2)
            if ~isempty(trls{jj})
                for k = 1:size(trls{jj}.sampleinfo,1)
                    thisTrl = ismember(trls{jj}.sampleinfo(k,1):trls{jj}.sampleinfo(k,2),times);
                    overlap = any(thisTrl);
                end
                postPsd{jj}.relPow = psdTrls{jj}.relPow(:,:,~overlap);
                postCoh{jj}.normBandCoh = coh{jj}.normBandCoh(:,:,~overlap);
            end
        end
    else
        postPsd = psdTrls;
        postCoh = coh;
    end
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\postDrink2\',fNames{ii}(1:end-13),'.mat'],'postPsd','postCoh')
end
%% Power or Coh: PreBinge-Binge-PostBinge
% Channels: SL,SR,CL,CR; SLSR,SLCL,SLCR,SRCL,SRCR,CLCR
% Freq: delta,theta,alpha,beta,lgamma,hgamma
clear chan pair freq
chan = 1;
pair = [];
freq = 1;
feat = 'd';
loc = 'sl';
% Get preBinge and notBinge data
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper3\'...
    'preDrinkCombined2']);
% Preallocate
[prePow,preCoh] = deal(zeros(size(files,2),63));
[notPow,notCoh] = deal(zeros(1,size(files,2)));
prePow(prePow==0) = NaN;
preCoh(preCoh==0) = NaN;
notPow(notPow==0) = NaN;
notCoh(notCoh==0) = NaN;
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'psdTrls')
    else
        load(files{ii},'coh')
    end
    for t = 1:size(psdTrls,2)
        if ~isempty(psdTrls{1,t})
            if isempty(pair)
                prePow(ii,t) = mean(psdTrls{t}.relPow(freq,chan,:),'omitnan');
            else
                preCoh(ii,t) = mean(coh{t}.rel(pair,freq,:),'omitnan');
            end
        end
    end
    if isempty(pair)
        notPow(ii) = mean(psdTrls{1,end}.relPow(freq,chan,:),'omitnan');
    else
        notCoh(ii) = mean(coh{1,end}.rel(pair,freq,:));
    end
end

if isempty(pair)
    mPre = mean(prePow,1,'omitnan');
    sPre = std(prePow,[],1,'omitnan');
    mNot = mean(mean(notPow,1,'omitnan'));
else
    mPre = mean(preCoh,1,'omitnan');
    sPre = std(preCoh,[],1,'omitnan');
    mNot = mean(mean(notCoh,1,'omitnan'));
end

% Get drink data
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\firstLastDrink.mat')
if isempty(pair)
    mDrink = [mean(firstPsd(freq,chan,:)) mean(lastPsd(freq,chan,:))];
    sDrink = [std(firstPsd(freq,chan,:)) std(lastPsd(freq,chan,:))];
else
    mDrink = [mean(firstCoh(pair,chan,:)) mean(lastCoh(pair,chan,:))];
    sDrink = [std(firstCoh(pair,chan,:)) std(lastCoh(pair,chan,:))];
end

% Get PostBinge data
files = fileSearch(['C:\Users\Pythia\Documents\GreenLab\data\paper3\'...
    'postDrink2'],'.mat');
% Preallocate
[postPow,postCoh] = deal(zeros(size(files,2),61));
postPow(postPow==0) = NaN;
postCoh(postCoh==0) = NaN;
for ii = 1:length(files)
    if isempty(pair)
        load(files{ii},'postPsd')
    else
        load(files{ii},'postCoh')
    end
    for t = 1:61
        if isempty(pair)
            postPow(ii,t) = mean(postPsd{t}.relPow(freq,chan,:),'omitnan');
        else
            postCoh(ii,t) = mean(postCoh{t}.rel(pair,freq,:),'omitnan');
        end
    end
end
if isempty(pair)
    mPost = mean(postPow,1,'omitnan');
    sPost = std(postPow,[],1,'omitnan');
else
    mPost = mean(postCoh,1,'omitnan');
    sPost = std(postCoh,[],1,'omitnan');
end

% Plot
figure
hold on
shadedErrorBar(1:63,fliplr(mPre.*100),fliplr(sPre.*100))
scatterErr(64,mDrink(1)*100,sDrink(1)*100,0,'col',[0 0.45 0.74])
scatterErr(70,mDrink(2)*100,sDrink(2)*100,0,'col',[0 0.45 0.74])

shadedErrorBar(71:131,mPost(1:61).*100,sPost(1:61).*100)
plot(1:131,ones(1,131).*mNot*100,'k')
plot(1:131,ones(1,131).*mean(mDrink)*100,'--','color',[0 0.45 0.74])
xlim([1 131])
set(gca,'XTick',[1:10:51,63.5,70.5,81:10:131],...
    'XTickLabel',[-62.5:10:-12.5,0,0,12.5:10:62.5])
title([loc,' ',feat]);
ylabel(['% ',feat])
text(162,mean(mDrink)*100,'Drink','color',[0 0.45 0.74])
text(162,mNot*100,'Other')
xlabel('Time')
box off
%%
chan = 1;
pair = [];
freq = 1;

psd = zeros(44,32);
csd = zeros(44,32);

files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper3\periDrink\','.mat');
for ii = 1:size(files,2)
    load(files{ii})
    for jj = 1:size(trls,2)
        allTimes = [hist.eventTs.t{1,6},hist.eventTs.t{1,7}];
        for k = 1:size(trls{jj}.sampleinfo,1)
            % Get index of the event this trial is based on
            ind = nearest_idx3(trls{jj}.sampleinfo(k,1)/400+hist.eoi{jj,2}(1),hist.eventTs.t{1,6});
            % Build time vector of all other drinking excluding above trial
            times = allTimes;
            times(ind,:) = [];
            % Convert to samples
            times = times.*400;
            % Check for overlap and remove if any exists
            if any(ismember(trls{jj}.sampleinfo(k,1):trls{jj}.sampleinfo(k,2),times))
                psdTrls{jj}.relPow(:,:,k) = NaN;
                coh{jj}.normBandCoh(:,:,k) = NaN;
            end
        end
        if ~isempty(psdTrls{jj})
            if isempty(pair)
                psd(ii,jj) = mean(psdTrls{jj}.relPow(chan,freq,:));
            else
                csd(ii,jj) = mean(coh{jj}.normBandCoh(pair,freq,:));
            end
        else
            if isempty(pair)
                psd(ii,jj) = NaN;
            else
                csd(ii,jj) = NaN;        
            end
        end
    end
end
%%
chan = 1;
pair = [];
freq = 1;
pow = zeros(44,32);
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper3\periDrink\','.mat');
for ii = 1:size(files,2)
    load(files{ii})
    for jj = 1:size(psdTrls,2)
        if ~isempty(psdTrls{jj})
            pow(ii,jj) = mean(psdTrls{jj}.relPow(chan,freq,:));
        else
            pow(ii,jj) = NaN;
        end
    end
end
load('C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\drinkNotRaw.mat')
mDrink = mean(catData{1,1}(:,(chan-1)*6+freq));
mNot = mean(catData{1,2}(:,(chan-1)*6+freq));
figure
hold on
shadedErrorBar(12:16,fliplr(mean(pow(:,1:5),1,'omitnan')),fliplr(std(pow(:,1:5),[],1,'omitnan')),{'color',[0 0.45 0.74]})
shadedErrorBar(1:12,fliplr(mean(pow(:,5:16),1,'omitnan')),fliplr(std(pow(:,5:16),[],1,'omitnan')))
% shadedErrorBar(14:16,fliplr(mean(pow(:,1:3),1,'omitnan')),fliplr(std(pow(:,1:3),[],1,'omitnan')),{'color',[0 0.45 0.74]})
% shadedErrorBar(1:14,fliplr(mean(pow(:,3:16),1,'omitnan')),fliplr(std(pow(:,3:16),[],1,'omitnan')))
shadedErrorBar(17:21,mean(pow(:,17:21),1,'omitnan'),std(pow(:,17:21),[],1,'omitnan'),{'color',[0 0.45 0.74]})
shadedErrorBar(21:32,mean(pow(:,21:32),1,'omitnan'),std(pow(:,21:32),[],1,'omitnan'))
% shadedErrorBar(17:19,mean(pow(:,17:19),1,'omitnan'),std(pow(:,17:19),[],1,'omitnan'),{'color',[0 0.45 0.74]})
% shadedErrorBar(19:32,mean(pow(:,19:32),1,'omitnan'),std(pow(:,19:32),[],1,'omitnan'))
set(gca,'XTick',[2:2:16,17:2:32],'XTickLabel',[-11.5:2:2.5,-2.5:2:12.5])
xtickangle(90)
plot([1 32],[mDrink mDrink],'--','color',[0 0.45 0.74])
plot([1 32],[mNot mNot],'--k')
