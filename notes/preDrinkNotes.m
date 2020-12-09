%%
% First get subindices for animals with more channels
waterFeat = names({'ILL','CA1L','PL','SL','PR','CA1R','ILR','SR'},...
    {'d','t','a','b','lg','hg'});
alcFeat = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
% Find overlap
waInds = zeros(1,60);
for ii = 1:60
    ind = logicFind(alcFeat{ii},waterFeat,'==');
    if ~isempty(ind)
        waInds(ii) = ind;
    end
end
% Doesn't find PR-SL coherence since in waterFeat it is SL-PR
waInds(43:48) = 157:162;
% Determine overlap between alcFeat and waterAlcFeat (same features,
% different order)
waterAlcFeat = names({'PL','SL','PR','SR'},{'d','t','a','b','lg','hg'});
inds2 = zeros(1,60);
for ii = 1:60
    ind = logicFind(alcFeat{ii},waterAlcFeat,'==');
    if ~isempty(ind)
        inds2(ii) = ind;
    end
end
% Doesn't find PR-SL coherence since in waterFeat it is SL-PR
inds2(43:48) = 43:48;
files = fileSearch('F:\paper3\periDrink\','.mat')';
fStuff = cellfun(@(x) strsplit(x,'_'),files,'UniformOutput',0);
nFiles = numel(files);
watInds = [25:36,52:61,67:69];
for ii = 1:nFiles
    % Get water, alcohol, and not data/samps
    if ~any(ismember(ii,watInds))
       load(['F:\paper3\drinkNot\',fStuff{ii}{1},'_',fStuff{ii}{2},...
           '_drink_vs_~drink.mat'])
       if ~isempty(trls{1})
           alcSamps{ii} = trls{1}.sampleinfo;
       else
           alcSamps{ii} = [];
       end
       notSamps{ii} = trls{2}.sampleinfo;
       for jj = 1:numel(psdTrls)
           if ~isempty(trls{jj})
               [b,c,t] = size(psdTrls{1,jj}.relPow);
               thisPow = reshape(psdTrls{1,jj}.relPow,b*c,t)';
               [cmb,b,t] = size(coh{1,jj}.normBandCoh);
               thisCoh = reshape(permute(coh{1,jj}.normBandCoh,...
                   [2,1,3]),cmb*b,t)';
               if jj == 1
                  alcData{ii} = [thisPow,thisCoh]; 
               else
                  notData{ii} = [thisPow,thisCoh];
               end
           end
       end
    else
        load(['F:\paper3\waterAlcohol\notDrink\',fStuff{ii}{1},'_',...
            fStuff{ii}{2},'_',fStuff{ii}{3},'_~Both.mat'])
        neitherSamps{ii} = trls{1}.sampleinfo;
        [b,c,t] = size(psdTrls{1,1}.relPow);
        thisPow = reshape(psdTrls{1,1}.relPow,b*c,t)';
        [cmb,b,t] = size(coh{1,1}.normBandCoh);
        thisCoh = reshape(permute(coh{1,1}.normBandCoh,[2,1,3]),cmb*b,t)';
        neitherData{ii} = [thisPow,thisCoh];
        if exist(['F:\paper3\waterAlcohol\processed\',fStuff{ii}{1},'_',...
                fStuff{ii}{2},'_',fStuff{ii}{3},'_water_vs_alcohol.mat'])
            load(['F:\paper3\waterAlcohol\processed\',fStuff{ii}{1},'_',...
                fStuff{ii}{2},'_',fStuff{ii}{3},'_water_vs_alcohol.mat'])
            if ~isempty(trls{1})
                watSamps{ii} = trls{1}.sampleinfo;
                [b,c,t] = size(psdTrls{1,1}.relPow);
                thisPow = reshape(psdTrls{1,1}.relPow,b*c,t)';
                [cmb,b,t] = size(coh{1,1}.normBandCoh);
                thisCoh = reshape(permute(coh{1,1}.normBandCoh,...
                    [2,1,3]),cmb*b,t)';
                theseData = [thisPow,thisCoh];
                theseData = theseData(:,inds2);
                waterData{ii} = theseData;
            else
                waterData{ii} = [];
                watSamps{ii} = [];
            end
            if ~isempty(psdTrls{1,2})
                alcSamps{ii} = trls{1,2}.sampleinfo;
                [b,c,t] = size(psdTrls{1,2}.relPow);
                thisPow = reshape(psdTrls{1,2}.relPow,b*c,t)';
                [cmb,b,t] = size(coh{1,2}.normBandCoh);
                thisCoh = reshape(permute(coh{1,2}.normBandCoh,...
                    [2,1,3]),cmb*b,t)';
                theseData = [thisPow,thisCoh];
                theseData = theseData(:,inds2);
                alcData{ii} = theseData;
            else
                alcData{ii} = [];
                alcSamps{ii} = [];
            end
            % Combine neither and water into not
            notSamps{ii} = [neitherSamps{ii};watSamps{ii}];
            notData{ii} = [neitherData{ii};waterData{ii}];
        else
            notSamps{ii} = neitherSamps{ii};
            notData{ii} = neitherData{ii};
        end
    end
end
%% Go through preSamps and remove corresponding notData
% Reset inds and go through files
inds = [];
past = cell(69,26);
allPre = []; allPreID = []; allNot = []; allNotID = [];
for k = 1:nFiles
    load(['F:\paper3\periDrink\',files{k}])
    %     for ii = 1:size(trls{1,31}.sampleinfo,1)
    if any(k == [1:24,37:51,62:66])
        alcInd = 6;
    else
        alcInd = 4;
    end
    % Get alc and water samples for this recordings
    theseAlcSamps = [];
    if ~isempty(alcSamps{k})
        for ii = 1:size(alcSamps{k},1)
            theseAlcSamps = [theseAlcSamps,...
                alcSamps{k}(ii,1):alcSamps{k}(ii,2)];
        end
    end
    theseWaterSamps = [];
    if ~isempty(watSamps{k})
        for ii = 1:size(watSamps{k},1)
            theseWaterSamps = [theseWaterSamps,...
                watSamps{k}(ii,1):watSamps{k}(ii,2)];
        end
    end
    % Go through each possible drink event
    for ii = 1:size(hist.eventTs.t{1,alcInd},1)
        % And go back through time, checking if any data came from that
        % epoch
        for jj = 1:26
            if ~isempty(trls{1,27-jj})
                if jj > 1 && overlapPre{k,28-jj}(ii) == 1
                    overlapPre{k,27-jj}(ii) = 1;
                else
                    % Get index of corresponding pre-window for each drink
                    % event going backwards through time
                    dummy = logicFind(...
                        nearest_idx3(hist.eventTs.t{1,alcInd}(ii),...
                        LFPTs.tvec)-(2000+400*(jj-1)),...
                        trls{1,27-jj}.sampleinfo(:,1),'==');
                    if isempty(dummy)
                        inds{k}(ii,jj) = 0;
                    else
                        inds{k}(ii,jj) = dummy;
                        % Get samples for pre window
                        preSamps = trls{27-jj}.sampleinfo(...
                            inds{k}(ii,jj),1):...
                            trls{27-jj}.sampleinfo(inds{k}(ii,jj),2);
                    end
                    % Check for overlap with notSamps (remove pre-drinking
                    % data from notData)
                    for m = 1:size(notSamps{k},1)
                        if any(ismember(preSamps,notSamps{k}(m,1):...
                                notSamps{k}(m,2)))
                            overlap{k}(m,27-jj,ii) = 1;
                        else
                            overlap{k}(m,27-jj,ii) = 0;
                        end
                    end
                    % Check for overlap with other alcSamps (remove
                    % pre-drinking data that overlaps with previous
                    % drinking)
                    if any(ismember(preSamps,theseAlcSamps))
                        overlapPre{k,27-jj}(ii) = 1;
                    else
                        if any(ismember(preSamps,theseWaterSamps))
                            overlapPre{k,27-jj}(ii) = 2;
                        else
                            overlapPre{k,27-jj}(ii) = 0;
                        end
                    end
                end
            else
                inds{k}(ii,jj) = 0;
                overlapPre{k,27-jj}(ii) = 3;
            end
            % Grab pre-drinking data and store in past
            if overlapPre{k,27-jj}(ii) == 0 && inds{k}(ii,jj) ~= 0
                [b,c,~] = size(psdTrls{1,27-jj}.relPow);
                thisPow = reshape(psdTrls{1,27-jj}.relPow(:,:,...
                    inds{k}(ii,jj)),b*c,1)';
                [cmb,b,~] = size(coh{1,27-jj}.normBandCoh);
                thisCoh = reshape(permute(coh{1,27-jj}.normBandCoh(:,:,...
                    inds{k}(ii,jj)),[2,1,3]),cmb*b,1)';
                theseData = [thisPow,thisCoh];
                % If a waterAlcohol recording, get rid of extra channels
                % and reorder features
                if alcInd == 4
                theseData = theseData(:,waInds);
                end
                if jj == 1
                    allPre = [allPre;theseData];
                    allPreID = [allPreID;repmat(fStuff{k}(1),...
                        size(theseData,1),1)];
                else
                    past{k,27-jj} = [past{k,jj};thisPow,thisCoh];
                end
            end
        end
    end
    % Add to pastID
    pastID{k} = fStuff{k}{1};
    % Add not data
    allNot = [allNot;...
        notData{k}(~any(any(overlap{k},3),2),:)];
    allNotID = [allNotID;fStuff{k}(1)];
end
%%
files = fileSearch('F:\paper3\periDrink\','.mat');
fStuff = cellfun(@(x) strsplit(x,'_'),files,'UniformOutput',0);
ids = cellfun(@(x) x{1},fStuff,'UniformOutput',0);
allPre = []; preID = [];
allNot = []; notID = [];
nFiles = numel(files);
[past,pastID] = deal(cell(nFiles,60));
preSamps = [];
for ii = 1:nFiles
    cd('F:\paper3\periDrink\')
    load(files{ii})
    for jj = 1:26
        [b,c,t] = size(psdTrls{1,jj}.relPow);
        thisPow = reshape(psdTrls{1,jj}.relPow,b*c,t)';
        [cmb,b,t] = size(coh{1,jj}.normBandCoh);
        thisCoh = reshape(permute(coh{1,jj}.normBandCoh,[2,1,3]),cmb*b,t)';
        if jj == 26
            allPre = [allPre;thisPow,thisCoh]; %#ok
            preID = [preID;repmat(ids(ii),size(thisPow,1),1)]; %#ok
        else
            past{ii,jj} = [thisPow,thisCoh];
            preSamps = [preSamps;trls{jj}.sampleinfo]; %#ok
        end
    end
    
    % Get samps for alcohol and all other (not drink and water)
    if ii <= 24
        load(['F:\paper3\drinkNot\',fStuff{ii}{1},'_',fStuff{ii}{2},'_drink_vs_~drink.mat'])
        if ~isempty(trls{1})
            alcSamps = trls{1}.sampleinfo;
        else
            alcSamps = [];
        end
        notSamps = trls{2}.sampleinfo;
    else
       load(['F:\paper3\waterAlcohol\notDrink\',fStuff{ii}{1},'_',fStuff{ii}{2},'_~Both.mat'])
       neitherSamps = trls{1}.sampleinfo;
       load(['F:\paper3\waterAlcohol\processed\',fStuff{ii}{1},'_',fStuff{ii}{2},'_water_vs_alcohol.mat'])
       watSamps = trls{1}.sampleinfo;
       alcSamps = trls{1}.sampleinfo;
       % Combine neither and water into not
       notSamps = [nethierSamps;watSamps];
    end
    % Turn start and stop values into complete vectors and remove repeats
    allNot = []; allAlc = [];
    for k = 1:size(notSamps,1)
        allNot = [allNot,notSamps(k,1):notSamps(k,2)];
    end
    allNot = unique(allNot);
    for k = 1:size(alcSamps,1)
        allAlc = [allAlc,alcSamps(k,1):alcSamps(k,2)];
    end
    allAlc = unique(allAlc);
    % Compare to samples associated with preDrink behavior and remove 
    
end
%% Get fileNames of preDrink data
fNames = fileSearch('F:\paper3\preDrink\','.mat');
for ii = 1:size(fNames,2)
    names(ii,:) = strsplit(fNames{ii},'_');
end
%% Combine predrink data with non-overlapping non-drink data
for ii = 1:size(names,1)
    disp(num2str(ii))
    load(['F:\paper3\drinkNot\',names{ii,1},'_',names{ii,2},'_drink_vs_~drink.mat'])
    notCoh = coh;
    notPow = psdTrls;
    notTrls = trls;
    load(['F:\paper3\mat\',names{ii,1},'_',names{ii,2},'.mat'],'eventTs')
    ind = logicFind(1,cell2mat(cellfun(@(x) strncmpi('drink (s',x,8),eventTs.label,'UniformOutput',0)),'==');
    times{ii} = [eventTs.t{ind}*hist.adfreq,eventTs.t{ind+1}*hist.adfreq];
    load(['F:\paper3\preDrink\',names{ii,1},'_',names{ii,2},'_',names{ii,3}])
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
            trls{1,k}.eoi = hist.eoi{k,2};
            
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
    save(['F:\paper3\preDrinkCombined3\',names{ii,1},'_',names{ii,2},'.mat'],'trls','psdTrls','coh')
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
