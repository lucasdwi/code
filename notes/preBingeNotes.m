%% Get time from eoi
for ii = 1:61
   eoi(ii,:) = {'binge (s',[-4-ii 1-ii]}; 
end
test = cat(1,eoi{:,2});
time = test(:,2)-2.5;
%%
load('preBingeModelData.mat')
newTrainY = cell(20,1);
newTrainX = cell(20,1);
for ii = 1:size(allTrainY,1)
    for jj = 1:size(allTrainY,2)
        inds = logicFind(1,allTrainY{ii,jj},'==');
        newTrainY{ii} = [newTrainY{ii};allTrainY{ii,jj}(inds,1).*jj];
        newTrainX{ii} = [newTrainX{ii};allTrainX{ii,jj}(inds,:)];
    end
end
newTestY = cell(20,1);
newTestX = cell(20,1);
for ii = 1:size(allTestY,1)
    for jj = 1:size(allTestY,2)
        inds = logicFind(1,allTestY{ii,jj},'==');
        newTestY{ii} = [newTestY{ii};allTestY{ii,jj}(inds,1).*time(jj)];
        newTestX{ii} = [newTestX{ii};allTestX{ii,jj}(inds,:)];
    end
end
%% Get fileNames of preBinge data
fNames = fileSearch('C:\Users\Lucas\Desktop\GreenLab\data\paper2\preBinge\','.mat');
for ii = 1:size(fNames,2)
    names(ii,:) = strsplit(fNames{ii},'_');
end
%% Combine prebinge data with non-overlapping non-binge data
for ii = 1:size(names,1)
    disp(num2str(ii))
%     load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\preBinge4\',names{ii,1},'_',names{ii,2},'_binge (s.mat'])
%     preCoh = [preCoh,coh];
%     prePow = [prePow,psdTrls];
%     preTrls = [preTrls,trls];
    load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\binge_notbinge\',names{ii,1},'_',names{ii,2},'_binge_vs_notbinge.mat'])
    notCoh = coh;
    notPow = psdTrls;
    notTrls = trls;
    load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\mat\',names{ii,1},'_',names{ii,2},'.mat'],'eventTs')
    ind = logicFind(1,cell2mat(cellfun(@(x) strncmpi('binge (s',x,8),eventTs.label,'UniformOutput',0)),'==');
    times{ii} = [eventTs.t{ind}*hist.adfreq,eventTs.t{ind+1}*hist.adfreq];
    load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\preBinge\',names{ii,1},'_',names{ii,2},'_',names{ii,3}])
    preCoh = coh;
    prePow = psdTrls;
    preTrls = trls;
    clear coh psdTrls trls
    %% Check that pre trials do no overlap with binges
    for k = 1:size(prePow,2)
        for jj = 1:size(preTrls{1,k}.sampleinfo,1)
            for m = 1:size(times{ii},1)
                thisTrl = ismember(preTrls{1,k}.sampleinfo(jj,1):preTrls{1,k}.sampleinfo(jj,2),round(times{ii}(m,1)):round(times{ii}(m,2)));
                thisOverlap(jj,m) = any(thisTrl(:) == 1);
            end
        overlap(jj) = any(thisOverlap(jj,:)==1);
        end
%         trls{1,k}.label = preTrls{1,k}.label;
%         trls{1,k}.fsample = preTrls{1,k}.fsample;
        trls{1,k}.trial = preTrls{1,k}.trial(:,:,~overlap);
        trls{1,k}.time = preTrls{1,k}.time(~overlap);
        trls{1,k}.sampleinfo = preTrls{1,k}.sampleinfo(~overlap,:);
        
        psdTrls{1,k}.Pow = prePow{1,k}.Pow(:,:,~overlap);
        psdTrls{1,k}.F = hist.powF;
        psdTrls{1,k}.hammSize = hist.hammSize;
        psdTrls{1,k}.bandPow = prePow{1,k}.bandPow(:,:,~overlap);
        psdTrls{1,k}.totPow = prePow{1,k}.totPow(:,:,~overlap);
        psdTrls{1,k}.relPow = prePow{1,k}.relPow(:,:,~overlap);
        psdTrls{1,k}.Pow = prePow{1,k}.Pow(:,:,~overlap);
        
        coh{1,k}.Cxy = preCoh{1,k}.Cxy(:,:,~overlap);
        coh{1,k}.rel = preCoh{1,k}.rel(:,:,~overlap);
        coh{1,k}.band = preCoh{1,k}.band(:,:,~overlap);
        coh{1,k}.freq = preCoh{1,k}.freq;
        coh{1,k}.winSize = preCoh{1,k}.winSize;
        coh{1,k}.total = preCoh{1,k}.total(:,:,~overlap);
        % Reset overlap
        overlap = [];
    end
    %% Find notTrls that do not overlap with preTrls
    for jj = 1:size(notTrls{1,2}.sampleinfo,1)
        thisTrl = ismember(trls{1,1}.sampleinfo,notTrls{1,2}.sampleinfo(jj,1):notTrls{1,2}.sampleinfo(jj,2));
        overlap(jj) = any(thisTrl(:) == 1);
    end
    %% Extract non-overlapping trials
    trls{1,k+1}.label = notTrls{1,2}.label;
    trls{1,k+1}.fsample = notTrls{1,2}.fsample;
    trls{1,k+1}.trial = notTrls{1,2}.trial(:,:,~overlap);
    trls{1,k+1}.time = notTrls{1,2}.time(~overlap);
    trls{1,k+1}.sampleinfo = notTrls{1,2}.sampleinfo(~overlap,:);
    
    psdTrls{1,k+1}.Pow = notPow{1,2}.Pow(:,:,~overlap);
    psdTrls{1,k+1}.F = hist.powF;
    psdTrls{1,k+1}.hammSize = hist.hammSize;
    psdTrls{1,k+1}.bandPow = notPow{1,2}.bandPow(:,:,~overlap);
    psdTrls{1,k+1}.totPow = notPow{1,2}.totPow(:,:,~overlap);
    psdTrls{1,k+1}.relPow = notPow{1,2}.relPow(:,:,~overlap);
    psdTrls{1,k+1}.Pow = notPow{1,2}.Pow(:,:,~overlap);
    
    coh{1,k+1}.Cxy = notCoh{1,2}.Cxy(:,:,~overlap);
    coh{1,k+1}.rel = notCoh{1,2}.rel(:,:,~overlap);
    coh{1,k+1}.band = notCoh{1,2}.band(:,:,~overlap);
    coh{1,k+1}.freq = notCoh{1,2}.freq;
    coh{1,k+1}.winSize = notCoh{1,2}.winSize;
    coh{1,k+1}.total = notCoh{1,2}.total(:,:,~overlap);
    overlap = [];
    %%
    save(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\preBingeCombined\',names{ii,1},'_',names{ii,2},'.mat'],'trls','psdTrls','coh')
end
