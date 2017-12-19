files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\toProcess','base','in');
for ii = 1:length(files)
   load(files{ii},'eventTs')
   len{ii} = eventTs.t{8}-eventTs.t{7};
end
%%
fNames = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\paper2\postBinge\','(e','in');
for ii = 1:size(fNames,2)
    names(ii,:) = strsplit(fNames{ii},'_');
end
%%
for ii = 1:size(names,1)
    disp(num2str(ii))
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\binge_notbinge\',names{ii,1},'_',names{ii,2},'_binge_vs_notbinge.mat'])
    notCoh = coh;
    notPow = psdTrls;
    notTrls = trls;
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\mat\',names{ii,1},'_',names{ii,2},'.mat'],'eventTs')
    ind = logicFind(1,cell2mat(cellfun(@(x) strncmpi('binge (s',x,8),eventTs.label,'UniformOutput',0)),'==');
    times{ii} = [eventTs.t{ind}*hist.adfreq,eventTs.t{ind+1}*hist.adfreq];
    load(['C:\Users\Pythia\Documents\GreenLab\data\paper2\postBinge\',names{ii,1},'_',names{ii,2},'_',names{ii,3}])
    postCoh = coh;
    postPow = psdTrls;
    postTrls = trls;
    clear coh psdTrls trls
    %% Check that times do not extend beyond the start of binge nor into
    % the next binge
    for k = 1:size(postPow,2)
        % Reset overlap
        overlap = [];
        for jj = 1:size(postTrls{1,k}.sampleinfo,1)
            eInd = nearest_idx3(postTrls{1,k}.sampleinfo(jj,1)-hist.eoi{k,2}(1,1)*hist.adfreq,times{ii}(:,2));
            if eInd == size(eventTs.t{1,8},1)
                overlap(jj) = (gt(eventTs.t{1,7}(eInd),eventTs.t{1,8}(eInd)+hist.eoi{k,2}(1,1)) | lt(eventTs.t{1,1}(end),eventTs.t{1,8}(eInd)+hist.eoi{k,2}(1,2)+60));
            else
                overlap(jj) = (gt(eventTs.t{1,7}(eInd),eventTs.t{1,8}(eInd)+hist.eoi{k,2}(1,1)) | lt(eventTs.t{1,7}(eInd+1),eventTs.t{1,8}(eInd)+hist.eoi{k,2}(1,2)+60));
            end
        end
        trls{1,k}.trial = postTrls{1,k}.trial(:,:,~overlap);
        trls{1,k}.time = postTrls{1,k}.time(~overlap);
        trls{1,k}.sampleinfo = postTrls{1,k}.sampleinfo(~overlap,:);
        
        psdTrls{1,k}.Pow = postPow{1,k}.Pow(:,:,~overlap);
        psdTrls{1,k}.F = hist.powF;
        psdTrls{1,k}.hammSize = hist.hammSize;
        psdTrls{1,k}.bandPow = postPow{1,k}.bandPow(:,:,~overlap);
        psdTrls{1,k}.totPow = postPow{1,k}.totPow(:,:,~overlap);
        psdTrls{1,k}.relPow = postPow{1,k}.relPow(:,:,~overlap);
        psdTrls{1,k}.Pow = postPow{1,k}.Pow(:,:,~overlap);
        
        coh{1,k}.Cxy = postCoh{1,k}.Cxy(:,:,~overlap);
        coh{1,k}.rel = postCoh{1,k}.rel(:,:,~overlap);
        coh{1,k}.band = postCoh{1,k}.band(:,:,~overlap);
        coh{1,k}.freq = postCoh{1,k}.freq;
        coh{1,k}.winSize = postCoh{1,k}.winSize;
        coh{1,k}.total = postCoh{1,k}.total(:,:,~overlap);
    end
    save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\postBinge\',names{ii,1},'_',names{ii,2},'test.mat'],'trls','psdTrls','coh')
end
