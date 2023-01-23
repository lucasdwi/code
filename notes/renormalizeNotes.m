files = fileSearch('E:\blast\processed\control\','.mat');
for fI = 1:numel(files)
    disp(fI)
    load(['E:\blast\processed\control\',files{fI}])
    % normalize power by average power per channel across bands across all time
    bInd = bandIndices(hist.bands,psdTrls{1}.f);
    mPow = repmat(mean(psdTrls{1}.Overall(:,bInd(1,1):bInd(end,end),:,:),...
        [2,3,4]),1,6)';
    psdTrls{1}.meanRelPow = psdTrls{1}.bandPow./mPow;
    % normalize coherence by average coherence per channel combination across
    % bands across time
    bInd = bandIndices(hist.bands,coh{1}.f);
    mCoh = repmat(mean(coh{1}.Cxy(:,bInd(1,1):bInd(end,end),:,:),...
        [2,3,4]),1,6);
    coh{1}.meanRelCoh = coh{1}.mBandCoh./mCoh;
    save(['E:\blast\processed\control\',files{fI}],'-append',...
        'psdTrls','coh')
end