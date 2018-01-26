function [coh,cohPlots] = cohCompMatOld(LFPTs,chans,trls,foi,winSize,overlap,adfreq,bands,zeroedChannel,eoi)
%% Cycle through trials and calculate coherence
cmb = nchoosek(1:chans,2);
nCmb = size(cmb,1);
% Determine combinations indices to skip
skipInd = [];
for skipi = 1:length(zeroedChannel)
    skipInd = [skipInd,logicFind(zeroedChannel(skipi),cmb(:,1),'=='),logicFind(zeroedChannel(skipi),cmb(:,2),'==')]; %#ok<AGROW>
end
% Use unique subset
skipInd = unique(skipInd);
% Set colors to be used
cols = mat2cell(distinguishable_colors(size(cmb,1)),(ones(size(cmb,1),1)),3);
% Create foi vector
foiV = (foi(1):foi(2):foi(3));
cohPlots{1} = figure;
% Preallocate coh structure
coh = cell(1,size(eoi,1));
for e = 1:size(eoi,1)
    Cxy = zeros(size(cmb,1),length(foiV),size(trls{1,e}.trial,3)); 
    CxyFull = zeros(size(cmb,1),length(foiV));
    for c = 1:size(cmb,1)
        tic
        % Check if channel should be skipped
        if ismember(c,skipInd)
            disp(['Skipping combination ',num2str(c)])
            % Set Cxy to NaN
            Cxy(c,:,:) = NaN;
        else
            disp(['Calculating combination ',num2str(c)])
            for t = 1:size(trls{1,e}.trial,3)
                [Cxy(c,:,t),F1] = mscohere(...
                    trls{1,e}.trial(cmb(c,1),:,t),...
                    trls{1,e}.trial(cmb(c,2),:,t),...
                    hamming(winSize),[],foiV,adfreq);
                notchInd = [nearest_idx3(57.5,F1);...
                    nearest_idx3(62.5,F1)];
                Cxy(c,notchInd(1):notchInd(2),t) = NaN;
                Cxy(c,notchInd(1):notchInd(2),t) = interp1(...
                    find(~isnan(Cxy(c,:,t))),...
                    Cxy(c,~isnan(Cxy(c,:,t)),t),...
                    find(isnan(Cxy(c,:,t))),'linear');
            end
            % Whole data set
            % Check for completeness of signal, if NaN exists, skip full-coherence
            nanChk = sum(sum(LFPTs.data));
            if ~isnan(nanChk)
                [CxyFull(c,:),F1] = mscohere(LFPTs.data(cmb(c,1),:),LFPTs.data(cmb(c,2),:),hamming(1024),[],foiV,adfreq);
                % Interpolate over notch filter data
                notchInd = [nearest_idx3(57.5,F1);nearest_idx3(62.5,F1)];
                CxyFull(c,notchInd(1):notchInd(2)) = NaN;
                CxyFull(c,notchInd(1):notchInd(2)) = interp1(find(~isnan(CxyFull(c,:))),CxyFull(c,~isnan(CxyFull(c,:))),find(isnan(CxyFull(c,:))),'linear');
            end
        end
        toc
    end
    %% Get population statistics for all trial coherence
    mCxy = squeeze(nanmean(Cxy,3));
    sdCxy = squeeze(std(Cxy,[],3,'omitnan'));
    %% Calculate average coherence per frequency band
    % Set frequency band limits
    %bands.thet.limit = [4,7]; bands.alph.limit = [8,13]; bands.bet.limit = [15,30]; bands.lgam.limit = [45,65]; bands.hgam.limit = [70,90];
    % Find indices corresponding to frequency bands
    bInd = bandIndices(bands,F1);
    %% Get average coherence within each band per channel pair and 
    % normalize by average coherence of channel pair (above)
    bandCoh = zeros(size(cmb,1),size(bInd,1),size(trls{1,e}.trial,3));
    totalCoh = zeros(size(bandCoh));
    for iC = 1:size(cmb,1)
        for iT = 1:size(trls{1,e}.trial,3)
            for iB = 1:size(bInd,1)
                bandCoh(iC,iB,iT) = mean(Cxy(iC,bInd(iB,1):bInd(iB,2),iT));
            end
%             totalCoh(iC,:,iT) = repmat(mean(Cxy(iC,bInd(1,1):bInd(end,2),iT)),1,nCmb);
%             totalCoh(iC,:,iT) = repmat(mean(Cxy(iC,bInd(1,1):bInd(end,2),iT)),1,size(bands,1));
        end
        totalCoh(iC) = mean(mCxy(iC,bInd(1,1):bInd(end,2)));
    end
    totalCoh = repmat(mean(mCxy(:,bInd(1,1):bInd(end,2)),2),1,size(bands,1),size(trls{1,e}.trial,3));
    relCoh = bandCoh./totalCoh;
    %%
%     for iC = 1:size(cmb,1)
%         for iB = 1:size(bInd,1)
%             avgBand(iC,iB) = mean(mCxy(iC,bInd(iB,1):bInd(iB,2)));
%         end
%         avgTotal(iC,:) = repmat(mean(mCxy(iC,bInd(1,1):bInd(end,2))),1,size(bInd,1));
%     end
%     avgRelCoh = avgBand./avgTotal;
    % Get mean coherence per band and normalize by average across bands
%     bandCoh = zeros(size(bInd,1),size(cmb,1));
%     relCoh = zeros(size(bandCoh));
%     for j = 1:size(bInd,1)
%         for k = 1:size(cmb,1)
%             bandCoh(j,k) = mean(mCxy(k,bInd(j,1):bInd(j,2),:));
%             relCoh(j,k) = bandCoh(j,k)./mean(mCxy(k,bInd(1,1):bInd(end,2),:));
%         end
%     end
    %% Plot full and mean with error coherence, if full exists
    cohPlots{1} = figure;
    if ~isnan(nanChk)
        ax1 = subplot(1,size(eoi,1)+1,1); hold on;
        for c = 1:size(cmb,1)
            plot(F1,CxyFull(c,:));
        end
        title('Full Recording Coherence'); xlabel('Frequency (Hz)'); ylabel('Coherence');
        ax2 = subplot(1,size(eoi,1)+1,e+1); hold on;
        for c = 1:size(cmb,1)
            shadedErrorBar(F1,mCxy(c,:),sdCxy(c,:),{'color',cols{c,:},'linewidth',2},1);
        end
        title(['Average Trialized Coherence: ',eoi{e,1}]); xlabel('Frequency (Hz)'); ylabel('Coherence');
        linkaxes([ax1,ax2],'xy');
        % Otherwise just plot mean coherence with error
    else
        subplot(1,size(eoi,1),e)
        hold on;
        for c = 1:size(cmb,1)
            shadedErrorBar(F1,mCxy(c,:),sdCxy(c,:),{'color',cols{c,:},'linewidth',2},1);
        end
        title(['Average Trialized Coherence: ',eoi{e,1}]); xlabel('Frequency (Hz)'); ylabel('Coherence');
    end
    %% Setup output structure 'coh'
    coh{1,e}.Cxy = Cxy;
    if ~isnan(nanChk)
        coh{1,e}.CxyFull = CxyFull;
    end
    coh{1,e}.mCxy = mCxy;
    coh{1,e}.sdCxy = sdCxy;
    coh{1,e}.rel = relCoh;
    coh{1,e}.band = bandCoh;
    coh{1,e}.freq = F1;
    coh{1,e}.winSize = 1024;
    coh{1,e}.total = totalCoh;
%     coh{1,e}.avgRel = avgRelCoh;
end
