function [coh,cohPlots] = cohCompMat(LFPTs,chans,trls,foi,winSize,overlap,adfreq)
%% Cycle through trials and calculate coherence
cmb = nchoosek(1:chans,2);
% Create foi vector
foiV = (foi(1):foi(2):foi(3));
Cxy = zeros(size(cmb,1),length(foiV),length(trls{1,1}.trial)); CxyFull = zeros(size(cmb,1),length(foiV));
for c = 1:size(cmb,1)
    tic
    disp(['Calculating combination ',num2str(c)])
    for t = 1:length(trls{1}.trial)
        [Cxy(c,:,t),F1] = mscohere(trls{1}.trial{t}(cmb(c,1),:),trls{1}.trial{t}(cmb(c,2),:),hanning(winSize),winSize*overlap,foiV,adfreq);
        % Interpolate over notch filter data
        notchInd = [nearest_idx3(57.5,F1);nearest_idx3(62.5,F1)];
        Cxy(c,notchInd(1):notchInd(2),t) = NaN;
        Cxy(c,notchInd(1):notchInd(2),t) = interp1(find(~isnan(Cxy(c,:,t))),Cxy(c,~isnan(Cxy(c,:,t)),t),find(isnan(Cxy(c,:,t))),'linear');
    end
    % Whole data set
    % Check for completness of signal, if NaN exists, skip full-coherence
    nanChk = sum(sum(LFPTs.data));
    if ~isnan(nanChk)
        [CxyFull(c,:),F1] = mscohere(LFPTs.data(cmb(c,1),:),LFPTs.data(cmb(c,2),:),hanning(winSize),winSize*overlap,foiV,adfreq);
        % Interpolate over notch filter data
        notchInd = [nearest_idx3(57.5,F1);nearest_idx3(62.5,F1)];
        CxyFull(c,notchInd(1):notchInd(2)) = NaN;
        CxyFull(c,notchInd(1):notchInd(2)) = interp1(find(~isnan(CxyFull(c,:))),CxyFull(c,~isnan(CxyFull(c,:))),find(isnan(CxyFull(c,:))),'linear');
    end
    toc
end
%% Get population statistics for all trial coherence
mCxy = sq(mean(Cxy,3));
sdCxy = sq(std(Cxy,0,3));
%% Plot full and mean with error coherence, if full exists
    cols = {[0 1 1];[1 0 0];[0 1 0]; [0 0 1]; [0 0 0]; [1 0 1]};
    cohPlots{1} = figure;
if ~isnan(nanChk)
    ax1 = subplot(1,2,1); hold on;
    for c = 1:size(cmb,1)
        plot(F1,CxyFull(c,:),'color',cols{c});
    end
    title('Full Recording Coherence'); xlabel('Frequency (Hz)'); ylabel('Coherence');
    ax2 = subplot(1,2,2); hold on;
    for c = 1:size(cmb,1)
        shadedErrorBar(F1,mCxy(c,:),sdCxy(c,:),{'color',cols{c}},1);
    end
    title('Average Trialized Coherence'); xlabel('Frequency (Hz)'); ylabel('Coherence');
    linkaxes([ax1,ax2],'xy');
% Otherwise just plot mean coherence with error
else
    hold on;
    for c = 1:size(cmb,1)
        shadedErrorBar(F1,mCxy(c,:),sdCxy(c,:),{'color',cols{c}},1);
    end
    title('Average Trialized Coherence'); xlabel('Frequency (Hz)'); ylabel('Coherence');
end

%% Setup output structure 'coh'
coh.Cxy = Cxy;
if ~isnan(nanChk)
    coh.CxyFull = CxyFull;
end
coh.mCxy = mCxy;
coh.sdCxy = sdCxy;
