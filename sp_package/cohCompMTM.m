function [coh] = cohCompMTM(trls,adfreq,NW,eoi,bands,chans,zeroedChannel,winSize,overlap,vis)
% Get channel combinations
cmbs = nchoosek(1:chans,2);
% Determine combinations to skip due to zeroed channels
if ~isempty(zeroedChannel)
    skipInd = logicFind(1,any(cmbs == zeroedChannel),2,'==');
    % Remove skipInds from cmbs
    cmbs(skipInd,:) = [];
end
% Count number of combinations and bands
nCmbs = size(cmbs,1);
nBands = size(bands,1);
% Preallocate output structure
coh = cell(1,size(eoi,1));
for ei = 1:size(eoi,1)
    % Get number of trials
    nTrls = size(trls{1,ei}.trial,3);
    % Preallocate 'Cxy' and 'mBandCoh'
%     Cxy = zeros(nCmbs,length(foiV),nTrls);
    mBandCoh = zeros(nCmbs,nBands,nTrls);
    for ci = 1:nCmbs
        disp(['Event ',num2str(ei),' of ',num2str(size(eoi,1)),...
            ': Calculating coherence for channel pair ',...
            num2str(cmbs(ci,1)),'-',num2str(cmbs(ci,2))])
        for ti = 1:nTrls
            if discrete
            % Set up the two signals
            x = trls{1,ei}.trial(cmbs(ci,1),:,ti);
            y = trls{1,ei}.trial(cmbs(ci,2),:,ti);
            % Calculate coherence
            [f,thisCxy,~,~,~] = cmtm(x,y,1/adfreq,NW,0,0,0);
            else
                % Convert winSize from seconds to samples
                winSize = winSize*adfreq;
                % Set up counter and start
                c = 0;
                start = 1;
                while start+winSize < size(trls{1,ei}.trial,2)
                    start = 1+c*winSize*overlap;
                    stop = start+winSize;
                    x = trls{1,ei}.trial(cmbs(ci,1),start:stop,ti);
                    y = trls{1,ei}.trial(cmbs(ci,2),start:stop,ti);
                    [f,thisCxy(:,c)] = cmtm(x,y,1/adfreq,NW,0,0,0);
                    c = c+1;
                end
            end
            % Truncate frequeny vector to 100 Hz
            f = f(1:nearest_idx3(100,f));
            % Truncate signal and store
            Cxy(ci,:,ti) = thisCxy(1:nearest_idx3(100,f));
            % Set up notch for 60 Hz line noise
            notchInd = [nearest_idx3(57.5,f);nearest_idx3(62.5,f)];
             % Interpolate over notch
            samp = [notchInd(1)-1,notchInd(2)+1];
            v = [Cxy(ci,notchInd(1)-1,ti),Cxy(ci,notchInd(2)+1,ti)];
            sampQ = notchInd(1):notchInd(2);
            Cxy(ci,notchInd(1):notchInd(2),ti) = interp1(samp,v,sampQ);
        end
    end
    % Get mean coherence across frequencies per trial per channel pair;
    % replicate to match dimension of 'mBandCoh'
    mtCxy = repmat(mean(Cxy,2),1,size(bands,1),1);
    % Get indices for frequency bands
    bandInd = bandIndices(bands,f);
    % Get mean coherence of frequency band per trial per channel pair
    for bi = 1:nBands
       mBandCoh(:,bi,:) = mean(Cxy(:,bandInd(bi,1):bandInd(bi,2),:),2);
    end
    % Normalize 'mBandCoh' by 'mtCxy'
    normBandCoh = mBandCoh./mtCxy;
    % Setup output structure 'coh'
    coh{ei}.Cxy = Cxy;
    coh{ei}.mtCxy = mtCxy;
    coh{ei}.mBandCoh = mBandCoh;
    coh{ei}.normBandCoh = normBandCoh;
    coh{ei}.f = f;
end
%% Plot
if isequal(vis,1)
    % Set colors to be used
    cols = mat2cell(distinguishable_colors(size(cmbs,1)),...
        ones(size(cmbs,1),1),3);
    % Set up first figure - mean coherence across trials with 95% conf
    cohPlots{1} = figure;
    % Plot each event as it's own subplot
    for ei = 1:size(eoi,1)
        h = zeros(1,nCmbs);
        leg = zeros(1,nCmbs);
        subplot(1,size(eoi,1),ei)
        hold on
        % Plot mean and 95% confidence interval of coherence for each pair
        for ci = 1:nCmbs
            h(ci) = shadedErrorBar(coh{ei}.f,squeeze(mean(coh{ei}.Cxy(ci,:,:),3)),...
                conf(squeeze(coh{ei}.Cxy(ci,:,:)),0.95),{'color',cols{ci,:}},1);
            % Grab legend info
            leg(ci) = h(ci).mainLine;
        end
        % Populate title and axes
        title([eoi{ei,1},' - average'])
        xlabel('Frequency (Hz)')
        ylabel('Coherence')
        % Add legend
        legend(leg,bands(:,1),'Orientation','horizontal','Location',...
            'southoutside')
        % Clear out h
        clear h
    end
    % Set up second figure - normalized band coherence across trials
    cohPlots{2} = figure;
    % Plot each event as it's own subplot
    for ei = 1:size(eoi,1)
        subplot(1,size(eoi,1),ei)
        hold on
        for ci = 1:nCmbs
            h(ci) = scatterErr(1:6,mean(coh{ei}.normBandCoh(ci,:,:),3),...
                conf(squeeze(coh{ei}.normBandCoh{ei}(ci,:,:)),0.95),0,...
                'col',cols{ci,:});
        end
        % Populate title and axes
        title([eoi{ei,1},' - normalized'])
        set(gca,'XTickLabel',bands(:,1))
        ylabel('Ratio of Average Coherence')
        % Add legend
        legend(h,bands(:,1))
    end
end