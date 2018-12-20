function [coh,cohPlots] = cohComp(trls,adfreq,eoi,bands,zeroedChannel,foi,nFilt,vis,method,varargin)
%% Calculates coherence using either a single hamming window and mscohere
% or the multitaper method proposed by Percival and Walden.
%__________________________________________________________________________
% INPUTS
% trls = structure containing data to be analyzed; requires a field
%   'trials' which contains raw data channel X time X trial
% adfreq = sampling rate; format = integer
% eoi = events of interest; format = cell array, first column has name of
%   event and second column has timing around event
% bands = cell array of frequency bands; format = cell array, first column
%   has name of bands and second column has inclusive start and stop of
%   bands
% zeroedChannel = index of any channels to be skipped due to noise; format
%   = integer
% foi = frequencies of interest; format = [low frequency, step, high
%   frequency]
% nFilt = frequencies to interpolate over to account for notch filter; if
%   no notch filter was applied, leave empty; format = [low high]
% vis = whether or not to visualize results; format = 0 (n) or 1 (y)
% method = which method to use for calculating coherence; format = either
%   'mat' for mscohere or 'mtm' for multitaper method
% NW = number of windows to use for multitaper method; format = integer,
%   default is 8
% overlap = percent overlap between windows for mscohere; format = decimal,
%   default is 0.5 (50%)
%__________________________________________________________________________
% OUTPUTS
% coh = data structure with one cell per event with the following fields
%   Cxy = coherence; format = channel pair X frequency X trial
%   mtCxy = mean coherence over trials; format = channel pair X band X
%       trial (replicated values across bands)
%   mBandCoh = mean coherence per band; format = channel pair X band X
%       trial
%   normBandCoh = normalized band coherence; format = channel pair X band X
%       trial (mBandCoh./mtCxy)
%   f = frequency vector corresponding to 2nd dimension of Cxy; format = Hz
%__________________________________________________________________________
% LLD 2017-2018
% Uses cmtm.m from Peter Huybers to implement multitaper method
%% Parse inputs
p = inputParser;
p.CaseSensitive = false;
addRequired(p,'trls');
addRequired(p,'adfreq',@isnumeric);
addRequired(p,'eoi');
addRequired(p,'bands');
addRequired(p,'zeroedChannel');
addRequired(p,'foi');
addRequired(p,'filter');
addRequired(p,'vis');
addRequired(p,'method');
% Needed for 'mtm' method
addParameter(p,'NW',8);
% Needed for 'mat' method
addParameter(p,'overlap',0.5);
parse(p,trls,adfreq,eoi,bands,zeroedChannel,foi,nFilt,vis,method,...
    varargin{:});
%%
% Count channels using first non-empty trls struct
empt = cellfun(@isempty,trls);
chans = size(trls{logicFind(0,empt,'==','first')}.trial,1);
% Get channel pairs for coherence
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
if strcmpi(method,'mat')
    % Create foi vector
    foiV = foi(1):foi(2):foi(3);
    % Calculate window size by finding nearest power of two for 1 second
    [~,winSize] = nearestPow2(adfreq);
    % Set up window
    window = hamming(winSize);
    % Calculate overlap samples
    oSamp = round(p.Results.overlap*winSize);
end
% Preallocate coh
coh = cell(1,size(eoi,1));
% Go through non-empty trls
for ei = logicFind(0,empt,'==')
    % Get number of trials
    nTrls = size(trls{1,ei}.trial,3);
    % Preallocate mBandCoh
    mBandCoh = zeros(nCmbs,nBands,nTrls);
    % Preallocate Cxy depending on method being used
    if strcmpi(method,'mat')
        Cxy = zeros(nCmbs,length(foiV),nTrls);
    % If using mtm, then use formula used to calculate frequency range in
    % cmtm.m
    elseif strcmpi(method,'mtm')
        n = size(trls{1,ei}.trial,2);
        dt = 1/adfreq;
        s = 0:1/(n*dt):1/dt-1/(n*dt);
        pls = 2:(n+1)/2+1;
        s = s(pls);
        Cxy = zeros(nCmbs,nearest_idx3(foi(end),s),nTrls);
    end
    
    for ci = 1:nCmbs
        disp(['Event ',num2str(ei),' of ',num2str(size(eoi,1)),...
            ': Calculating coherence for channel pair ',...
            num2str(cmbs(ci,1)),'-',num2str(cmbs(ci,2))])
        for ti = 1:nTrls
            % Set up the two signals
            x = trls{1,ei}.trial(cmbs(ci,1),:,ti);
            y = trls{1,ei}.trial(cmbs(ci,2),:,ti);
            % Calculate coherence
            if strcmpi(method,'mat')
                [Cxy(ci,:,ti),f] = mscohere(x,y,window,oSamp,foiV,adfreq);
            elseif strcmpi(method,'mtm')
                [f,thisCxy,~,~,~] = cmtm(x,y,1/adfreq,p.Results.NW,0,0,0);
                % Truncate frequeny vector to end of foi
                f = f(1:nearest_idx3(foi(end),f));
                % Truncate signal and store
                Cxy(ci,:,ti) = thisCxy(1:nearest_idx3(foi(end),f));
            end
            % Check if notch filter needs to be interpolated over
            if ~isempty(nFilt)
                % Set up notch for 60 Hz line noise
                notchInd = [nearest_idx3(nFilt(1),f);...
                    nearest_idx3(nFilt(2),f)];
                % Interpolate over notch
                samp = [notchInd(1),notchInd(2)];
                v = [Cxy(ci,notchInd(1),ti),Cxy(ci,notchInd(2),ti)];
                sampQ = notchInd(1):notchInd(2);
                Cxy(ci,sampQ,ti) = interp1(samp,v,sampQ,'linear');
            end
        end
    end
    % Get indices for frequency bands
    bandInd = bandIndices(bands,f);
    % Get mean coherence across frequencies per trial per channel pair;
    % replicate to match dimension of 'mBandCoh'
%     mtCxy = repmat(mean(Cxy,2),1,size(bands,1),1);
    % Only normalize by frequencies in range of analysis
    mtCxy = repmat(mean(Cxy(:,bandInd(1,1):bandInd(end,2),:),2),1,...
        size(bands,1),1);
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
if strcmpi(vis,'y')
    % Set colors to be used
    cols = mat2cell(distinguishable_colors(size(cmbs,1)),...
        ones(size(cmbs,1),1),3);
    % Set up first figure - mean coherence across trials with 95% conf
    cohPlots{1} = figure;
    set(gcf,'Position',[930,130,710,560]);
    % Plot each event as it's own subplot
    for ei = 1:size(eoi,1)
        clear h
        leg = zeros(1,nCmbs);
        subplot(1,size(eoi,1),ei)
        hold on
        % Plot mean and 95% confidence interval (if 'mat' method) of
        % coherence for each pair - only plot if data exist
        if ~empt(ei)
            for ci = 1:nCmbs
                if strcmpi(method,'mat')
                    h(ci) = shadedErrorBar(coh{ei}.f,squeeze(mean(...
                        coh{ei}.Cxy(ci,:,:),3)),conf(squeeze(...
                        coh{ei}.Cxy(ci,:,:)),0.95),{'color',cols{ci,:}},1);
                    % Grab legend info
                    leg(ci) = h(ci).mainLine;
                else
                    h(ci) = plot(coh{ei}.f,squeeze(mean(...
                        coh{ei}.Cxy(ci,:,:),3)),'color',cols{ci,:});
                    % Grab legend info
                    leg(ci) = h(ci);
                end
                legLabel{ci} = [num2str(cmbs(ci,1)),'-',...
                    num2str(cmbs(ci,2))];
            end
            % Add legend on outside if last subplot
            if ei == size(eoi,1)
                legend(leg,legLabel,'Location','eastoutside')
                set(legend,'Position',[0.25 0 0.5 0.04],...
                    'Orientation','horizontal');
            end
        else
            text(50,0.5,'No trials exist!','HorizontalAlignment','center')
        end
        % Populate title
        title([eoi{ei,1},': average'])
        % Set x axis limits
        xlim([0 100])
        % Label axes
        if ei == 1
            ylabel('Coherence')
        end
        xlabel('Frequency (Hz)')
    end
    % Set up second figure - normalized band coherence across trials
    cohPlots{2} = figure;
    % Plot each event as it's own subplot
    for ei = 1:size(eoi,1)
        clear h
        subplot(1,size(eoi,1),ei)
        hold on
        if ~empt(ei)
            for ci = 1:nCmbs
                % Check for only one trial
                if size(coh{ei}.normBandCoh,3) == 1
                    h(ci) = plot(1:6,coh{ei}.normBandCoh(ci,:),'-o',...
                        'col',cols{ci,:},'LineWidth',1);
                else
                    h(ci) = plot(1:6,mean(coh{ei}.normBandCoh(ci,:,:),3)...
                        ,'-o','col',cols{ci,:},'LineWidth',1);
                    scatterErr(1:6,mean(coh{ei}.normBandCoh(ci,:,:),3),...
                        conf(squeeze(coh{ei}.normBandCoh(ci,:,:)),0.95),0,...
                        'col',cols{ci,:});
                end
            end
        else
            text((size(bands,1)+0.5)/2,0.5,'No trials exist!',...
                'HorizontalAlignment','center')
        end
        % Populate title and axes
        title([eoi{ei,1},': norm.'])
        set(gca,'Xtick',1:size(bands,1),'XTickLabel',latinToGreek(bands(:,1)))
        ylabel('Ratio of Average Coherence')
        xlim([0.5 size(bands,1)+0.5])
        % Add legend
        if ei == size(eoi,1)
            legend(h,legLabel)
        end
    end
else
    cohPlots = [];
end