function [coh,cohPlots] = cohComp(trls,adfreq,eoi,bands,zeroedChannel,foi,vis,method,varargin)
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
% vis = whether or not to visualize results; format = 0 (n) or 1 (y)
% foi = frequencies of interest; format = [low frequency, step, high
%   frequency]
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
addRequired(p,'vis');
addRequired(p,'method');
% Needed for 'mtm' method
addParameter(p,'NW',8);
% Needed for 'mat' method
addParameter(p,'overlap',0.5);
parse(p,trls,adfreq,eoi,bands,zeroedChannel,vis,method,varargin{:});
%%
% Count channels and check that there are the same number of channels in
% all trls
chans = unique(cellfun(@(x) size(x.trial,1),trls));
if length(chans) > 1 
    error('There is a inequal number of channels in trls.')
end
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
    overlapSamp = round(overlap*winSize);
end
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
            % Set up the two signals
            x = trls{1,ei}.trial(cmbs(ci,1),:,ti);
            y = trls{1,ei}.trial(cmbs(ci,2),:,ti);
            % Calculate coherence
            if strcmpi(method,'mat')
                [Cxy(ci,:,ti),f] = mscohere(x,y,window,overlapSamp,foiV,...
                    adfreq);
            else if strcmpi(method,'mtm')
                    [f,thisCxy,~,~,~] = cmtm(x,y,1/adfreq,p.Results.NW,...
                        0,0,0);
                    % Truncate frequeny vector to end of foi
                    f = f(1:nearest_idx3(foi(end),f));
                    % Truncate signal and store
                    Cxy(ci,:,ti) = thisCxy(1:nearest_idx3(foi(end),f));
                end
            end
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
        clear h
        leg = zeros(1,nCmbs);
        subplot(1,size(eoi,1),ei)
        hold on
        % Plot mean and 95% confidence interval (if 'mat' method) of
        % coherence for each pair
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
            legLabel{ci} = [num2str(cmbs(ci,1)),'-',num2str(cmbs(ci,2))];
        end
        % Populate title and axes
        title([eoi{ei,1},' - average'])
        xlabel('Frequency (Hz)')
        ylabel('Coherence')
        % Add legend
        legend(leg,legLabel,'Orientation','horizontal','Location',...
            'southoutside')
    end
    % Set up second figure - normalized band coherence across trials
    cohPlots{2} = figure;
    % Plot each event as it's own subplot
    for ei = 1:size(eoi,1)
        clear h
        subplot(1,size(eoi,1),ei)
        hold on
        for ci = 1:nCmbs
            h{ci} = scatterErr(1:6,mean(coh{ei}.normBandCoh(ci,:,:),3),...
                conf(squeeze(coh{ei}.normBandCoh(ci,:,:)),0.95),0,...
                'col',cols{ci,:});
        end
        % Populate title and axes
        title([eoi{ei,1},' - normalized'])
        set(gca,'XTickLabel',bands(:,1))
        ylabel('Ratio of Average Coherence')
        % Add legend
        legend(h,bands(:,1))
        clear h
    end
end