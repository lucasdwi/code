function trial = ddtTrials(eventTs,LFPTs,vis)
%% NB: two forced trials (of any kind) delineate start
%%
% Count number of tones
nTones = size(eventTs.t{7},1);
% Replace eventTs labels
eventTs.label = {'levIM','levDEL','head','feeder','ltDEL','ltIM','tone',...
    'TrialInit'};
% Preallocate structure
trial = struct('outcome',NaN,'start',NaN,'stop',NaN,'blockN',NaN,...
    'trialN',NaN,'leverPress',NaN,'feeder',NaN,'headEntry',NaN,...
    'tone_nosePoke_latency',NaN,'nosePoke_lever_latency',NaN,...
    'lever_headEntry_latency',NaN,'blockPercentDelay',NaN);
trial = repmat(trial,nTones,1);
% Set up counters
forcedCounter = 0;
blockCounter = 1;
trialCounter = 1;
% Set skipTone to false
skipTone = 0;
% Set tone to NaN so that first trial doesn't fail
tone = NaN;
% Use tone as start time for trial
for ii = 1:nTones
    % Check if previous tone is the same as this one
    if tone == eventTs.t{7}(ii)
        skipTone = 1;
    end
    if skipTone == 0
        % Tone
        tone = eventTs.t{7}(ii);
        % Next tone
        if ii ~= nTones
            nextTone = eventTs.t{7}(ii+1);
        else
            nextTone = NaN;
        end
        % Check for double-tone
        if nextTone - tone < 30
            % Skip next tone
            skipTone = 1;
        end
        % Add blockN
        trial(ii).blockN = blockCounter;
        % Find nearest (subsequent) nose poke unless nose poke and tone are
        % the same
        if eventTs.t{8}(nearest_idx3(tone,eventTs.t{8})) == tone
            nearestNosePoke = eventTs.t{8}(nearest_idx3(tone,...
                eventTs.t{8}));
        else
            nearestNosePoke = eventTs.t{8}(nearest_idx3(tone,...
                eventTs.t{8},1));
        end
        % If nearestNosePoke occurs before tone, set to NaN
        if nearestNosePoke < tone
            nearestNosePoke = NaN;
        end
        % Determine if nose poke occurs within 49 seconds of tone or at all
        if nearestNosePoke >= tone+49 || isnan(nearestNosePoke)
            % Determine if timeout was during forced or free trial
            if forcedCounter < 2
                trial(ii).outcome = 'timeout_forced';
            else
                trial(ii).outcome = 'timeout_free';
                % Set trialN
                trial(ii).trialN = trialCounter;
                % Add one to trial counter
                trialCounter = trialCounter + 1;
            end
        else
            % Calculate latency to nose poke after tone
            trial(ii).tone_nosePoke_latency = nearestNosePoke - tone;
            % Get nearest feeder
            thisFeed = eventTs.t{4}(nearest_idx3(tone,eventTs.t{4},1));
            % Check that thisFeed occurs AFTER tone
            if thisFeed < tone
               thisFeed = NaN; 
            end
            % Check if last feed, if not get next feeder
            if nearest_idx3(tone,eventTs.t{4},1) == numel(eventTs.t{4})...
                    || ii == nTones
                nextFeed = NaN;
            else
                nextFeed = eventTs.t{4}(nearest_idx3(tone,eventTs.t{4},...
                    1)+1);
            end
            % Determine if forced or free
            if forcedCounter < 2
                % Check for failed forced trial (i.e., no feeder before
                % next tone)
                if thisFeed >= nextTone
                    trial(ii).outcome = 'failed_forced';
                else
                    thisReward = whichLever(eventTs,tone);
                    trial(ii).outcome = [thisReward,'_forced'];
                    % Add one to forcedCounter
                    forcedCounter = forcedCounter + 1;
                    % Add feeder
                    trial(ii).feeder = thisFeed;
                    % Grab lever press time, using thisReward to determine
                    % lever
                    trial = lpt(thisReward,eventTs,trial,ii,tone);
                    % Get head entry time post lever press
                    trial(ii).headEntry = eventTs.t{3}(nearest_idx3(...
                        trial(ii).leverPress,eventTs.t{3},1));
                    % Calculate latency to lever press from nose poke
                    trial = npl(trial,ii,nearestNosePoke);
                    % Calculate latency to head entry from lever press
                    trial = lhe(trial,ii);
                end
                % Free
            else
                % Determine if any lever is pressed (i.e., did feeder
                % activate after this tone or after next tone and the next
                % tone isn't a skip)
                if isnan(thisFeed) || thisFeed >= nextTone && ~skipTone 
                    trial(ii).outcome = 'failed_free';
                else
                    % Otherwise, success
                    thisReward = whichLever(eventTs,tone);
                    trial(ii).outcome = [thisReward,'_free'];
                    trial(ii).feeder = thisFeed;
                    % Grab lever press time, using thisReward to determine
                    % lever
                    trial = lpt(thisReward,eventTs,trial,ii,tone);
                    % Calculate latency to lever press from nose poke
                    trial = npl(trial,ii,nearestNosePoke);
                    % Get head entry time post lever press
                    trial(ii).headEntry = eventTs.t{3}(nearest_idx3(...
                        trial(ii).leverPress,eventTs.t{3},1));
                    % Calculate latency to head entry from lever press
                    trial = lhe(trial,ii);
                end
                % Set trialN
                trial(ii).trialN = trialCounter;
                % Add one to trial counter
                trialCounter = trialCounter + 1;
            end
        end
        % Get start and stop time for trial
        trial(ii).start = LFPTs.tvec(nearest_idx3(tone,LFPTs.tvec));
        % Use next tone to find stop time, unless last tone
        if ii == nTones
            trial(ii).stop = LFPTs.tvec(end);
        else
            if skipTone == 1
                % Check if skipped tone is last
                if ii+1 == nTones
                    trial(ii).stop = LFPTs.tvec(end);
                else
                    trial(ii).stop = LFPTs.tvec(nearest_idx3(...
                        eventTs.t{7}(ii+2),LFPTs.tvec)-1);
                end
            else
                trial(ii).stop = LFPTs.tvec(nearest_idx3(...
                    eventTs.t{7}(ii+1),LFPTs.tvec)-1);
            end
        end
        % If a block is finished (10 trials; trialCoutner = 11), move to
        % next block, reset trialCounter, reset forcedCounter
        if trialCounter == 11
            blockCounter = blockCounter + 1;
            trialCounter = 1;
            forcedCounter = 0;
        end
    else
        % Reset skipTone to false
        skipTone = 0;
        % Set outcome to 'doubleTone'
        trial(ii).outcome = 'doubleTone';
        % Set tone
        tone = eventTs.t{7}(ii);
    end
end
% Determine block statistics - percent delay choices per block
% Convert structure to cell
trialCell = struct2cell(trial);
% Convert blockN cells to array
blocks = cell2mat(trialCell(4,:));
% Calculate delays from lever press to feeder
delays = round(extractfield(trial,'feeder')-...
    extractfield(trial,'leverPress'),1);
% Cycle through each block to add percent delay choice and block delay
for ii = 1:5
    % Percent delay choice
    % Find indices of this block
    inds = blocks == ii & ~isnan(cell2mat(trialCell(5,:)));
    % Calculate percent delay choices
    percentDelay(ii) = sum(strcmp('delay_free',trialCell(1,inds)))/...
        (sum(strcmp('delay_free',trialCell(1,inds)))+...
        sum(strcmp('immediate_free',trialCell(1,inds))));
    % Input percent delay into each trial of this block
    for jj = logicFind(1,inds,'==')
        trial(jj).blockPercentDelay = percentDelay(ii);
    end
    % Delay time
    delayInds = blocks == ii & strcmp(trialCell(1,:),'delay_forced');
    theseDelays(ii) = delays(delayInds);
    % Input delay into each trial of this block
    for jj = logicFind(1,blocks == ii,'==')
        trial(jj).delay = theseDelays(ii);
    end
end
[trial(:).auc] = deal(sum(percentDelay)*100);
[trial(:).trapzAUC] = deal(trapz(theseDelays,percentDelay));
% Visualization
if vis
    figure
    hold on
    for ii = 1:8
        plot(eventTs.t{ii},ones(1,numel(eventTs.t{ii}))*ii,'.')
    end
    % Get block numbers
    uBlocks = unique(blocks(~isnan(blocks)));
    % Find last trial of each block
    bInds = zeros(numel(uBlocks),1);
    for ii = 1:numel(uBlocks)
        bInds(ii) = logicFind(uBlocks(ii),blocks,'==','last');
        plot([trial(bInds(ii)).stop-2 trial(bInds(ii)).stop-2],...
            [0.5 8.5],'--k')
    end
    for ii = 1:numel(eventTs.t{7})
        text(eventTs.t{7}(ii),7.5,num2str(ii),'HorizontalAlignment',...
            'center')
    end
    set(gca,'yticklabel',eventTs.label(1:8),'ytick',(1:8))
    ylim([0.5 8.5])
    xlabel('Time (s)')
end
end
%%
function thisReward = whichLever(eventTs,tone)
% Get next lever presses after tone
nextImmediate = eventTs.t{1}(nearest_idx3(tone,eventTs.t{1},1));
nextDelay = eventTs.t{2}(nearest_idx3(tone,eventTs.t{2},1));
% Make sure nextImmediate and nextDelay occur after tone
if nextImmediate < tone
    % If not, set to Inf so that it will not be smallest
    nextImmediate = Inf;
end
if nextDelay < tone
    nextDelay = Inf;
end
if nextDelay < nextImmediate
    thisReward = 'delay';
else
    thisReward = 'immediate';
end
end
%%
function trial = lpt(thisReward,eventTs,trial,ii,tone)
% Determines Lever Press Time; i.e., when lever was pressed; either
% 'delay' or 'immediate'
if strcmp(thisReward,'delay')
    trial(ii).leverPress = eventTs.t{2}(...
        nearest_idx3(tone,eventTs.t{2},1));
else
    trial(ii).leverPress = eventTs.t{1}(...
        nearest_idx3(tone,eventTs.t{1},1));
end
end
%%
function trial = npl(trial,ii,nearestNosePoke)
% Determine Nose Poke to Lever press latency
trial(ii).nosePoke_lever_latency = trial(ii).leverPress - nearestNosePoke;
end
%%
function trial = lhe(trial,ii)
% Determine Lever to Head Entry latency
trial(ii).lever_headEntry_latency = trial(ii).headEntry - ...
    trial(ii).leverPress;
end