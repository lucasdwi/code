chans = 1;
f = 130;
rep = 0;
params = [200, -200, 90, 90, 15];
baseTime = 600;
firstTime = 600;
stimOn = 55;
stimOff = 5;
totStim = 3600;
nStim = (stimOn+stimOff)/totStim;

if size(chan,2) ~= size(f,2)
    error(['There is an imbalance between the number of channels being'...
        'stimulated and the number of frequencies given.'])
end
% Figure out mapping between 'channels' and stim channels
% Initialize stim
PS_InitAllStim
% Get number of available stimulators
stimN = PS_GetNStim;
% Get maximum number of channels for stimN
% [NChan, err] = PS_GetNChannels(7);
% [NChan, err] = PS_GetNChannels(8);

% Monitor channel 7 from stimulator 1 
PS_SetMonitorChannel(1,7)
% Setup each channel for stim
for ii = 1:size(chan,2)
    % Set frequency of stim
    PS_SetRate(1,chan(ii),f(ii))
    % Set repeititions
    PS_SetRepetitions(1,chan(ii),rep(ii))
    % Set pattern for biphasic rectangular pulse
    PS_SetPatternType(1,chan(ii),0)
    % Set parameters of biphasic rectnagular pulse
    PS_SetRectParams(1,chan(ii),params(ii,:))
    % Load channel
    PS_LoadChannel(1,chan(ii))
end
disp(['Setup complete; waiting ',num2str(baseTime),...
    ' seconds, before stimulating.'])
% Wait until baseTime has passed
pause(baseTime)
% Turn on stim for first stim chunk
PS_StartStimChannel(1,chan);
% for ii = 1:size(chan,2)
%     PS_StartStimChannel(1,chan(ii));
% end
pause(firstTime)
% Turn off for one stimOff interval
% for ii = 1:size(chan,2)
%    PS_StopStimChannel(1,chan(ii)); 
% end
PS_StopStimChannel(1,chan);
disp('First stimulation chunk complete; moving on to on/off protocol.')
% Start stim on/off protocol
for ii = 1:nStim
    PS_StartStimChannel(1,chan)
    pause(stimOn)
    PS_StopStimChannel(1,chan)
    pause(stimOff)
end
PS_CloseAllStim
