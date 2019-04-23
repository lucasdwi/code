%does stimulation on channel 1 (stimulator 1) using rectangular pulse with
%new parameters

res = PS_InitAllStim
[NStim, err] = PS_GetNStim
[NChan, err] = PS_GetNChannels(1)

res = PS_SetMonitorChannel(1, 1)
PS_SetRepetitions(1, 1, 0)

%check parameters for the rectangular pulse
[Values, err] = PS_GetRectParam(1, 1)

%set new parameters for the rectangular pulse
pattern = [100, -100, 25, 25, 25];
PS_SetRectParam(1, 1, pattern);

%check if new parameters for the rectangular pulse has been set up
Values = PS_GetRectParam(1, 1)


err = PS_LoadChannel(1,1)
err = PS_StartStimChannel(1,1);

input('Press Return to Stop');
err = PS_StopStimChannel(1,1);
err = PS_CloseStim(1)
