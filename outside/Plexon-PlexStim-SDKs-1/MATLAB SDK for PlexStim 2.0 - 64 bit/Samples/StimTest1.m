%does stimulation on channel 1 (stimulator 1) using predefined rectangular pulse

err = PS_InitAllStim
[NStim, err] = PS_GetNStim
[NChan, err] = PS_GetNChannels(1)

err = PS_SetMonitorChannel(1, 1)
err = PS_SetRepetitions(1, 1, 0)
err = PS_LoadChannel(1,1)
err = PS_StartStimChannel(1,1);

input('Press Return to Stop');
err = PS_StopStimChannel(1,1);
err = PS_CloseStim(1)
