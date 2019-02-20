%does stimulation on channel 2 (stimulator 1) using pattern from file
%'3_pulse_burst_fixed.pat'

res = PS_InitAllStim
NStim = PS_GetNStim

res = PS_SetMonitorChannel(1, 2)
PS_SetRepetitions(1, 2, 0)

PS_SetPatternType(1, 2, 1);
res = PS_LoadArbPattern(1, 2, 'c:\PlexonData\Stim-2\Waveform pattern files\3_pulse_burst_fixed.pat')

%get points to draw pattern
NPoints = PS_GetNPointsArbPattern(1, 2)
PointsArrayX = PS_GetArbPatternPointsX(1, 2)
PointsArrayY = PS_GetArbPatternPointsY(1, 2)
plot(PointsArrayX, PointsArrayY);

res = PS_LoadChannel(1,2)
res = PS_StartStimChannel(1,2);

input('Press Return to Stop');
res = PS_StopStimChannel(1,2);
res = PS_CloseStim(1)
