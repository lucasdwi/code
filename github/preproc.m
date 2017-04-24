function [output] = preproc(LFPTs,cfg)
%cfg
% 1. Line noise filter
% 2. NaNs noise artifacts based on threshold
% 3. Smooth data for detection of sleep events with SleepGUI

chans = size(LFPTs.data,1);
[LFPTs] = lineFilt(LFPTs);
[LFPTs,~] = threshFilt(LFPTs,thresh,onset,offset,minInt,adfreq,1,chans);
[smInstAmp,LFPTs] = sleepDetect(LFPTs,minInt,thresh,onset,offset,adfreq,smMthod);

