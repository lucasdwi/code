function [histBin] = histBuilder(x,y,nBins)

%License:

%  This software is distributed under the "Creative Commons Attribution
%  Noncommercial-Share Alike License"
%
%     Version 3.0, available at
%
%         http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
%     You are free:
%
%         To Share -- To copy, distribute and transmit the work
%
%         To Remix -- To adapt the work
%
%     Under the following conditions:
%
%         Attribution -- You must attribute the work in the manner specified
%            by the author or licensor (but not in any way that suggests that
%            they endorse you or your use of the work).
%
%         Noncommercial -- You may not use this work for commercial purposes.
%
%         Share Alike -- If you alter, transform, or build upon this work,
%            you may distribute the resulting work only under the same or
%            similar license to this one.
%
%     See the above link for the full text of the license.
%     _______________________________________________________________________
%
%     Disclaimer
%
%     This software is provided 'as-is', without any express or implied
%     warranty. In no event will the author be held liable for any damages
%     arising from the use of this software.
if mod(360,nBins)~=0
    error('The number of bins must be an integer factor of 360')
end
phaseInt = 360/nBins;
histBin = zeros(nBins,size(x,2),size(x,2),size(x,3));
for trialIndex = 1:size(x,3)
    phaseSeries = x(:,:,trialIndex);
    phaseSeries = (180*phaseSeries)./pi; %convert the phase series to degrees
    phaseSeries = phaseSeries+180; %set the bound to [0,360]
    amplitude = y(:,:,trialIndex);
    phaseSeriesFreqIndex = ceil(phaseSeries / phaseInt);
    %for rare cases in which phaseSeriesFreqIndex is exactly zero
    %phaseSeriesFreqIndex(phaseSeriesFreqIndex==0)=nBins;
    %Java function called here
    histBin(:,:,:,trialIndex) = crossfrequencycoupling.CFCalc.fastbuilderjava(phaseSeriesFreqIndex, amplitude, nBins);
end
end