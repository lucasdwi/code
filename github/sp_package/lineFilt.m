function [LFPTs] = lineFilt(LFPTs,start,stop)
% Filter out 60Hz noise with given start and stop frequencies
% INPUTS
% LFPTs = data structure; format LFPTs.data = channel x time matrix
% start = beginning edge of stop-band; format = Hz
% stop = ending edge of stop-band; foamt = Hz
%%
[b,a] = cheby1(2,.5,[start stop]*2/adfreq,'stop');
% Plot filter
%fvtool(b,a,'Fs',adfreq);
%figure; plot(LFPTs.data(1,1:20000));
% Setup filtered LFPT structure
LFPTsFilt = LFPTs;
% Apply filter to all chans
for i = 1:chans
    LFPTsFilt.data(i,:) = filtfilt(b,a,LFPTsFilt.data(i,:));
end
% Overwrite LFPTs with LFPTsFilt
LFPTs = LFPTsFilt;