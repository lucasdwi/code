function [LFPTs,adfreq] = dwnSample(LFPTs,dsf,adfreq,chans)
%% Decimates LFPTs data and downsamples time vector
% Inputs:
% LFPTs = data structure
% dsf = downsample factor; format = integer 
% adfreq = sampling rate; format = Hz
% chans = number of channels; format = integer

% Outputs:
% LFPTs = downsampled data
% adfreq = new sampling rate
%%
LFPTsDwn.data = [];

for i = 1:chans
    LFPTsDwn.data(i,:) = decimate(LFPTs.data(i,:),dsf);
end
LFPTsDwn.tvec = downsample(LFPTs.tvec,dsf);
% Overwrite old LFPTs with LFPTsdwn; for further analysis
LFPTs.data = LFPTsDwn.data;
LFPTs.tvec = LFPTsDwn.tvec;
adfreq = adfreq/dsf;
