function [LFPTs,adfreq] = dwnSample(LFPTs,dsf,adfreq,chans)
%% Decimates LFPTs data and downsamples time vector
% Inputs:
% LFPTs = LFP data structure
% dsf = downsample factor 
% adfreq = sampling rate
% chans = number of channels 
% Outputs:
% LFPTs = downsampled LFP data
% adfreq = new sampling rate
%%
dwnLFPTs.data = [];

for i = 1:chans
    dwnLFPTs.data(i,:) = decimate(LFPTs.data(i,:),dsf);
end
dwnLFPTs.tvec = downsample(LFPTs.tvec,dsf);
% Overwrite old LFPTs with dwnLFPTs; for further analysis
LFPTs.data = dwnLFPTs.data;
LFPTs.tvec = dwnLFPTs.tvec;
adfreq = adfreq/dsf;
