function [LFPTs,adfreq] = dwnSample(LFPTs,dsf,adfreq)
%% Decimates LFPTs data and downsamples time vector
%__________________________________________________________________________
% INPUTS:
% LFPTs = data structure, with .data as field with data in it
% dsf = downsample factor, needs to be a divisor of adfreq; format: integer 
% adfreq = sampling rate; format: integer
%__________________________________________________________________________
% OUTPUTS:
% LFPTs = decimated data
% adfreq = adjusted sampling rate
%__________________________________________________________________________
% USE:
% [LFPTs,adfreq] = dwnSample(LFPTs,5,2000);
% Decimates data in LFPTs which was sampled at 2 kHz by a factor of 5 and
% adjusts sampling rate accordingly (2000/5=400)
%__________________________________________________________________________
% LLD 2016-17
%% Preallocate counts
% Count number of channels
chans = size(LFPTs.data,1);
% Count number of data points
tStamps = size(LFPTs.data,2);
%% Check appropriate dsf with data
% Check that dsf is divisor of adfreq
if mod(adfreq,dsf) ~= 0
    error([num2str(adfreq) ' (adfreq) is not divisible by ',num2str(dsf),' (dsf).'])
end
% Check that dsf is divisor of data length, if it isn't then truncate data
if mod(tStamps,dsf) ~= 0
    LFPTs.data = LFPTs.data(:,1:end-mod(tStamps,dsf));
    LFPTs.tvec = LFPTs.tvec(:,1:end-mod(tStamps,dsf));
    tStamps = size(LFPTs.data,2);
end
%% Set up output and downsample 
% Set up output data
dwnData = zeros(chans,tStamps/dsf);
for iC = 1:chans
    dwnData(iC,:) = decimate(LFPTs.data(iC,:),dsf);
end
dwnTvec = downsample(LFPTs.tvec,dsf);
% Overwrite old LFPTs with LFPTsdwn; for further analysis
LFPTs.data = dwnData;
LFPTs.tvec = dwnTvec;
% Adjust sampling rate
adfreq = adfreq/dsf;
% Record that data was downsampled
% LFPTs.cfg.history.mfun{end+1} = 'dwnSample';
