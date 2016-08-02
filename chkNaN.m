function [nNaN,indSkp] = chkNaN(LFPTsNaN,chans,NaNcutoff)
%% Checks if any of the channels are much noisier than others
% If channel has more NaNs than allowed by cutoff; NaNs entire channel and
% provides indices of channels to skip.

% Inputs:
% LFPTsNaN = LFPTsNaN data structure; from threshFilt.m
% chans = number of channels; from threshFilt.m
% NaNcutoff = point at which channel has too many NaNs; format: standard
%   deviations (e.g. 1.5)

% Outputs:
% nNaN = matrix of number of NaNs per channel (row 1) and number of
%   standard deviations from mean (row2)
% indSkp = indices of those channels with number of NaNs >= cutoff
%   N.B. used by threshFilt.m
%% Find channels with more NaNs than cutoff allows
nNaN = zeros(2,chans);
indSkp = [];
for i = 1:chans
    nNaN(1,i) = sum(isnan(LFPTsNaN.data(i,:))); %counts number of NaNs per channel
end
mNaN = mean(nNaN(1,:)); %Mean number of NaNs per channel
sNaN = std(nNaN(1,:)); %Standard deviation of number of NaNs per channel
for i = 1:chans
    nNaN(2,i) = abs(mNaN - nNaN(1,i))/sNaN;
    if nNaN(2,i) >= NaNcutoff
        %LFPTsNaN.data(i,:) = NaN;
        indSkp = horzcat(indSkp,i);
    end
end