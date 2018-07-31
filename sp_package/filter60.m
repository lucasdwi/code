function [filtData] = filter60(data,fs,display,varargin)
%% Uses second order Chebychev type I notch filter to remove 60 cycle line 
% noise
%__________________________________________________________________________
% INPUTS:
% data = data to filter; can either be plain matrix or in structure in
%   .data field; format: in either case chan X time
% fs = sampling frequency of data; format: integer 
% display = whether to plot filter; format: either 1 or 0
%__________________________________________________________________________
% OUTPUTS:
% filtData = filtered data in same format as input; format: channel X time
%__________________________________________________________________________
% USE:
% [filtData] = filter60(LFPTs,2000,'off')
% Removes 60 cycle line noise from data within LFPTs structure which was
% sampled at 2 kHZ and does not display filter used.
%__________________________________________________________________________
% LLD 2016-17
%% Parse inputs
p = inputParser;
p.CaseSensitive = false;
addRequired(p,'data');
addRequired(p,'fs');
addRequired(p,'display');

addParameter(p,'low',59);
addParameter(p,'high',61);
addParameter(p,'order',2);
addParameter(p,'ripple',05);
parse(p,data,fs,display,varargin{:})
%% Check if input data is structure or matrix and put into toFilt variable
if isa(data,'struct')
   toFilt = data.data; 
elseif isa(data,'double')
   toFilt = data;
end
%%
% Create second order Chebychev Type I stopband filter from 59 to 61 Hz
% with 0.5 dB of ripple
[b,a] = cheby1(p.Results.order,p.Results.ripple,[p.Results.low p.Results.high]*2/fs,'stop');
% Plot filter
if display == 1
    fvtool(b,a,'Fs',fs);
end
% Count number of channels
chans = size(toFilt,1);
% Preallocate filtData
filtData = zeros(size(toFilt));
% Apply filter to all chans
for iC = 1:chans
    filtData(iC,:) = filtfilt(b,a,toFilt(iC,:));
end
% If a data structure (i.e. LFPTs) record that data was filtered
if isa(data,'struct')
    if ~isfield(data,'cfg')
        data.cfg = [];
    else
        data.cfg.history.mfun{end+1} = 'filter60';
    end
end
