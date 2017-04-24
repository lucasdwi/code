function [corrNames] = powCorrNames(bands,chan)
%% Generates a string vector of power correlation combination names
% INPUTS
% bands = string vector of band names
% chan = string vector of channel names

% OUTPUTS
% corrNames = column vector of power correlation names

% EXAMPLE
% [corrNames] = powCorrNames({'d','t','a','b','l','h'},{'1','2','3','4'})
% Gets names of data using 6 frequency bands and 4 channels
%%
% bands = {'d','t','a','b','l','h'};
% chan = {'1','2','3','4'};
% Get number of bands and channels
nBands = length(bands);
nChan = length(chan);
% Determine number of correlations within frequency band
sameBlock = nChan^2/2 - nChan/2;
% Determine number of correlations with next frequency bands
diffBlock = nChan*nChan;
% Determine overall number of correlations
values = sameBlock*nBands + diffBlock*(nBands-1);
% Preallocate name vector
corrNames = cell(values,1);
% Start of value counter
vi = 1;
% Cycle through each band
for bi = 1:nBands
    % Cycle through each channel for above band
   for ci = 1:nChan
       % Check if on last band, if so then only use that band 
       if bi == nBands
           biEnd = bi;
       else
           % Otherwise, use this band and the next
           biEnd = bi+1;
       end
       for bi2 = bi:biEnd
           % Check if using two bands and currently at second, fi so set
           % ciStart to 0 so that all channels are used
           if length(bi:biEnd) == 2 && bi2 == biEnd
              ciStart = 0;
           else
               % Otherwise, start at current channel to avoid repeated
               % comparisons (e.g. 1v2 == 2v1)
               ciStart = ci;
           end
           for ci2 = ciStart+1:nChan
               % Concatenate band1-chan1-band2-chan2
               corrNames{vi} = strcat(bands{bi},chan{ci},bands{bi2},chan{ci2});
               % Add to value counter
               vi = vi + 1;
           end
       end
   end
end