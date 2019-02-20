% Syntax
%   err = PS_AbortAll()
% Description
%   Causes all stimulation to cease immediately even if there is a pulse or arbitrary waveform in progress. This is in contrast to stopping stimulation by calling PS_StopStimChannel or PS_StopStimAllChannels.
%   Returns: 
%       -1 - invalid argument(s)
%       0 – OK, 
%       1 - device error.
