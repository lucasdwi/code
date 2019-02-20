% Syntax
%   err = PS_StartStimChannel(StimN, ChN)
% Description
%   Starts stimulation for channel ChN with previously configured parameters.
%   StimN  - stimulator number to configure (starts from 1)
%   ChN  - channel number to start stimulation (starts from 1)
%   Returns:
%       -1 - invalid argument(s)
%       0 - OK
%       1 - device error
%       4 - wrong trigger mode (trigger mode doesn't set to PS_TRIG_SOFT (value is 0))
