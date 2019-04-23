% Syntax
%   err = PS_StartStimAllChannels(StimN)
% Description
%   Starts stimulation for all channels for the stimulator StimN.
%   StimN  - stimulator number to configure (starts from 1)
%   Returns:
%       -1 - invalid argument(s)
%       0 - OK
%       1 - device error
%       4 - wrong trigger mode (trigger mode isn't set to PS_TRIG_SOFT (value is 0))


