% Syntax
%   [Mode, err] = PS_GetTriggerMode(StimN)
% Description
%   Gets mode for a stimulator StimN to start stimulation 
%   StimN  - stimulator number to query (starts from 1)
% 	Returns:
% 		Mode: can have the following values:
%           0 - sets stimulator StimN to start stimulation from software;
%           1 - stimulation begins when the digital input for the channel transitions from low (~0V) to high (~5V);
%           2 -  stimulation begins when the digital input transitions from low (~0V) to high (~5V), but if the digital input is still high when the stimulation protocol completes then the stimulation protocol will begin again.
%       err:
%           -1 - invalid argument(s)
%           0 - OK;
