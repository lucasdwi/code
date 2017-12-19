function data = unpack(data)
%% Unpacks cell array if it only has one entry.
%__________________________________________________________________________
% INPUTS
% data = cell array to be unpacked
%__________________________________________________________________________
% OUTPUTS
% uData = unpacked or unchanged data
%__________________________________________________________________________
% LLD 2017
%%
dim = size(data);
if dim == [1,1]
   data = data{1,1};
end