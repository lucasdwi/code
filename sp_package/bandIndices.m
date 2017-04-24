function [bInd] = bandIndices(bands,f,varargin)
%% Finds indices within frequency vector corresponding to frequency band 
% starts and stops. Uses index that is >= for start and <= for stop. 
% N.B.: Overlapping frequency bands will give warning.
%
% USE:
% [bandInd] = bandIndices({'delta',[1 4];'low gamma',[45 65]},1:2:100)
% Will return:
% [1 2;23 33]
%__________________________________________________________________________
% INPUTS:
% bands = structure with bands of interest starting with lowest; format:
%   {'band1',[lower upper];'band2',[lower upper];...}
% f = frequency vector; format: row or column vector of numbers
% varargin{1} = whether or not to check for overlapping; 1 = check, 0 =
%   don't check
%__________________________________________________________________________
% OUTPUTS:
% bInd = matrix of indices at which the lower and upper (start/stop) points
%   from 'bands' fall within 'f'
%__________________________________________________________________________
% LLD 2017
%%
bInd = zeros(size(bands));
for iB = 1:size(bands,1)
    bInd(iB,1) = logicFind(bands{iB,2}(1),f,'>=','first');
    bInd(iB,2) = logicFind(bands{iB,2}(2),f,'<=','last');
end 
if ~isequal(numel(bInd),numel(unique(bInd))) && varargin{1} == 1
    disp(['Warning: Band indices overlap in frequency vector! This may '...
        'lead to over-representing data.'])
    disp('Press any key to continue or Ctrl+C to quit.')
    pause
end  