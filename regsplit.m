function [parts] = regsplit(string)
%% Splits string into parts starting with capital letter using regexp 
% Includes delimiter in parts

% INPUTS
% string = string to be split

% OUTPUT
% parts = cell array of parts, each cell starting at delimited character
%%
inds = regexp(string,'[A-Z]');
parts = cell(1,numel(inds));
for ii = 1:numel(inds)
    if ii == numel(inds)
        parts{1,ii} = string(inds(ii):end);
    else
        parts{1,ii} = string(inds(ii):inds(ii+1)-1);
    end
end