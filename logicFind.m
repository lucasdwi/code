function [inds] = logicFind(x,data,op)
%% Uses logical indexing to find values within data that have the
% relationship to x defined by op. Requires data to be a row vector
%% INPUTS
% x = value in reference to which to look for
% data = row vector through which to look
% op = operator defining relation with x to look for

% Example: logicFind(30,[10,20,30,40,50],'>=') will return the indices of all
% values in this vector that are greater than or equal to 30, so [3,4,5]

allIdx = 1:numel(data);
if strcmp(op,'==')
    inds = allIdx(data == x);
end
if strcmp(op,'<=')
    inds = allIdx(data <= x);
end
if strcmp(op,'>=')
    inds = allIdx(data >= x);
end
if strcmp(op,'~=')
    inds = allIdx(data ~= x);
end
