function [inds] = logicFind(x,data,op,position)
%% Uses logical indexing to find values within data that have the
% relationship to x defined by op. Requires data to be a row vector
% Can also be used on strings in cells using strcmp as operator
%__________________________________________________________________________
% INPUTS:
% x = value in reference to which to look for
% data = row vector through which to look
% op = operator defining relation with x to look for
% position = which indices to report: first, last, or all; deafult is all
%__________________________________________________________________________
% OUTPUTS:
% inds = indices within 'data' at which 'x' has 'op' relationship with
% 'data'
%__________________________________________________________________________
% USE: 
% inds = logicFind(30,[10,20,30,40,50],'>=') 
% Will return the indices of all values in this vector that are greater
% than or equal to 30, so [3,4,5]
%__________________________________________________________________________
% LLD 2016-17
%%
% If empty position, set to default 'all'
if nargin == 3
    position = 'all';
end
% Check if position is a valid input
if ~strcmpi(position,'all')&&~strcmpi(position,'first')&&~strcmpi(position,'last')
    error([position,' is not a valid input for position.'])
end

allIdx = 1:numel(data);
if strcmp(op,'==')
    if isa(data,'cell')
       inds = allIdx(strcmpi(x,data)); 
    else
        inds = allIdx(data == x);
    end
end
if strcmp(op,'<=')
    inds = allIdx(data <= x);
end
if strcmp(op,'>=')
    inds = allIdx(data >= x);
end
if strcmp(op,'>')
    inds = allIdx(data > x);
end
if strcmp(op,'<')
    inds = allIdx(data < x);
end
if strcmp(op,'~=')
    if isa(data,'cell')
       inds = allIdx(~strcmp(x,data)); 
    else
        inds = allIdx(data ~= x);
    end
end
if strcmpi(position,'first')
    if isempty(inds)
        inds = [];
    else
        inds = inds(1);
    end
elseif strcmpi(position,'last')
    if isempty(inds)
        inds = [];
    else
        inds = inds(end);
    end
end
