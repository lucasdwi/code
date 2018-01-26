function [tierStart] = tier(x)
%% Finds tiers within data by finding when the data becomes significantly 
% different from the first value in the tier. For example, the first tier
% are all the values which are not significantly different from the first
% data entry. 
%__________________________________________________________________________
% INPUTS
% x = data array to be tiered; format = samples X variable. First dimension
% has all samples of the variable. Needs to have many values to be able to
% run t-tests. N.B.: x NEEDS to be sorted in descending order before
% running.
%__________________________________________________________________________
% OUTPUTS
% tierStart = integer vector with the start indices of each tier
%__________________________________________________________________________
% LLD 2018
%%
% Set compare index to 1
compInd = 1;
% Start teirs with 1
tierStart = 1;
% Go through values until end of data matrix
while compInd<size(x,2)
    % Cycle through data running ttests - puts zeros from 1:compInd if
    % compInd ~= 1
    for ii = compInd:size(x,2)
        [~,p(ii)] = ttest2(x(:,compInd),x(:,ii));
    end
    % Replace 0s with NaN to avoid triggering on those indices
    p(p==0) = NaN;
    % Apply Bonferroni correction but only take into condsideration actual
    % number of comparisons
    p = p.*(size(x,2)-compInd);
    % Find when data becomes significantly different from compInd, update
    % compInd
    compInd = logicFind(0.05,p,'<=','first');
    % Add teir start value - but if empty, then end by setting compInd to
    % size(x,2)
    if isempty(compInd)
        compInd = size(x,2);
    else
        tierStart = [tierStart,compInd]; %#ok<AGROW>
    end
    % Clear p for next iteration
    clear p
end
