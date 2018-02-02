function [greek] = latinToGreek(latin)
%% Swaps Latin spellings of Greek letters with Matlab special characters
%__________________________________________________________________________
% INPUTS
% latin = cell string array with Latin spellings of Greek letters; format =
%   cell array in which each letter is separated
%__________________________________________________________________________
% OUTPUTS
% greek = Matlab special character equivalent of Greek letters
%__________________________________________________________________________
% LLD 2018
%% Cycle through each element of Latin and use the corresponding Greek
% Preallocate greek
greek = cell(size(latin));
for ii = 1:length(latin)
    % Use capital delta
    if strcmpi(latin{ii},'delta')
        greek{ii} = '\Delta';
    % Add '\' in front of theta, alpha, and beta
    elseif any(strcmpi(latin{ii},{'theta','alpha','beta'}))
        greek{ii} = ['\',latin{ii}];
    % Detect low and high gamma, split off 'l' or 'h', add '\', then add
    % back 'gamma'
    elseif any(strcmpi(latin{ii},{'lgamma','hgamma'}))
        greek{ii} = [latin{ii}(1),'\',latin{ii}(2:end)];
    end
end
    