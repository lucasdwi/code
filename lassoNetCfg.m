function [cfg] = lassoNetCfg(naive,naiveType,rand,normalize,foldGen,cvIterations,minTerm,weights)
%% Generates cfg used in lassoNet.m
%__________________________________________________________________________
% INPUTS: 
% naive = percent of data to withhold as naive test set, if not using a
%   naive set, leave empty []; format: decimal
% naiveType = whether or not to use completely random subset as naive or to
%   generate a balanced naive set with same proportion of binomial
%   responses, only set to 'balanced' in binomial data sets; format:
%   string, either 'balanced' or 'random'
% rand = whether to randomize full data set; used for permutation testing;
%   format: string, 'y' or 'n'
% normalize = whether to normalize data; format: 'y' or 'n'; typically set
%   to 'y'
% foldGen = whether or not to manually generate fold assignments; only
%   needed in cases of binomial data with small ns to make sure that the
%   folds have both cases; format: 'y' or 'n'
% cvIterations = number of cross-validation iterations; format: integer
% minTerm = which error term should be used to pick best model; format:
%   string, either '1se' for lambda+1 standard deviation or 'min' for
%   minimum lambda
% weights = weights to be assigned to predictors; format: column vector
%   with the same number of rows as input data
%__________________________________________________________________________
% OUTPUTS:
% cfg = config structure used in lassoNet.m
%__________________________________________________________________________
%% LLD 2017
%%
if iscell(naive)
    cfg.naive.testX = naive{1,1};
    cfg.naive.testY = naive{1,2};
else
    cfg.naive = naive;
end
cfg.naiveType = naiveType;
cfg.rand = rand;
cfg.normalize = normalize;
cfg.foldGen = foldGen;
cfg.cvIterations = cvIterations;
cfg.minTerm = minTerm;
cfg.weights = weights;
