function [foldid] = foldPerm(kfolds,response)
%% Generates pseudorandom k-fold indices, useful for cross validation of small datasets with logistic or nominal reponses
% N.B.: code was built for binary/logisitc responses; if number of levels
% >2 consider validity of this method

% INPUTS:
% nfolds = the number of folds; format = integer
% response = vector of response variables, used to verfiy that each fold
%   has at least one of each group; format = column vector

% OUTPUTS: 
% foldid = vector of indices (1:kfolds) of the same size as response where
%   indices refer to which group each response belongs to

%%
N = size(response,1);
population = cat(2, repmat(1:kfolds, 1, floor(N/kfolds)), 1:mod(N,kfolds));
repeat = 'y';
count = 0;
while strcmp(repeat,'y')
    randperm(length(population),N);
    foldid = population(randperm(length(population), N));
    for ii = 1:kfolds
        chksize(ii) = size(unique(response(foldid==ii)),1);
    end
    if all(chksize==2)
        repeat = 'n';
    end
    count = count +1;
end