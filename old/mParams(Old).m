function [allDev] = mParams(T)
%%
% Get conditions from T
cond = fieldnames(T.Base);
% Set up alpha range to test
alph = [.001:.001:1];
% Extract predictors and normalize to be >0
predict = table2array(T.Base(:,8:end))+1;
% Columns of response variables
for j = 2:5
    % Extract response
    response = table2array(T.Base(:,j));
    % Cycle through alphas
    tic
    for ii = 1 %:length(alph)
        for r = 1:20
            tic
            % Fit cross validated (3-fold due to size) logistic regression
            [B,fitInfo] = lassoglm(predict,response,'binomial','Alpha',alph(ii),'CV',5);
            % Get index of lambda with lowest misclassification error
            ind = find(fitInfo.Lambda == fitInfo.LambdaMinDeviance);
            % Store each minimum error in matrix of cell array
            allDev{1,j-1}(ii,r) = fitInfo.Deviance(ind);
            toc
        end
    end
    toc            
end