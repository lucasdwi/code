function [allDev] = params(T)
%%
% index for each predictor columns of T
% mnrfit(X,Y) where X = predictors (variables), Y = response (group)
%% Initialize
% Get conditions from T
cond = fieldnames(T);
% Set up alpha range to test
alph = [.001:.001:1];
%% Test alphas
for c = 1:length(cond)
    % Extract predictors and normalize to be >0
    predict = table2array(T.(cond{c})(:,8:end))+1;
    % Columns of response variables
    for j = 2:5 
        % Extract response
        response = table2array(T.(cond{c})(:,j));
        % Cycle through alphas
        for ii = 1:length(alph)
            tic
            % Setup 'opts' structure with alpha value
            opts.alpha = alph(ii);
            for r = 1:20
                % Fit cross validated (3-fold due to size) logistic regression
                CVerr = cvglmnet(predict,response,'binomial','opts','class',3);
                % Get index of lambda with lowest misclassification error
                ind = find(CVerr.lambda == CVerr.lambda_min);
                % Store each minimum error in matrix of cell array
                allDev{c,j}(ii,r) = CVerr.cvm(ind);
            end
            toc
        end
    end
end

%% Extract minimum binomial deviance
%DevMin(j,ii,c) = 
