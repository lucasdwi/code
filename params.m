function [allDev] = params(T)
%%
% index for each predictor columns of T
% mnrfit(X,Y) where X = predictors (variables), Y = response (group)
%% Initialize
% Double table
% row2 = T.Base.Properties.RowNames;
% row2 = strcat(row2,'_2');
% T2.Base = T.Base;
% T2.Base.Properties.RowNames = row2;
% Tx2.Base = vertcat(T.Base,T2.Base);
%%
% Get conditions from T
cond = fieldnames(T.Base);
% Set up alpha range to test
alph = [.001:.001:1];
%% Test alphas
%for c = 1:length(cond)
    % Extract predictors and normalize to be >0
    predict = table2array(T.Base(:,8:end))+1;
    % Columns of response variables
    for j = 2:5
            % Extract response
            response = table2array(T.Base(:,j));
            % Cycle through alphas
            tic
            for ii = 1:length(alph)
                % Setup 'opts' structure with alpha value
                opts.alpha = alph(ii);
                for r = 1:20
                    % Fit cross validated (3-fold due to size) logistic regression
                    CVerr = cvglmnet(predict,response,'binomial','opts','class',3);
                    % Get index of lambda with lowest misclassification error
                    ind = find(CVerr.lambda == CVerr.lambda_min);
                    % Store each minimum error in matrix of cell array
                    allDev{1,j-1}(ii,r) = CVerr.cvm(ind);
                end
            end
            toc            
    end
%end


%% Get average for each alpha (row 2), find minimum(-a) (row 3) and use index to get first best alpha value (row 4)
for r = 1:length(allDev)
    allDev{2,r} = mean(allDev{1,r},2);
    allDev{3,r} = min(allDev{2,r});
    allDev{4,r} = alph(find(allDev{2,r} == allDev{3,r},1));
end
%% Use alpha found above to tune lambda
for c = 3:length(allDev)
    opts.alpha = allDev{4,c};
    response = table2array(T.Base(:,c+1));
    tic
    for ii = 1:1000
        CVerr = cvglmnet(predict,response,'binomial','opts','class',3);
        % Save minimum lambda
        allLam{c}(ii,1) = CVerr.lambda_min;
        % Save lamba +1 SE from min
        allLam{c}(ii,2) = CVerr.lambda_1se;
        % Save index of +1SE lambda
        allLam{c}(ii,3) = find(CVerr.lambda == CVerr.lambda_1se);
        % Save misclassification error for +1SE lambda
        allLam{c}(ii,4) = CVerr.cvm(allLam{c}(ii,3));
        % Save betas at lambda +1 SE from min
        allBeta{c}(ii,:) = CVerr.glmnet_fit.beta(:,allLam{c}(ii,3))';
    end
    thisc = allLam{c}(:,4);
    % Find index of least non-zero misclassifcation error
    minLam(c,1) = find(allLam{c}(:,4) == min(thisc(thisc>0)),1); 
    % Save +1SE lambda with least non-zero misclassification error
    minLam(c,2) =  allLam{c}(minLam(c,1),2);
    toc
end
%%
for c = 1:size(allBeta,2)
    for ii = 1:size(allBeta{1,c},2)
        % Get mean beta
        allBeta{2,c}(1,ii) = mean(allBeta{1,c}(:,ii));
        % Get survival rate
        allBeta{2,c}(2,ii) = sum(allBeta{1,c}(:,ii)~=0)/1000;
    end
end
% %% Manual filling - alpha
% for j = 3
%     z = find(not(allDev{1,j-1}),1);
%     s = size(allDev{1,j-1},1);
%     while  s < 1000 
%         try
%         response = table2array(T.Base(:,j));
%         % Cycle through alphas
%         for ii = s:length(alph)
%             % Setup 'opts' structure with alpha value
%             opts.alpha = alph(ii);
%             % Start at beginning of row, first column
%             for r = 1:20
%                 % Fit cross validated (3-fold due to size) logistic regression
%                 CVerr = cvglmnet(predict,response,'binomial','opts','class',3);
%                 % Get index of lambda with lowest misclassification error
%                 ind = find(CVerr.lambda == CVerr.lambda_min);
%                 % Store each minimum error in matrix of cell array
%                 allDev{j-1}(ii,r) = CVerr.cvm(ind);
%             end
%         end
%         catch err
%             if strcmp(err.identifier, 'Index exceeds matrix dimensions.')
%                 s = size(allDev{1,j-1},1);
%             end
%         end
%         s = size(allDev{1,j-1},1);
%         disp(s)
%     end
% end
% %% Manual filling - lambda
% for j = 5
%     z = find(not(allLam{1,j-1}),1);
%     s = size(allLam{1,j-1},1);
%     while  s < 1000 
%         try
%         response = table2array(T.Base(:,j));
%         % Cycle through alphas
%         for ii = s:1000
%             CVerr = cvglmnet(predict,response,'binomial','opts','class',3);
%             % Save minimum lambda
%             allLam{c}(ii,1) = CVerr.lambda_min;
%             % Save lamba +1 SE from min
%             allLam{c}(ii,2) = CVerr.lambda_1se;
%             % Save index of +1SE lambda
%             allLam{c}(ii,3) = find(CVerr.lambda == CVerr.lambda_1se);
%             % Save misclassification error for +1SE lambda
%             allLam{c}(ii,4) = CVerr.cvm(allLam{c}(ii,3));
%             % Save betas at lambda +1 SE from min
%             allBeta{c}(ii,:) = CVerr.glmnet_fit.beta(:,allLam{c}(ii,3))';
%         end
%         catch err
%             if strcmp(err.identifier, 'Index exceeds matrix dimensions.')
%                 s = size(allLam{1,j-1},1);
%             end
%         end
%         s = size(allLam{1,j-1},1);
%         disp(s)
%     end
% end