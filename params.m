function [masterDev,masterLam,masterBeta,masterMinLam,masterSurvBeta,masterMeanBeta,masterOR,masterStdBeta,TOrig,betaNames] = params(T,z,m,r,method,resps,surv,k)
%% Tunes and runs elasticnet 
% Inputs:
% T = data table to run elastic net on; format = table (from
%   tabulateData.m)
%   N.B.: reponses variables MUST be the first n columns
% z = normalize with z-score; format = string, 'y' or 'n'
% m = number of master iterations to run whole program through
% r = randomize per iteration; format = string, 'y' or 'n'
% method = type of regression to run; format = either 'binom' or 'gauss'
% resps = columns of responses; format = vector i.e. [1:2]
% surv = survival percent required
%% Initialize
% Set up alpha range to test
alph = (.001:.001:1);
% Save original T (for randomization purposes)
TOrig = T;
% Check method
if ~strcmpi(method,'binom') && ~strcmpi(method,'gauss')
    error('Method not supported; choose either "binom" or "gauss".')
end
%% Run optimizing and testing m times
for mi = 1:m
    % Randomize
    if strcmp(r,'y')
        for c = resps(end)+1:width(T.Base)
            thisPerm = randperm(height(T.Base))';
            T.Base(:,c) = TOrig.Base(thisPerm',c);
        end
    end        
    %% Tune alphas
    % Extract predictors and normalize if wanted
    if strcmp(z,'y')
        predict = zscore(table2array(T.Base(:,resps(end)+1:end)));
    else if strcmp(z,'n')
            predict = table2array(T.Base(:,resps(end)+1:end));
        end
    end
    %%
    % PA allDev; columns = number of models (responses)
    %allDev = cell(4,4);
    allDev = cell(4,length(resps));
    % Columns of response variables
    for j = resps
        % Extract response
        response = table2array(T.Base(:,j));
        % Cycle through alphas
        tic
        for ii = 1:length(alph)
            % Setup 'opts' structure with alpha value
            opts.alpha = alph(ii);
            for r = 1:20
                % Fit cross validated (5-fold due to size) logistic regression
                if strcmpi(method,'binom')
                    CVerr = cvglmnet(predict,response,'binomial','opts','class',k);
                else
                    % Fit cross validated (k-fold) regression
                    if strcmpi(method,'gauss')
                        CVerr = cvglmnet(predict,response,'gaussian','opts','mse',k);
                    end
                end
                % Get index of lambda with lowest misclassification error
                ind = find(CVerr.lambda == CVerr.lambda_min);
                % Store each minimum error in matrix of cell array
                %allDev{1,j-1}(ii,r) = CVerr.cvm(ind);
                allDev{1,j}(ii,r) = CVerr.cvm(ind);
            end
        end
        toc            
    end
    %% Get average for each alpha (row 2), find minimum (row 3) and use index to get first best alpha value (row 4)
    for r = 1:size(allDev,2)
        allDev{2,r} = mean(allDev{1,r},2);
        allDev{3,r} = min(allDev{2,r});
        allDev{4,r} = alph(find(allDev{2,r} == allDev{3,r},1));
    end
    %% Use alpha found above to tune lambda
    % PA allLam, allBeta, and minLam
    allLam = cell(1,size(allDev,2)); allBeta = cell(1,size(allDev,2)); minLam = zeros(size(allDev,2),3);
    for c = 1:size(allDev,2)
        opts.alpha = allDev{4,c};
        %response = table2array(T.Base(:,c+1)); %Weird case because of
        %non-used first column
        response = table2array(T.Base(:,c));
        tic
        for ii = 1:1000
            if strcmpi(method,'binom')
                CVerr = cvglmnet(predict,response,'binomial','opts','class',k);
            else
                if strcmpi(method,'gauss')
                    CVerr = cvglmnet(predict,response,'gaussian','opts','mse',k);
                end
            end
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
        % Find index of least misclassifcation error
        minLam(c,1) = find(allLam{c}(:,4) == min(thisc),1);
        % Save + 1SE lambda with least misclassification error
        minLam(c,2) =  allLam{c}(minLam(c,1),2);
        % Get misclassification error
        minLam(c,3) = allLam{1,c}(minLam(c,1),4);
        toc
    end
    %%
    for c = 1:size(allBeta,2)
        for ii = 1:size(allBeta{1,c},2)
            % Get mean beta
            allBeta{2,c}(1,ii) = mean(nonzeros(allBeta{1,c}(:,ii)));
            % Get std beta
            allBeta{2,c}(2,ii) = std(nonzeros(allBeta{1,c}(:,ii)));
            % Get survival rate
            allBeta{2,c}(3,ii) = sum(allBeta{1,c}(:,ii)~=0)/1000;
        end    
    end
    %% Get survival rates and mean values for all betas
    survBeta = []; mBeta = []; stdBeta = [];
    for c = 1:size(allBeta,2)
        survBeta = vertcat(survBeta,allBeta{2,c}(3,:));
        mBeta = vertcat(mBeta,allBeta{2,c}(1,:));
        stdBeta = vertcat(stdBeta,allBeta{2,c}(2,:));
    end
    %%
    % Set survivals <= surv to zero
    survBeta(survBeta <= surv) = 0;
    % Set means and STD with 0 survival to NaN
    mBeta(survBeta == 0) = NaN;
    stdBeta(survBeta ==0) = NaN;
    %% Transform mean beta into odds ratio (e^beta)
    mOR = exp(mBeta);
    %% Save variables into master structures
    masterDev{mi} = allDev;
    masterLam{mi} = allLam;
    masterBeta{mi} = allBeta;
    masterMinLam{mi} = minLam;
    masterSurvBeta{mi} = survBeta;
    masterMeanBeta{mi} = mBeta;
    masterOR{mi} = mOR;
    masterStdBeta{mi} = stdBeta;
    disp(mi)
    for ii = 1:size(masterSurvBeta{mi},1)
        inds{mi,ii} = logicFind(0,masterSurvBeta{mi}(ii,:),'~=');
        inds{mi,ii} = inds{mi,ii}+length(resps);
        betaNames{mi,ii} = T.Base.Properties.VariableNames(inds{mi,ii});
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