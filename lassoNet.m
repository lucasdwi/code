function [allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x,y,family,type,alph,kfolds,repeats,cfg)
%% Runs elastic net using glmnet functions from the Tibshirani group. Check
% input documentation for details on different options for analysis.
%__________________________________________________________________________
% INPUTS
% x = matrix of predictors; format: observations (animal) x predictor
%     (feature)
% y = matrix of response variables; format: observations (animal) x
%     response; if more than one goes through each separately (perhaps add
%     multinomial later)
% family = type of regression; format: string: 'binomial', 'multinomial',
%          or 'gaussian'
% type = type of error term to use; format: string, 'mse' and 'mae' can be
%        used by all, 'auc' and 'class' can only be used by binomial
% alph = elastic net tuning parameter alpha (0 = ridge, 1 = lasso); format:
%        either integer or vector of values to tune through
% kfolds = number of folds to be used for k-fold CV; format: integer
% repeats = number of times to repeat whole procedure; format: integer
% cfg = configuration structure, includes the following fields
%       naive = percent of data to be left out entirely for final model
%               testing; format: decimal (0.9 for 90%) or two structures
%               with first being x and second being y
%       rand = whether to randomize data; format: string 'y' or 'n',
%              default is 'n'
%       normalize = whether to normalize data; format: string 'y' or 'n',
%                   default is 'n'
%       foldGen = whether to generate foldid manually to ensure each fold
%                 has at least one of each category (used for small
%                 logistic datasets); format: string 'y' or 'n', default
%                 is 'n'
%       cvIterations = number of iterations for cross-validating alpha and
%                      lambda; format: integer, default is 100
%       minTerm = which lambda term to use for optimization; either 'min'
%                 or '1se'
%       weight = weights to apply to observations; format: column vector
%                with same number of rows as x
%__________________________________________________________________________
% OUTPUTS
% allAlpha = structure of alpha tuning data
%          bestAlpha = 
%          meanAlpha = 
%          minAlphaErr = 
% allLambda = structure of lambda tuning data
%           lamErr = error associated with tuned lambdas (either smallest
%               lambda or smallest lambda + 1 S.E.); nx1 where n = number
%               of CV iterations
%           lamErrAvg = mean of lamErr
%           lamErrSd = standard deviation of lamErr
%           minLam = values of tuned lambdas (either smallest or +1 S.E.);
%               nx1 where n = number of CV iterations
%           minLamErr = error associated with bestLambda
%           bestLambda = lambda with smallest error in minLam (max if auc)
%           bestLambdaInds = index of bestLambda
%           allErr = average of the mean cross-validated errors per
%               iteration; 1 value comes from each iteration which is the
%               mean of that iterations cvms
%           allErrAvg = mean of allErr
%           allErrSd = standard deviation of allErr
% allBeta = structure of beta data
%         betas = 
%         survBeta = 
% CVfits = cell array of data structure outputs of cvglmnet while tuning
%        lambda; contains fit object needed for using cvglmnetPredict
% acc = 
% hist = history structure containing information on function inputs
%% Consider using an alpha of 1-epsilon where epsilon is near zero term rather than 1 for the case of lasso
%% Check inputs and fill out cfg as necessary
% Check for equal numbers of independent variables in x and y
if size(x,1) ~= size(y,1)
    error('Inequal predictor and response matrices.')
end
% Check regression family
if ~strcmpi(family,'binomial') && ~strcmpi(family,'gaussian') && ~strcmpi(family,'multinomial') && ~strcmpi(family,'poisson')
    error(['Family, ',family,', not supported; choose either "binomial", "multinomial", "gaussian", or "poisson".'])
end
% Check type of error matches regression family
if (strcmpi(type,'auc') || strcmpi(type,'class')) && ...
        (~strcmpi(family,'binomial') && ~strcmpi(family,'multinomial') ...
        && ~strcmpi(family,'poisson'))
   error(['Can not use type ',type,' with family ',family,'.'])
end
% Check nfolds
if kfolds <=3
    warning('It is not recommended to run cross-validation with less than 4 folds, consider resetting. Press "Ctrl+C" to quit, or any key to continue.')
    %pause
end
% Set cfg.rand to default if unset: 'n'
if ~isfield(cfg,'rand') || isempty(cfg.rand)
    cfg.rand = 'n';
end
% Set cfg.normalize to default if unset: 'y' 
if ~isfield(cfg,'normalize') || isempty(cfg.normalize)
    cfg.norm = 'y';
end
% Set cfg.foldGen to default if unset: 'n'
if ~isfield(cfg,'foldGen') || isempty(cfg.foldGen)
    cfg.foldGen = 'n';
end
% Set cfg.cvIterations to default if unset: 100
if ~isfield(cfg,'cvIterations') || isempty(cfg.cvIterations)
    cfg.Iterations = 100;
end
% Set cfg.naiveType to default if unset and cfg.naive is set: 'random'
if (~isfield(cfg,'naiveType') || isempty(cfg.naiveType)) && (isfield(cfg,'naive') || ~isempty(cfg.naive))
   cfg.naiveType = 'random'; 
end
% Set cfg.minTerm to defualt if unset: '1se'
if ~isfield(cfg,'minTerm') || isempty(cfg.minTerm)
    cfg.minTerm = '1se';
end
%% Split data into naive-test set and training set
% Use CV to train best model, then get error from predicting naive-test set
for rep = 1:repeats
    if ~isempty(cfg.naive)
        % Check if using a numeric, percent of complete dataset to use
        if isnumeric(cfg.naive)
            if strcmpi(cfg.naiveType,'balanced')
                % Maintain distribution of response in test and train sets
                % Get percent of ones in data
                nOne = sum(y==1);
                oneDist = nOne/(numel(y));
                % Determine size of naive set, rounding up
                nNaive = ceil(numel(y)*cfg.naive);
                % Calculate number of ones to include in test set
                nNaiveOnes = ceil(nNaive*oneDist);
                % Grab random subset of ones for naive testing
                oneInd = logicFind(1,y,'==');
                rng('shuffle')
                randOneInd = randperm(nOne,nNaiveOnes);
                naiveOneInd = oneInd(randOneInd);
                % Grab random subset of zeros to complete naive test set
                nZero = sum(y==0);
                nNaiveZero = nNaive-nNaiveOnes;
                zeroInd = logicFind(0,y,'==');
                rng('shuffle')
                randZeroInd = randperm(nZero,nNaiveZero);
                naiveZeroInd = zeroInd(randZeroInd);
                % Combine indices
                naiveInd = [naiveOneInd';naiveZeroInd'];
            elseif strcmpi(cfg.naiveType,'random')
                % Determine indices of random subset for naive testing; rounds down
                rng('shuffle')
                naiveInd = randperm(size(x,1),floor(size(x,1)*cfg.naive));
            end
            % Use naiveInd to build test and train sets
            testX = x(naiveInd,:);
            testY = y(naiveInd,:);
            if ~isempty(cfg.weights)
                testWeights = cfg.weights(naiveInd)';
            end
            % Set up index vector
            inds  = 1:size(x,1);
            % Grab train data using those indices not in naiveInd
            trainX = x(~ismember(inds,naiveInd),:);
            trainY = y(~ismember(inds,naiveInd),:);
            if ~isempty(cfg.weights)
                trainWeights = cfg.weights(~ismember(inds,naiveInd))';
            end
        % Otherwise, use dataset provided as test set and complete data in
        % x as training set
        elseif isstruct(cfg.naive)
            % Set up training set
            trainX = x;
            trainY = y;
            % Set up testing set
            testX = cfg.naive.testX;
            testY = cfg.naive.testY;
            % Set up empty naiveInd
            naiveInd = [];
        end
    else
        % If not using naive set, rename x and y and create empty test matrices
        trainX = x; testX = [];
        trainY = y; testY = [];
        if ~isempty(cfg.weights)
            trainWeights = cfg.weights;
        end
        % Create empty predY and naiveInd
        predY = []; naiveInd = [];
    end
    % Determine number of response (dependent) variables - number of models
    nResp = size(trainY,2);
    % Determine number of independent variables
    nVar = size(trainX,2);
    % Determine number of observations - number of animals
    nObv = size(trainY,1);
    %% Shuffle data assignment to response by randomizing y
    if strcmpi(cfg.rand,'y')
        rng('shuffle')
        thisPerm = randperm(nObv);
        trainY = trainY(thisPerm',:);
        if ~isempty(cfg.weights)
            trainWeights = trainWeights(thisPerm);
        end
    end
    %% Normalize (z-score) data
    if strcmpi(cfg.normalize,'y')
        trainX = zscore(trainX);
        testX = zscore(testX);
    end
    %% For ENET generate folds s.t. the same divisions are used in tuning alpha and lambda
    foldid = zeros(cfg.cvIterations,nObv,nResp);
    for r = 1:nResp
        for ii = 1:cfg.cvIterations
            if strcmpi(cfg.foldGen,'y')
                foldid(ii,:,r) = foldPerm(kfolds,trainY(:,r));
            else
                N = nObv;
                population = cat(2, repmat(1:kfolds, 1, floor(N/kfolds)), 1:mod(N,kfolds));
                rng('shuffle')
                foldid(ii,:,r) = population(randperm(length(population), N));
            end
        end
    end
    %% Tune alpha if given a range of values
    if length(alph) ~= 1
        % Preallocate
        minAlphaErr = zeros(length(alph),cfg.cvIterations,nResp);
        meanMinAlpha = zeros(length(alph),nResp);
        minMeanAlphaErr = zeros(1,nResp);
        bestAlpha = zeros(1,nResp);
        % Cycle through each response variable
        for r = 1:nResp
            % Cycle through each value of alpha
            for c = 1:size(foldid,1)
                disp(['Tuning alpha through cross-validation: ',num2str(c),' out of ',num2str(size(foldid,1)),' for response variable ',num2str(r),'...'])
                for a = 1:length(alph)
                    opts.alpha = alph(a);
                    if ~isempty(cfg.weights)
                        opts.weights = trainWeights;
                    end
                    CVerr = cvglmnet(trainX,trainY(:,r),family,opts,type,[],foldid(c,:,r));
                    % Get index of lambda with lowest misclassification error
                    % (using min rather than min + 1 SE since size of lambda is
                    % not an issue at this step.)
                    ind = logicFind(CVerr.lambda_min,CVerr.lambda,'==');
                    % Store each minimum error in matrix of cell array
                    minAlphaErr(a,c,r) = CVerr.cvm(ind);
                end
            end
            % Get mean errors for each alpha across CVs
            meanMinAlpha(:,r) = mean(minAlphaErr(:,:,r),2);
            % Find minimum error and index thereof
            [minMeanAlphaErr(1,r),minAlphaInd] = min(meanMinAlpha(:,r));
            % Get value of alpha based on above index
            bestAlpha(1,r) = alph(minAlphaInd);
        end
    end
    %% Use k-fold CV; accumulate error across folds.
    % Standard cvglmnet reports average and sd of CV error. Gives mean,
    % upper/lower bounds for each lambda - default of 100 lambda but ends
    % earlier if no significant difference in error.
    
    % Use either min lambda or min + 1 SE lambda (min+1SE gives larger lambda
    % thus the greater the shrinkage and the less complex the model)
    tic
    disp('Tuning lambda...')
    % Preallocate
    CVfits = cell(cfg.cvIterations,nResp);
    bestLambda = zeros(1,nResp);
    bestLambdaInds = zeros(1,nResp);
    bestLamErr = zeros(1,nResp);
    betas = zeros(cfg.cvIterations,nVar,nResp);
    survBeta = zeros(nResp,nVar);
    signBeta = zeros(nResp,nVar);
    % Cycle through each response variable
    for r = 1:nResp
        disp(['Model ',num2str(r),' out of ',num2str(nResp),'...'])
        % Set alpha to either given value, or tuned value from above
        if length(alph) == 1
            opts.alpha = alph;
        else
            opts.alpha = bestAlpha(1,r);
        end
        % If weights exist, add to opts structure
        if ~isempty(cfg.weights)
           opts.weights = trainWeights; 
        end
        % Cycle through cross-validation folds
        for c = 1:cfg.cvIterations
            % Save fit structure for case of predicting naive data
            CVfits{c,r} = cvglmnet(trainX,trainY(:,r),family,opts,type,[],foldid(c,:,r));
        end
    end
    % Go through CVfits and extract errors
    if strcmpi(cfg.minTerm,'1se')
        lamErr = cellfun(@(exMinErr) exMinErr.cvm(logicFind(exMinErr.lambda_1se,exMinErr.lambda,'==')),CVfits);
    elseif strcmpi(cfg.minTerm,'min')
        lamErr = cellfun(@(exMinErr) exMinErr.cvm(logicFind(exMinErr.lambda_min,exMinErr.lambda,'==')),CVfits);
    end
    % Get average and standard deviation of min error for each model
    lamErrAvg = mean(lamErr,1);
    lamErrSd = std(lamErr,[],1);
    % Go through CVfits and extract lamda_1se with lowest misclassification
    % error for naive prediction
    if strcmpi(cfg.minTerm,'1se')
        minLam = cellfun(@(exLam) exLam.lambda_1se,CVfits);
    elseif strcmpi(cfg.minTerm,'min')
        minLam = cellfun(@(exLam) exLam.lambda_min,CVfits);
    end
    % Go through CVfits and extract mean of all error - cvm
    allErr = cellfun(@(exAllErr) mean(exAllErr.cvm),CVfits);
    % Get average and standard deviation of all errors for each model
    allErrAvg = mean(allErr,1);
    allErrSd = std(allErr,[],1);
    % Find lambda_1se (value and index) with smallest error - in case of
    % ties, choose larger lambda
    for r = 1:nResp
        % Get indices of best lamErr (minima for 'error' terms:
        % 'mae','mse', and 'class'; maxima for 'accuracy' term: 'auc')
        if strcmpi(type,'auc')
            bestInds = logicFind(max(lamErr(:,r)),lamErr(:,r),'==');
        else
            bestInds = logicFind(min(lamErr(:,r)),lamErr(:,r),'==');
        end
        % Get lambda value and index of largest if there is a tie
        [bestLambda(r),thisBestInd] = max(minLam(bestInds,r));
        bestLambdaInds(r) = bestInds(thisBestInd);
        % Get error of best lambda
        bestLamErr(r) = lamErr(bestLambdaInds(r),r);
    end
    % Get beta coefficients of lambda_1se (bestLam) - N.B. these values are
    % not comparable due to varying shrinkage, only used to compute
    % survival rates
    for r = 1:nResp
        for c = 1:cfg.cvIterations
            betas(c,:,r) = CVfits{c,r}.glmnet_fit.beta(:,logicFind(minLam(c,r),CVfits{c,r}.lambda,'=='));
        end
        survBeta(r,:) = mean(betas(:,:,r)~=0);
        % Also get sign of betas for predictive directionality
        signBeta(r,:) = sign(mean(betas(:,:,r)));
    end
    toc
    %% Use above models to predict naive dataset, if extistent
    if ~isempty(cfg.naive)
        nPred = size(testY,1);
        % Preallocate
        predY = zeros(1,nPred);
        acc = [];%zeros(1,nResp);
        % Use best CVfit for each model
        for r = 1:nResp
            % If family is gaussian (continuous), use defualt prediction
            % 'type' ('link')
            if strcmpi(family,'gaussian') || strcmpi(family,'poisson')
                [predY(r,:)] = cvglmnetPredict(CVfits{bestLambdaInds(r),r},testX,['lambda_',cfg.minTerm]);
                % Compare predY to testY: difference
                acc(r,:) = predY(r,:) - testY(:,r)';
            % Otherwise, if binomial or multinomial, use 'response' type to
            % get probability of group 1 assignment
            elseif strcmpi(family,'binomial') || strcmpi(family,'multinomial')
                [predY(r,:)] = cvglmnetPredict(CVfits{bestLambdaInds(r),r},testX,['lambda_',cfg.minTerm],'response');
                [~,~,~,acc(r)] = perfcurve(testY,predY(r,:),1);
            end
        end
    else
        acc = [];
    end
    %% Create output structures
    % Preallocate on first rep
    if rep == 1
        allAlpha = cell(1,repeats);
        allLambda = cell(1,repeats);
        allBeta = cell(1,repeats);
        accArray = cell(1,repeats);
        cvFitsArray = cell(1,repeats);
    end
    % Alpha structure
    if size(alph,2) ~= 1
        allAlpha{rep}.bestAlpha = bestAlpha;
        allAlpha{rep}.meanAlpha = meanMinAlpha;
        allAlpha{rep}.minMeanAlphaErr = minMeanAlphaErr;
        allAlpha{rep}.minAlphaErr = minAlphaErr;
    else
        allAlpha{rep} = [];
    end
    % Lambda structure
    allLambda{rep}.lamErr = lamErr;
    allLambda{rep}.lamErrAvg = lamErrAvg;
    allLambda{rep}.lamErrSd = lamErrSd;
    allLambda{rep}.minLam = minLam;
    allLambda{rep}.minLamErr = bestLamErr;
    allLambda{rep}.bestLambda = bestLambda;
    allLambda{rep}.bestLambdaInds = bestLambdaInds;
    allLambda{rep}.allErr = allErr;
    allLambda{rep}.allErrAvg = allErrAvg;
    allLambda{rep}.allErrSd = allErrSd;
    % Beta structure
    allBeta{rep}.betas = betas;
    allBeta{rep}.survBeta = survBeta;
    allBeta{rep}.signBeta = signBeta;
    % Create acc and CVfits array
    accArray{rep}.acc = acc;
    accArray{rep}.pred = predY;
    cvFitsArray{rep} = CVfits;
    % Create hist array
    hist.cfg = cfg;
    hist.family = family;
    hist.type = type;
    hist.nfolds = kfolds;
    hist.alpha = alph;
    hist.repeats = repeats;
    hist.testInd = naiveInd;
end