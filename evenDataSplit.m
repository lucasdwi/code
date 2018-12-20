function [all,each,rnd,weights] = evenDataSplit(data,trainSize,testSize,method,rep)
%% Splits data into training and testing sets. Forces the training set into 
% a 50/50 split, while the testing set exhibits the distribution of classes
% given. N.B. first column will be treated as minority class.
%__________________________________________________________________________
% INPUTS:
% data = input data; format: cell array where rows represent subsets (e.g.
%   animals or days), columns represent group where column one = 1, column
%   two = 0 (N.B. column one is expected to be minority class)
% trainSize = number of samples to put in training data for each subset
%   (e.g. per animal); format: integer
% testSize = number of samples to put in testing data for each subset
%   (e.g. per animal); format: integer
% method = imputation method to use in case of minority class being too
%   small; options: 'ADA' for ADASYN method (see ADASYN.mat) or 'gauss' for
%   manual gaussian distribution imputation
% rep = repetitions or number of times to split data (for creating many 
%   datasets to control for any potential effects of imputing and/or of 
%   partitioning); format: integer
%__________________________________________________________________________
% OUTPUTS:
% all = data structure with the following datasets corresponding to data
%   collapsed across all subsets (e.g. animals); all have format: 1 x rep 
%   cell array rep is number of repetitions
%   trainX = training predictors
%   trainY = training outcomes
%   testX = testing predictors
%   testY = testing outcomes
% each = data structure with the following datasets correpsonding to data
%   within each subset (e.g. animal); all have format: n x rep where  n is 
%   the number of subsets in original data and rep are repetitions
%   trainX = training predictors
%   trainY = training outcomes
%   testX = testing predictors
%   testY = testing outcomes
% rnd = randomized data structure; has format corresponding to mirror
%   dataset in either all or each structure (see above)
%   allTrainY = training predictors; randomized from all.trainY
%   allTestY = testing predictors; randomized from all.testY
%   eachTrainY = training predictors; randomized from each.trainY
%   eachTestY = testing predictors; randomized from each.testY
% w = weights; if not a 50/50% split, then calculates the necessary weights
%   for the minority class to be equal to the majority class; format:
%   cell array of weights with n columns for each subset (e.g. animal).
%   N.B. Sets majority class to have weight = 1.
%__________________________________________________________________________
% DEPENDENCIES:
% ADASYN.mat
%__________________________________________________________________________
% LLD 2017-2018
%% Initialize/Preallocate
eachTestX = cell(rep,size(data,1));
eachTestY = cell(rep,size(data,1));
xOneLeft = cell(rep,size(data,1));
xZeroLeft = cell(rep,size(data,1));
yOneLeft = cell(rep,size(data,1));
yZeroLeft = cell(rep,size(data,1));
newOneX = cell(rep,size(data,1));
newOneY = cell(rep,size(data,1));
eachTrainX = cell(rep,size(data,1));
eachTrainY = cell(rep,size(data,1));
eachTestYRand = cell(rep,size(data,1));
eachTrainYRand = cell(rep,size(data,1));
allTrainX = cell(1,rep); 
allTrainY = cell(1,rep);
allTestX = cell(1,rep);
allTestY = cell(1,rep);
allTrainYRand = cell(1,rep);
allTestYRand = cell(1,rep);
%% Split data
for r = 1:rep
    % Count number of ones and zeros in original dataset
    nOne = cell2mat(cellfun(@(x) size(x,1),data(:,1),'UniformOutput',0));
    nZero = cell2mat(cellfun(@(x) size(x,1),data(:,2),'UniformOutput',0));
    % Calculate distribution of ones to zeros
    d = nOne./(nOne+nZero);
    % Calculate number of ones and zeros to be put in training and test set
    % If whole numbers were given, then use 50/50% split for training;
    % replicate for indexing in for loop
    if rem(trainSize,1)== 0 && rem(testSize,1) == 0
        trainOneN = repmat(round(trainSize./2),1,size(data,1));
        trainZeroN = repmat(round(trainSize./2),1,size(data,1));
        testOneN = repmat(ceil(testSize.*d),1,size(data,1));
        testZeroN = repmat(floor(testSize.*(1-d)),1,size(data,1));
    % Otherwise use 'd' to maintain distribution in split given by decimals
    else 
       trainOneN = round(trainSize.*(nOne+nZero).*d);
       trainZeroN = round(trainSize.*(nOne+nZero).*(1-d));
       testOneN = round(testSize.*(nOne+nZero).*d);
       testZeroN = round(testSize.*(nOne+nZero).*(1-d));
    end
    % Calculate number of ones and zeros needed in test-set to replicate
    % distribution in original data; round ones up and zeros down
    
    for ii = 1:size(data,1)
        % Check if empty
        if ~isempty(data{ii,1})
            % Randomally generate indices of ones and zeros for test-sets
            thisTestOneInd = randperm(nOne(ii),testOneN(ii));
            thisTestZeroInd = randperm(nZero(ii),testZeroN(ii));
            % Pull ones and zeros found above into test-sets out so there
            % is no overlap between test and training sets
            eachTestX{r,ii} = [data{ii,1}(thisTestOneInd,:);...
                data{ii,2}(thisTestZeroInd,:)];
            eachTestY{r,ii} = [ones(testOneN(ii),1);...
                zeros(testZeroN(ii),1)];
            % Find leftover data indices
            xOneLeft{r,ii} = data{ii,1}(~ismember(1:nOne(ii),...
                thisTestOneInd),:);
            xZeroLeft{r,ii} = data{ii,2}(~ismember(1:nZero(ii),...
                thisTestZeroInd),:);
            % Pull leftover data out
            yOneLeft{r,ii} = ones(nOne(ii)-testOneN(ii),1);
            yZeroLeft{r,ii} = zeros(nZero(ii)-testZeroN(ii),1);
            % If imputing is desired, check if needed
            if ~isempty(method) && size(xOneLeft{r,ii},1) < trainSize/2
                if strcmpi(method,'ADA')
                    % Apply ADASYN to each set
                    [newOneX{r,ii},newOneY{r,ii}] = ADASYN([...
                        xOneLeft{r,ii};xZeroLeft{r,ii}],[yOneLeft{r,ii};...
                        yZeroLeft{r,ii}],1.5,5,5,0);
                    % Subsample new data s.t. new+old = trainSize/2
                    newOneX{r,ii}(trainSize/2-size(xOneLeft{r,ii},1)+1:...
                        end,:) = [];
                    newOneY{r,ii}(trainSize/2-size(xOneLeft{r,ii},1)+1:...
                        end,:) = [];
                elseif strcmpi(method,'gauss')
                    % Apply manual
                    n = size(xZeroLeft{r,ii},1)-size(xOneLeft{r,ii},1);
                    newOneX{r,ii} = normOversample(xOneLeft{r,ii},n);
                    newOneY{r,ii} = ones(n,1);
                end
                % Randomally generate indices of ones and zeros to be used
                % in training set
                thisTrainOneInd = randperm(size(xOneLeft{r,ii},1)+...
                    size(newOneX{r,ii},1),trainOneN(ii));
                thisTrainZeroInd = randperm(size(xZeroLeft{r,ii},1),...
                    trainZeroN(ii));
                % Combine leftover data with imputed data
                thisCatX = [xOneLeft{r,ii};newOneX{r,ii}];
                % If no imputing is needed, then just use leftover data
            else
                thisCatX = [xOneLeft{r,ii}];
                % Create indices from leftover data
                thisTrainOneInd = randperm(size(xOneLeft{r,ii},1),...
                    trainOneN(ii));
                thisTrainZeroInd = randperm(size(xZeroLeft{r,ii},1),...
                    trainZeroN(ii));
            end
            % Pull data out of combined dataset using indices found above
            eachTrainX{r,ii} = [thisCatX(thisTrainOneInd,:);...
                xZeroLeft{r,ii}(thisTrainZeroInd,:)];
            eachTrainY{r,ii} = [ones(trainOneN(ii),1);...
                zeros(trainZeroN(ii),1)];
            % Randomize test and training Ys
            eachTestYRand{r,ii} = eachTestY{r,ii}(randperm(testOneN(ii)+...
                testZeroN(ii),testOneN(ii)+testZeroN(ii)));
            eachTrainYRand{r,ii} = eachTrainY{r,ii}...
                (randperm(trainOneN(ii)+trainZeroN(ii),trainOneN(ii)+...
                trainZeroN(ii)));
        end
    end
    % Concatenate 'each' sets into 'all'
    allTrainX{r} = cat(1,eachTrainX{r,:});
    allTrainY{r} = cat(1,eachTrainY{r,:});
    allTestX{r} = cat(1,eachTestX{r,:});
    allTestY{r} = cat(1,eachTestY{r,:});
    allTestYRand{r} = cat(1,eachTestYRand{r,:});
    allTrainYRand{r} = cat(1,eachTrainYRand{r,:});
end
%% Calculate weights if not evenly split 50/50%
if rem(trainSize,1) ~= 0 && rem(testSize,1) ~= 0
    % Recalculate distribution since it may have changed slightly due to
    % rounding
    newD = cellfun(@(x) sum(x)/numel(x),eachTrainY(1,:))';
    % Calculate weights for ones (minority class)
    w = 1./newD-1;
    % Create weight array 
    for ii = 1:size(nOne,1)
       weights{ii} = [ones(trainOneN(ii),1).*w(ii);ones(trainZeroN(ii),1)]; 
    end
% OTherwise set 'weights' to be empty
else
    weights = [];
end
%% Create ouput data structures
% All
all.trainX = allTrainX;
all.trainY = allTrainY;
all.testX = allTestX;
all.testY = allTestY;
% Each
each.trainX = eachTrainX;
each.trainY = eachTrainY;
each.testX = eachTestX;
each.testY = eachTestY;
% Random
rnd.allTrainY = allTrainYRand;
rnd.allTestY = allTestYRand;
rnd.eachTrainY = eachTrainYRand;
rnd.eachTestY = eachTestYRand;
