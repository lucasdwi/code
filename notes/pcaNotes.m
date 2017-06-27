load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\paper1data.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\indGroups.mat')
%%
y = y(:,[1,3]);
normX = x;
inds{1} = 1:60;
%%
for ii = 1:2
    for iO = 1:2:24
        trainX = normX(~ismember(1:24,iO:iO+1),:);
        trainY = y(~ismember(1:24,iO:iO+1),ii);
        testX = normX(iO:iO+1,:);
        testY = y(iO:iO+1,ii);
        % Run PCA on training data
        [coeff,scores,latent,~,explained] = pca(trainX,'Economy',false,'Centered',1);
        % Calculate cumulative percent of variance explained
        percVar = cumsum(explained);
        % Find number of PCs needed to explain 80% of the variance
        ind(iO) = logicFind(80,percVar,'>=','first');
        pcVar(iO) = percVar(ind(iO));
        % Project test data into PCA space - get scores
        testScores = testX/coeff';
        mdl = fitglm(scores(:,1:ind(iO)),trainY,'distribution','binomial');
        prob(ii,iO:iO+1) = predict(mdl,testScores(:,1:ind(iO)));
    end
    [~,~,~,a(ii)] = perfcurve(y(:,1),prob(ii,:)',1);
end
%%
probAllTest = []; probTrain = [];
% Cycle through each subset in inds/indName
for ii = 1%:size(inds,1)
    % Cycle through and take out each pair of two files from the same
    % animal
    for iO = 1:2:24       
        trainX = normX(~ismember(1:24,iO:iO+1),:);
        trainY = y(~ismember(1:24,iO:iO+1),:);
        testX = normX(iO:iO+1,:);
        testY = y(iO:iO+1,:);
        % Run PCA on training data
        [coeff,scores,latent,~,explained] = pca(trainX(:,inds{ii}),'Economy',false,'Centered',false);
        % Calculate cumulative percent of variance explained
        percVar = cumsum(explained);
        % Find number of PCs needed to explain 80% of the variance
        ind(ii) = logicFind(80,percVar,'>=','first');
        pcVar(ii) = percVar(ind(ii));
        % Project test data into PCA space - get scores
        testScores = testX/coeff';
        % Use those PCs in a logistic regression, with test set
        [~,~,~,~,~,~,probTest] = logPredict(scores(:,1:ind(ii)),trainY,testScores(:,1:ind(ii)),testY);
        probAllTest = [probAllTest;[probTest{1,1},probTest{1,2}]];
%         [~,~,~,~,~,~,probTrain{iO}] = logPredict(scores(:,1:ind(ii)),trainY,[],[]);
%         probAllTrain = [probAllTrain;[probTrain(1,1),probTrain{1,2}]];
        % Without test set
%         [AUC(ii,:,iO),xMat{ii,iO},yMat{ii,iO},dev{ii,iO},~] = logPredict(scores(:,1:ind(ii)),y,[],[]);
    end
end
for ii = 1:2
   [~,~,~,a(ii)] = perfcurve(y(:,ii),probAllTest(:,ii),1); 
end
%%
for ii = 1:size(inds,1)
%     for ri = 1:100
%         for ni = 1:size(normX,1)
%            thisPerm(ni,:) = randperm(size(normX,2));
%         end
% Randomized
%     [coeff,scores,latent,~,explained] = pca(normX(thisPerm(:,inds{ii})),'Economy',false);
% Go through different inds
    [coeff,scores,latent,~,explained] = pca(normX(:,inds{ii}),'Economy',false,'Centered',false);
% Use all x
%     [coeff,scores,latent,~,explained] = pca(normX,'Economy',false);
    % Plot number of prinicple components vs. % of variance explained
    percVar = cumsum(explained);
%     figure;
%     stairs(percVar)
%     title('Variance Explained by {\itn} Princinple Components')
%     xlabel('{\itn} Principal Components'); ylabel('% of Variance Explained')
    % Use first n pcs for regression
    ind(ii) = logicFind(80,percVar,'>=','first');
    pcVar(ii) = percVar(ind(ii));
    [AUC(ii,:),xMat{ii},yMat{ii},dev{ii},~,~,prob{ii}] = logPredict(scores(:,1:ind(ii)),y,[],[]);
%     [AUC(ii,:,ri),xMat{ii,ri},yMat{ii,ri},dev{ii,ri}] = logPredict(scores(:,1:ind(ii)),y);
%     for ri = 1:100
%         thisPerm = randperm(size(y,1));
%         yRand = y(thisPerm',:);
%         [randAUC(ii,:,ri),xMatRand{ii,ri},yMatRand{ii,ri},devRand{ii,ri}] = logPredict(scores(:,1:ind(ii)),yRand);
%     end
%     end
end
%% Get random AUCs and ROCs
aucRandMean = mean(randAUC,3);
aucRandSd = std(randAUC,[],3);
for xi = 1:size(xMatRand,1)
    for yi = 1:size(xMatRand,2)
        for zi = 1:2
            if size(xMatRand{xi,yi}{1,zi},1) < 25
                catRandX{zi}(xi,yi,:) = cell2mat([xMatRand{xi,yi}(1,zi);ones(25-size(xMatRand{xi,yi}{1,zi},1),1)]);
                catRandY{zi}(xi,yi,:) = cell2mat([yMatRand{xi,yi}(1,zi);ones(25-size(yMatRand{xi,yi}{1,zi},1),1)]);
            else
                catRandX{zi}(xi,yi,:) = cell2mat(xMatRand{xi,yi}(1,zi));
                catRandY{zi}(xi,yi,:) = cell2mat(yMatRand{xi,yi}(1,zi));
            end
        end
    end
end
xAvgRand = sq(mean(catRandX{2},2));
yAvgRand = sq(mean(catRandY{2},2));
%% Set up real AUC matrix
% aucComp = [AUC(1,:);AUC(2,:);AUC(3,:)];
% devComp = [dev{1,1};dev{1,2};dev{1,3}];
%% Create real AUC table
aucTab = table(AUC(:,1),AUC(:,2),AUC(:,3),AUC(:,4),'RowNames',indName,'VariableNames',{'SA','SS','CA','CS'});