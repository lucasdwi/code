n = 100;
for ii = 1:n
%    load([num2str(ii),'.mat'],'cvFitsArray')
   cvm(:,ii) = cvFitsArray{1,1}{ii,1}.cvm;
%    acc(ii) = a;
%    mInd = allLambda{1,1}.bestLambdaInds;
%    lamInd = logicFind(allLambda{1,1}.bestLambda,cvFitsArray{1,1}{allLambda{1,1}.bestLambdaInds,1}.lambda,'==')-10;
%    beta(:,ii) = cvFitsArray{1,1}{mInd,1}.glmnet_fit.beta(:,lamInd);
   lam(:,ii) = cvFitsArray{1,1}{ii,1}.lambda;
end
%%
for k = 1:20
    load([num2str(k),'.mat'])
    for ii = 1:100
        mse(:,ii) = cvFitsArray{1,1}{ii,1}.cvm;
        mseM(ii) = mean(mse(:,ii),1);
        lam(ii) = cvFitsArray{1,1}{ii,1}.lambda_1se;
    end
    [~,minMSEind] = min(mseM);
    model = cvFitsArray{1,1}{minMSEind,1};
    [predY] = cvglmnetPredict(model,testX,'lambda_1se','response');
    [~,~,~,a2(k)] = perfcurve(testY,predY,1);
    b2(k,:) = model.glmnet_fit.beta(:,logicFind(model.lambda_1se,model.lambda,'=='))';
    l2(k) = model.lambda_1se;
    a1(k) = a;
    b1(k,:) = allBeta{1,1}.betas(allLambda{1,1}.bestLambdaInds,:);
    l1(k) = allLambda{1,1}.bestLambda;
    clearvars -except k a1 a2 b1 b2 l1 l2
end
nB1 = sum(b1~=0,2);
nB2 = sum(b2~=0,2);
%%
for ii = 1:20
    load([num2str(ii),'.mat'])
    % Store x, y, a, and number of Beta for LASSO
    xLasso{ii} = x; yLasso{ii} = y; aLasso(ii) = a;
    nB(ii) = sum(allBeta{1,1}.betas(allLambda{1,1}.bestLambdaInds,:)~=0);
    [md1] = fitglm(trainX,trainY,'distribution','binomial');
    [prob] = predict(md1,testX);
    [xLog{ii},yLog{ii},~,aLog(ii)] = perfcurve(testY,prob,1);
    clearvars -except xLog yLog aLog xLasso yLasso aLasso nB
end
%%
% Run PCA on training data
[coeff,scores,latent,~,explained] = pca(trainX,'Economy',false,'Centered',true);
% Calculate cumulative percent of variance explained
percVar = cumsum(explained);
% Find number of PCs needed to explain 80% of the variance
ind = logicFind(80,percVar,'>=','first');
pcVar = percVar(ind);
% Project test data into PCA space - get scores
testScores = testX/coeff';
[~,~,~,~,~,~,probTest] = logPredict(scores(:,1:ind),trainY,testScores(:,1:ind),testY);
[pcaX,pcaY,~,a] = perfcurve(trainY,probTest,1);
% probAllTest = [probAllTest;[probTest{1,1},probTest{1,2}]];
[~,~,~,~,~,~,probTrain{iO}] = logPredict(scores(:,1:ind(ii)),trainY,[],[]);
%% Plot all ROCs
figure
hold on
for ii = 1:20
    plot(xLasso{ii},yLasso{ii},'k')
end
figure
hold on
for ii = 1:20
    plot(xLog{ii},yLog{ii},'r')
end
%% Manual probability caclulation
% coef = repmat(table2array(md1.Coefficients(:,1))',size(testX,1),1);
% logOdds = sum([ones(size(testX,1),1),testX].*coef,2);
% odds = exp(logOdds);
% prob = odds./(1+odds);
%% Check nTrial data against self prediction
for ii = 1:240
   load([num2str(ii),'.mat'],'accArray')
   a(ii) = accArray{1,1}.acc;
end
a = reshape(a,20,12);
aM = mean(a,1);
aS = std(a,[],1);
load('bingeNotBingeTrial.mat', 'trlDat')
for ii = 1:12
    nT(ii) = size(trlDat{1,ii},1);
    r(ii) = (sum(trlDat{1,ii}(:,1)==1)/nT(ii)).*100;
end
% Plot # trials vs. AUC
figure
scatter(nT,aM,'ko','filled')
hold on
for ii = 1:12
   plot([nT(ii),nT(ii)],[aM(ii)-aS(ii),aM(ii)+aS(ii)],'k') 
end
title('AUC vs. Number of Trials')
xlabel('Number of Trials')
ylabel('AUC (\mu \pm \sigma)')
% Plot % binge vs. AUC
figure
scatter(r,aM,'ko','filled')
hold on
for ii = 1:12
   plot([r(ii),r(ii)],[aM(ii)-aS(ii),aM(ii)+aS(ii)],'k') 
end
title('AUC vs. % Binge')
xlabel('% Binge')
ylabel('AUC (\mu \pm \sigma)')