testX = []; testY = []; trainX = []; trainY = [];
for ii = 1:size(data{1,1},1)
    nOne = size(data{1,1}{ii,1},1);
    nZero = size(data{1,1}{ii,2},1);
    dist = nOne/(nZero+nOne);
    nNaive = ceil((nOne+nZero)*.20);
    naiveOne = ceil(nNaive*dist);
    rng('shuffle')
    oneInd{ii} = randperm(nOne,naiveOne);
    naiveZero = nNaive-naiveOne;
    rng('shuffle')
    zeroInd{ii} = randperm(nZero,naiveZero);
    testX = [testX;data{1,1}{ii,1}(oneInd{ii},:);data{1,1}{ii,2}(zeroInd{ii},:)];
    testY = [testY;ones(naiveOne,1);zeros(naiveZero,1)];
    trainX = [trainX;data{1,1}{ii,1}(~ismember(1:nOne,oneInd{ii}),:);data{1,1}{ii,2}(~ismember(1:nZero,zeroInd{ii}),:)];
    trainY = [trainY;ones(nOne-naiveOne,1);zeros(nZero-naiveZero,1)];
end
cfg = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(trainX,trainY,'binomial','auc',1,5,1,cfg);

