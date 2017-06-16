% Uses ADASYN method to ensure 50-50 split in training set
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\bingeNotBingeTrial.mat','data')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\allbingeNotData.mat','data')
%%
di = 1;
for r = 1:20
    nOne = cell2mat(cellfun(@(x) size(x,1),data{1,di}(:,1),'UniformOutput',0));
    nZero = cell2mat(cellfun(@(x) size(x,1),data{1,di}(:,2),'UniformOutput',0));
    d = nOne./(nOne+nZero);
    trainSize = 500;
%     trainOneN = round(trainSize.*d);
%     trainZeroN = round(trainSize.*(1-d));
    testSize = .2*(trainSize/.8);
    testOneN = round(testSize.*d);
    testZeroN = round(testSize.*(1-d));
    for ii = 1:size(data{1,di},1)
        % Separate each animals test-set
        thisTestOneInd = randperm(nOne(ii),testOneN(ii));
        thisTestZeroInd = randperm(nZero(ii),testZeroN(ii));
        eachTestX{r,ii} = [data{1,di}{ii,1}(thisTestOneInd,:);data{1,di}{ii,2}(thisTestZeroInd,:)];
        eachTestY{r,ii} = [ones(testOneN(ii),1);zeros(testZeroN(ii),1)];
        % Randomize
%         eachTestY{r,ii} = eachTestY{r,ii}(randperm(testSize,testSize));
        % Get leftover data
        xOneLeft{r,ii} = data{1,di}{ii,1}(~ismember(1:nOne(ii),thisTestOneInd),:);
        xZeroLeft{r,ii} = data{1,di}{ii,2}(~ismember(1:nZero(ii),thisTestZeroInd),:);
        yOneLeft{r,ii} = ones(nOne(ii)-testOneN(ii),1);
        yZeroLeft{r,ii} = zeros(nZero(ii)-testZeroN(ii),1);
        % Apply ADASYN to each set
        [newOneX{r,ii},newOneY{r,ii}] = ADASYN([xOneLeft{r,ii};xZeroLeft{r,ii}],[yOneLeft{r,ii};yZeroLeft{r,ii}],1,5,5,0);
        % Randomly extract n (trainSize*.50) samples from combined new and old ones
        trainSize = 46;
        thisTrainOneInd = randperm(size(xOneLeft{r,ii},1)+size(newOneX{r,ii},1),round(trainSize/2));
        thisTrainZeroInd = randperm(size(xZeroLeft{r,ii},1),round(trainSize/2));
        thisCatX = [xOneLeft{r,ii};newOneX{r,ii}];
        eachTrainX{r,ii} = [thisCatX(thisTrainOneInd,:);xZeroLeft{r,ii}(thisTrainZeroInd,:)];
        eachTrainY{r,ii} = [ones(round(trainSize/2),1);zeros(round(trainSize/2),1)];
        % Randomize
%         eachTrainY{r,ii} = eachTrainY{r,ii}(randperm(trainSize,trainSize));
    end
    allTrainX{r} = cat(1,eachTrainX{r,:});
    allTrainY{r} = cat(1,eachTrainY{r,:});
    allTestX{r} = cat(1,eachTestX{r,:});
    allTestY{r} = cat(1,eachTestY{r,:});
end
%% Prebinge data
% Collate data
[data,samp] = collateData('C:\Users\Lucas\Desktop\GreenLab\data\paper2\preBingeCombined2\',{'base'},{'pow','coh'},'trl');
%% Split into training and testing sets
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\preBingeRaw.mat')
for di = 1:size(data{1,1},2)-1
    for r = 1:20
        for ii = 1:size(data{1,1},1)
            if ~isempty(data{1,1}{ii,di})
                nOne = size(data{1,1}{ii,di},1);
                if nOne > 1
                    nZero = size(data{1,1}{ii,end},1);
                    dist = nOne/(nOne+nZero);
                    nTestOne = ceil(nOne/2);
                    nTestZero = round((nTestOne*(1-dist))/dist);
                    nTrainOne = nOne-nTestOne;
                    nTrainZero = nTrainOne;
                    % Get all one and zero inds - without replacement
                    allOne = randsample(nOne,nTestOne+nTrainOne);
                    allZero = randsample(nZero,nTestZero+nTrainZero);
                    % Use indices to split into sets
                    eachTrainX{r,ii} = [data{1,1}{ii,di}(allOne(1:nTrainOne),:);data{1,1}{ii,end}(allZero(1:nTrainZero),:)];
                    eachTrainY{r,ii} = [ones(nTrainOne,1);zeros(nTrainZero,1)];
                    eachTestX{r,ii} = [data{1,1}{ii,di}(allOne(nTrainOne+1:end),:);data{1,1}{ii,end}(allZero(nTrainZero+1:end),:)];
                    eachTestY{r,ii} = [ones(nTestOne,1);zeros(nTestZero,1)];
                    % Randomize
                    eachTrainYRand{r,ii} = eachTrainY{r,ii}(randperm(nTrainOne+nTrainZero,nTrainOne+nTrainZero));
                    eachTestYRand{r,ii} = eachTestY{r,ii}(randperm(nTestOne+nTestZero,nTestOne+nTestZero));
                end
            end
        end
        allTrainX{r,di} = cat(1,eachTrainX{r,:});
        allTrainY{r,di} = cat(1,eachTrainY{r,:});
        allTestX{r,di} = cat(1,eachTestX{r,:});
        allTestY{r,di} = cat(1,eachTestY{r,:});
        
        allTrainYRand{r,di} = cat(1,eachTrainYRand{r,:});
        allTestYRand{r,di} = cat(1,eachTestYRand{r,:});
        eachTrainX = []; eachTrainY = []; eachTestX = []; eachTestY = []; eachTrainYRand = []; eachTestYRand = [];
    end
end
%% Balance 'all' training sets by subsampling 'each'
num = size(eachTrainX{1,1},1);
per = ceil(num/size(eachTrainX,2));
for ii = 1:20
    thisX = [];
    thisY = [];
    for jj = 1:11
        x = eachTrainX{ii,jj}([randperm(250,per/2);randperm(250,per/2)+250],:);
        y = [ones(21,1);zeros(21,1)];
        thisX = [thisX;x];
        thisY = [thisY;y];
    end
    balancedTrainX{ii} = thisX;
    balancedTrainY{ii} = thisY;
end
%% Create increasing training sets from 40 to 700
n = 40:20:trainSize;
for ni = 1:length(n)
    for r = 1:20
        for ii = 1:12
            thisTrainX{r,ii} = [eachTrainX{r,ii}(1:n(ni)/2,:);eachTrainX{r,ii}(trainSize/2+1:trainSize/2+n(ni)/2,:)];
            thisTrainY{r,ii} = [eachTrainY{r,ii}(1:n(ni)/2,:);eachTrainY{r,ii}(trainSize/2+1:trainSize/2+n(ni)/2,:)];
        end
    end
    trainX{ni} = thisTrainX;
    trainY{ni} = thisTrainY;
end
%% Rand dist and imputer (manual bootstrapping)
x = data{1,1}{1,1};
for ii = 1:60
   m = mean(x(:,ii));
   s = std(x(:,ii));
   pd = makedist('Normal','mu',m,'sigma',s);
   r(:,ii) = random(pd,1,598);
end
%% Splits data into equally sized train and test sets with original
% distributions maintained
for jj = 1:20
    for ii = 1:size(data,2)
        for k = 1:size(data{1,ii},1)
            nOne(ii,k) = size(data{1,ii}{k,1},1);
            nZero(ii,k) = size(data{1,ii}{k,2},1);
        end
    end
    nFull = nOne+nZero;
    dist = nOne./nFull;
    % Grab minTrain trials from each animal to populate training set
    minTrain = 300;
    trainOne = round(dist.*minTrain);
    trainZero = round((1-dist).*minTrain);
    % Get number of trials left
    oneLeft = nOne-trainOne;
    zeroLeft = nZero-trainZero;
    allLeft = oneLeft+zeroLeft;
    % Determine smallest size to keep distribution
    minTest = 75;%min(allLeft,[],2);
    testOne = round(dist.*minTest);
    testZero = round((1-dist).*minTest);
    % Baseline
    for k = 1:12
        oneIndsTrain{k} = randperm(nOne(1,k),trainOne(1,k));
        zeroIndsTrain{k} = randperm(nZero(1,k),trainZero(1,k));
        eachTrainX{jj,k} = [data{1,1}{k,1}(oneIndsTrain{k},:);data{1,1}{k,2}(zeroIndsTrain{k},:)];
        eachTrainY{jj,k} = [ones(trainOne(1,k),1);zeros(trainZero(1,k),1)];
    end
    allTrainX(:,:,jj) = cat(1,eachTrainX{jj,:});
    allTrainY(:,:,jj) = cat(1,eachTrainY{jj,:});
    for k = 1:12
        % Pick out n (testOne) indices that were not in the training set
        oneIndsVect = 1:nOne(1,k);
        oneIndsLeft = oneIndsVect(~ismember(oneIndsVect,oneIndsTrain{k}));
        oneIndsTest{k} = oneIndsLeft(randperm(oneLeft(1,k),testOne(1,k)));
        % Pick out n (testZero) indices that were not in the training set
        zeroIndsVect = 1:nZero(1,k);
        zeroIndsLeft = zeroIndsVect(~ismember(zeroIndsVect,zeroIndsTrain{k}));
        zeroIndsTest{k} = zeroIndsLeft(randperm(zeroLeft(1,k),testZero(1,k)));
        eachTestX{jj,k} = [data{1,1}{k,1}(oneIndsTest{k},:);data{1,1}{k,2}(zeroIndsTest{k},:)];
        eachTestY{jj,k} = [ones(testOne(1,k),1);zeros(testZero(1,k),1)];
    end
    allTestX(:,:,jj) = cat(1,eachTestX{jj,:});
    allTestY(:,:,jj) = cat(1,eachTestY{jj,:});
end
%% Sets up 80% train sets and equally sized test sets <= 20% using
% smallest test set
for jj = 1:20
    for ii = 1:size(data,2)
        for k = 1:size(data{1,ii},1)
            nOne(ii,k) = size(data{1,ii}{k,1},1);
            nZero(ii,k) = size(data{1,ii}{k,2},1);
        end
    end
    nFull = nOne+nZero;
    dist = nOne./nFull;
    fullTwentySplit = round(nFull.*.2);
    % Replace 0s with NaN
    fullTwentySplit(fullTwentySplit == 0) = NaN;
    % Find min test set size
    minTwentySplit = min(fullTwentySplit,[],2);
    % Get number of ones and zeros using minTwentySplit
    testMinOne = round(minTwentySplit(1,1).*dist);
    testMinZero = round(minTwentySplit(1,1).*(1-dist));
    % Use min test set size
    for k = 1:12
        oneIndsTest{k} = randperm(nOne(1,k),testMinOne(1,k));
        zeroIndsTest{k} = randperm(nZero(1,k),testMinZero(1,k));
        eachMinTestX{jj,k} = [data{1,1}{k,1}(oneIndsTest{k},:);data{1,1}{k,2}(zeroIndsTest{k},:)];
        eachMinTestY{jj,k} = [ones(testMinOne(1,k),1);zeros(testMinZero(1,k),1)];
    end
    allMinTestX{jj} = cat(1,eachMinTestX{jj,:});
    allMinTestY{jj} = cat(1,eachMinTestY{jj,:});
    oneLeft = nOne-testMinOne;
    zeroLeft = nZero-testMinZero;
    % Propagate min size to training set size
    minEightySplit = minTwentySplit.*4;
    trainOne = round(minEightySplit.*dist);
    trainZero = round(minEightySplit.*(1-dist));
    for k = 1:12
        oneIndsVect = 1:nOne(1,k);
        oneIndsLeft = oneIndsVect(~ismember(oneIndsVect,oneIndsTest{k}));
        oneIndsTrain{k} = oneIndsLeft(randperm(oneLeft(1,k),trainOne(1,k)));
        zeroIndsVect = 1:nZero(1,k);
        zeroIndsLeft = zeroIndsVect(~ismember(zeroIndsVect,zeroIndsTest{k}));
        zeroIndsTrain{k} = zeroIndsLeft(randperm(zeroLeft(1,k),trainZero(1,k)));
        eachMinTrainX{jj,k} = [data{1,1}{k,1}(oneIndsTrain{k},:);data{1,1}{k,2}(zeroIndsTrain{k},:)];
        eachMinTrainY{jj,k} = [ones(trainOne(1,k),1);zeros(trainZero(1,k),1)];
    end
end
%%
nOne = size(data{1,1}{1,1},1);
nZero = size(data{1,1}{1,2},1);
nFull = nOne+nZero;
dist = nOne./nFull;
minTrain = 300;
trainOne = round(dist.*minTrain);
trainZero = round((1-dist).*minTrain);
oneLeft = nOne-trainOne;
zeroLeft = nZero-trainZero;
allLeft = oneLeft+zeroLeft;
oneIndsTrain = randperm(nOne,trainOne);
zeroIndsTrain = randperm(nZero,trainZero);
eachTrainX = [data{1,1}{1,1}(oneIndsTrain,:);data{1,1}{1,2}(zeroIndsTrain,:)];
eachTrainY = [ones(trainOne,1);zeros(trainZero,1)];
cfg = lassoNetCfg([],[],'n','y','n',100,'1se');
[~,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(eachTrainX,eachTrainY,'binomial','deviance',1,10,1,cfg);
model = cvFitsArray{1,1}{allLambda{1,1}.bestLambdaInds,1};
%%
minTest = 75;
for ii = 1:20
    c = 1;
    for d = 0.05:0.01:0.8
        testOne = round(d.*minTest);
        testZero = round((1-d).*minTest);
        oneIndsVect = 1:nOne;
        oneIndsLeft = oneIndsVect(~ismember(oneIndsVect,oneIndsTrain));
        oneIndsTest = oneIndsLeft(randperm(oneLeft,testOne));
        zeroIndsVect = 1:nZero;
        zeroIndsLeft = zeroIndsVect(~ismember(zeroIndsVect,zeroIndsTrain));
        zeroIndsTest = zeroIndsLeft(randperm(zeroLeft,testZero));
        eachTestX{c,ii} = [data{1,1}{1,1}(oneIndsTest,:);data{1,1}{1,2}(zeroIndsTest,:)];
        eachTestY{c,ii} = [ones(testOne,1);zeros(testZero,1)];
        [predY] = cvglmnetPredict(model,eachTestX{c,ii},'lambda_1se','response');
        [~,~,~,eachA(c,ii)] = perfcurve(eachTestY{c,ii},predY,1);
        c = c+1;
    end
end
%%
dists = cellfun(@sum, eachTestY)./75;
figure
plot(dists(:,1),mean(eachA,2),'o')
figure
hold on
for ii = 1:20
   plot(dists(:,ii),eachA(:,ii),'ok')
end
%%
prob = rand(10000,1);
y = round(abs(prob-rand(10000,1)./2));
[rocx,rocy,~,a] = perfcurve(y,prob,1);
figure
plot(rocx,rocy)
hold on
for ii = 1:100
    inds = randperm(10000,100);
    thisProb = prob(inds);
    thisy = y(inds);
    [thisrocx{ii},thisrocy{ii},~,thisa(ii)] = perfcurve(thisy,thisProb,1);
    plot(thisrocx{ii},thisrocy{ii},'r-')
end
