load('/ihome/ldwiel/data/bingeNotBingeTrial.mat')
%%
for jj = 1:10
    dat = [];
    for ii = 1:12
        % Get distribution of ones to zeros
        nOne = size(data{1,1}{ii,1},1);
        nZero = size(data{1,1}{ii,2},1);
        dist = nOne/(nZero+nOne);
        % Get size of training and test sets
        full = nOne+nZero;
        trainN = ceil(full*0.8);
        testN = floor(full*0.2);
        % Randomize test set, maintaining original distribution
        testOne = round(testN*dist);
        testOneInds = randperm(nOne,testOne);
        testZero = testN-testOne;
        testZeroInds = randperm(nZero,testZero);
        testX{ii} = [data{1,1}{ii,1}(testOneInds,:);data{1,1}{ii,2}(testZeroInds,:)];
        testY{ii} = [ones(testOne,1);zeros(testZero,1)];
        % Find minimum train set (to have 20 ones)
        minT = ceil(20/dist);
        if minT < trainN
            % Randomize minimum train set from remaining ones and zeros
            trainOne = ceil(minT*dist);
            leftOne = logicFind(1,~ismember(1:nOne,testOneInds),'==');
            trainOneInds = leftOne(randperm(length(leftOne),trainOne));
            leftZero = logicFind(1,~ismember(1:nZero,testZeroInds),'==');
            trainZero = minT-trainOne;
            trainZeroInds = leftZero(randperm(length(leftZero),trainZero));
            nTrialData{ii,1} = [ones(size(trainOneInds,2),1),data{1,1}{ii,1}(trainOneInds,:);zeros(size(trainZeroInds,2),1),data{1,1}{ii,2}(trainZeroInds,:)];
            
            trainOneN = ceil(trainN*dist);
            % Find step needed for ~50 steps
            oneStep = round((trainOneN-trainOne)/50);
            oneUp = 1:oneStep:trainOneN-trainOne;
            % Find zero steps needed to keep distribution about the same
            zeroUp = round((1/dist-1).*oneUp);
            % Find remaining ones and zeros
            remainOne = logicFind(1,~ismember(1:nOne,[trainOneInds,testOneInds]),'==');
            remainZero = logicFind(1,~ismember(1:nZero,[trainZeroInds,testZeroInds]),'==');
            % Shuffle indices
            reaminOne = remainOne(randperm(length(remainOne)));
            reaminZero = remainZero(randperm(length(remainZero)));
            for k = 1:size(oneUp,2)-1
                nTrialData{ii,k+1} = [nTrialData{ii,k};ones(length(oneUp(k):oneUp(k+1)-1),1),data{1,1}{ii,1}(remainOne(oneUp(k):oneUp(k+1)-1),:);zeros(length(zeroUp(k):zeroUp(k+1)-1),1),data{1,1}{ii,2}(remainZero(zeroUp(k):zeroUp(k+1)-1),:)];
            end
        else
            nTrialData{ii,1} = [];
        end
        if sum(~cellfun(@isempty, nTrialData(ii,:))) > 1
            dat = [dat;nTrialData(ii,:)];
        else
            testX{ii} = [];
            testY{ii} = [];
        end
    end
    testX(cellfun(@isempty,testX)) = [];
    testY(cellfun(@isempty,testY)) = [];
    allData{jj} = dat;
    allTestX{jj} = testX;
    allTestY{jj} = testY;
end
%%