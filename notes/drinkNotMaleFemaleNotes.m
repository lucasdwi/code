%% Grab alcohol-not data from drinkNot
[data1,~,files1] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\drinkNot\'],{'.mat'},{'pow','coh'},'trl','rel');
%% Build datasets for lasso; even male/female split with ADASYN to impute 
% drinking
for jj = 1:100
    % Male
    maleDrink = []; maleNot = [];
    % AH12 = 6,7
    for ii = 6:7
        maleDrink = [maleDrink;...
            data1{1}{ii,1}(randperm(size(data1{1}{ii,1},1),9),:)];
        maleNot = [maleNot;...
            data1{1}{ii,2}(randperm(size(data1{1}{ii,2},1),25),:)];
    end
    % AH16 = 9,10,12
    allInds = [9,10,12];
    inds = allInds(randperm(3,2));
    for ii = inds
        maleDrink = [maleDrink;...
            data1{1}{ii,1}(randperm(size(data1{1}{ii,1},1),9),:)];
        maleNot = [maleNot;...
            data1{1}{ii,2}(randperm(size(data1{1}{ii,2},1),25),:)];
    end
    % Female
    femaleDrink = [];
    femaleNot = [];
    % First determine which females to build from
    inds = randperm(3,2);
    % AH3 = 27:30
    if any(inds==1)
        floof = 27:30;
        thisInd = randperm(4,2);
        for ii = floof(thisInd)
            femaleDrink = [femaleDrink;...
                data1{1}{ii,1}(randperm(size(data1{1}{ii,1},1),9),:)];
            femaleNot = [femaleNot;...
                data1{1}{ii,2}(randperm(size(data1{1}{ii,2},1),25),:)];
        end
    end
    % AH4 = 31:34
    if any(inds==2)
        floof = 31:34;
        thisInd = randperm(4,2);
        for ii = floof(thisInd)
            femaleDrink = [femaleDrink;...
                data1{1}{ii,1}(randperm(size(data1{1}{ii,1},1),9),:)];
            femaleNot = [femaleNot;...
                data1{1}{ii,2}(randperm(size(data1{1}{ii,2},1),25),:)];
        end
    end
    % AH5 = 36:40
    if any(inds==3)
        floof = 36:40;
        thisInd = randperm(4,2);
        for ii = floof(thisInd)
            femaleDrink = [femaleDrink;...
                data1{1}{ii,1}(randperm(size(data1{1}{ii,1},1),9),:)];
            femaleNot = [femaleNot;...
                data1{1}{ii,2}(randperm(size(data1{1}{ii,2},1),25),:)];
        end
    end
    % Impute drinking of both sexes up with ADASYN
    newMaleX = ADASYN([maleDrink;maleNot],[ones(36,1);...
        zeros(100,1)],1.2,5,5,0);
    maleDrink = [maleDrink;newMaleX(randperm(size(newMaleX,1),64),:)];
    newFemaleX = ADASYN([femaleDrink;femaleNot],[ones(36,1);...
        zeros(100,1)],1.2,5,5,0);
    femaleDrink = [femaleDrink;newFemaleX(randperm(size(newFemaleX,1),...
        64),:)];
    % Build datasets
    xMaleData(:,:,jj) = [maleDrink;maleNot];
    yMaleData(:,jj) = [ones(100,1);zeros(100,1)];
    xFemaleData(:,:,jj) = [femaleDrink;femaleNot];
    yFemaleData(:,jj) = [ones(100,1);zeros(100,1)];
end
%% Save above data
save(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'maleFemaleDrinkNotData.mat'],'xMaleData','yMaleData','xFemaleData',...
    'yFemaleData')
%% Run single feature models
load(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\'...
    'maleFemaleDrinkNotData.mat'])
for ii = 1:60
    disp(num2str(ii))
    for jj = 1:100
        testInds = randperm(200,20);
        trainInds = logicFind(1,~ismember(1:200,testInds),'==');
        % Male
        maleMdl = fitglm(xMaleData(trainInds,ii,jj),...
            yMaleData(trainInds,jj),'distribution','binomial',...
            'binomialsize',180);
        prob = predict(maleMdl,xMaleData(testInds,ii,jj));
        [mX(ii,jj,:),mY(ii,jj,:),~,mA(ii,jj)] = perfcurve(...
            yMaleData(testInds,jj),prob,1);
        % Female
        femaleMdl = fitglm(xFemaleData(trainInds,ii,jj),...
            yFemaleData(trainInds,jj),'distribution','binomial',...
            'binomialsize',180);
        prob = predict(femaleMdl,xFemaleData(testInds,ii,jj));
        [fX(ii,jj,:),fY(ii,jj,:),~,fA(ii,jj)] = perfcurve(...
            yFemaleData(testInds,jj),prob,1);
    end
end
%%
[mAS,mInd] = sort(mean(mA,2),'descend');
maleSortA = mA(mInd,:);
[fAS,fInd] = sort(mean(fA,2),'descend');
femaleSortA = fA(fInd,:);
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
mName = nameVect(mInd)';
fName = nameVect(fInd)';
%% Load lasso results
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper3\analyzed\',...
    'maleFemaleDrinkNot\'])
mErr = []; mRndErr = []; fErr = []; fRndErr = [];
for ii = 1:100
    load(['maleFemaleDrinkNot',num2str(ii),'.mat'])
    mErr = [mErr;male.allLambda{1}.allErr];
    mRndErr = [mRndErr;rndMale.allLambda{1}.allErr];
    fErr = [fErr;female.allLambda{1}.allErr];
    fRndErr = [fRndErr;rndFemale.allLambda{1}.allErr];
end
doubleHist((1-mErr).*100,(1-mRndErr).*100,'xlab','Accuracy (%)',...
    'main','Drink vs. Not: Male')
doubleHist((1-fErr).*100,(1-fRndErr).*100,'xlab','Accuracy (%)',...
    'main','Drink vs. Not: Female')
%%
n1 = cell2mat(cellfun(@(x) size(x,1),data1{1,1},'uniformoutput',0));
animal1 = []; day1 = [];
for ii = 1:size(data1{1},1)
    for jj = 1:size(data1{1},2)
        parts = strsplit(files1{1}{ii},'_');
        animal1{ii,jj} = repmat({parts{1}},n1(ii,jj),1);
        day1{ii,jj} = repmat({parts{2}},n1(ii,jj),1);
    end
end
dataTable = array2table([cat(1,data1{1}{:,1});cat(1,data1{1}{:,2})]);
% Grab alcohol and not data from paper3
% First get drinking
[data2,~,files2] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\waterAlcohol\processed\'],{'.mat'},{'pow','coh'},'trl',...
    'rel');
animal2 = []; day2 = [];
n2 = cell2mat(cellfun(@(x) size(x,1),data2{1,1}(:,2),'uniformoutput',0));
for ii = 1:size(data2{1},1)
    parts = strsplit(files2{1}{ii},'_');
    animal2{ii} = repmat({[parts{1},parts{2}]},n2(ii),1);
    day2{ii} = repmat({parts{4}},n2(ii),1);
end
% Add to dataTable
dataTable = [dataTable;array2table(cat(1,data2{1}{:,2}))];
% Then get not drinking
[data3,~,files3] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\waterAlcohol\notDrink\'],{'.mat'},{'pow','coh'},'trl',...
    'rel');
% Remove the three files that don't have equivalents in data2 (drinking)
data3{1} = data3{1}([1:14,16:18,20:26,28:end],:);
files3{1} = files3{1}([1:14,16:18,20:26,28:end]);
n3 = cell2mat(cellfun(@(x) size(x,1),data3{1,1},'uniformoutput',0));
animal3 = []; day3 = [];
for ii = 1:size(data3{1},1)
    parts = strsplit(files3{1}{ii},'_');
    animal3{ii} = repmat({[parts{1},parts{2}]},n3(ii),1);
    day3{ii} = repmat({parts{4}},n3(ii),1);
end
% Add to dataTable
dataTable = [dataTable;array2table(cat(1,data3{1}{:}))];


%% Add other variables
% Add animal and day
allAnimal = [cat(1,animal1{:,1});cat(1,animal1{:,2});cat(1,animal2{:});...
    cat(1,animal3{:})];
allDay = [cat(1,day1{:,1});cat(1,day1{:,2});cat(1,day2{:});...
    cat(1,day3{:})];
dataTable.ID = allAnimal;
dataTable.day = allDay;
% Add group
dataTable.group = [ones(sum(n1(:,1)),1);zeros(sum(n1(:,2)),1);...
    ones(sum(n2),1);zeros(sum(n3),1)];
%% Set up sex variable
m = {'AH11','AH12','AH16','AH17','AH20','AH38','AH39','AH87','AH88',...
    'AH90'};
f = {'AH1','AH3','AH4','AH5','AH8','AH26','AH29','AH30','AH80','AH81'};
for ii = 1:size(m,2)
    dataTable.sex(strcmp(m{ii},dataTable.ID)) = {'m'};
end
for ii = 1:size(f,2)
    dataTable.sex(strcmp(f{ii},dataTable.ID)) = {'f'};
end
%%
