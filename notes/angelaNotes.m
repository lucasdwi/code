%% Angela Data
naccMdl = fitglm(drink,naccDrink,'distribution','binomial');
pfcMdl = fitglm(drink(~isnan(pfcDrink)),pfcDrink(~isnan(pfcDrink)),'distribution','binomial','link','logit');
y = glmval(table2array(naccMdl.Coefficients(:,1)),drink,'logit');
naccLog = sortrows([y,drink],1);
figure
plot(drink,y)
scatter(drink,nacc)
%% Reorder channels from Nick's animals (N<20)to be the same as Angela's (N>=20)
LFPTs.data = flipud(LFPTs.data);
LFPTs.label = fliplr(LFPTs.label);
%% NAcc Rest - Increase and Decrease
cd('C:\Users\Pythia\Documents\GreenLab\data\angela\')
load('NAccDecrease.mat')
realAcc = (1-allData.allLambda{1,1}.allErr).*100;
load('NAccDecreaseRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end
doubleHist(realAcc,randAcc,'Unit','%','Main','NAcc Rest: Decrease','Loc','northeast','xlab','Accuracy (%)');

load('NAccIncrease.mat')
realAcc = (1-allData.allLambda{1,1}.allErr).*100;
load('NAccIncreaseRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end
doubleHist(realAcc,randAcc,'Unit','%','Main','NAcc Rest: Increase','Loc','northeast','xlab','Accuracy (%)');
%% NAcc All - Increase and Decrease
cd('C:\Users\Pythia\Documents\GreenLab\data\angela\')
load('NAccDecreaseAll.mat')
realAcc = (1-allData.allLambda{1,1}.allErr)*100;
load('NAccDecreaseAllRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end
doubleHist(realAcc,randAcc,'Unit','%','Main','NAcc All: Decrease','Loc','northeast','xlab','Accuracy (%)');

load('NAccIncreaseAll.mat')
realAcc = (1-allData.allLambda{1,1}.allErr).*100;
load('NAccIncreaseAllRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end
doubleHist(realAcc,randAcc,'Unit','%','Main','NAcc All: Increase','Loc','northeast','xlab','Accuracy (%)');
%% NAcc All But One - Increase and Decrease
cd('C:\Users\Pythia\Documents\GreenLab\data\angela\')
load('NAccDecreaseAllBut.mat')
realAcc = (1-allData.allLambda{1,1}.allErr)*100;
load('NAccDecreaseAllButRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end
doubleHist(realAcc,randAcc,'Unit','%','Main','NAcc All but One: Decrease','Loc','northeast','xlab','Accuracy (%)');

load('NAccIncreaseAllBut.mat')
realAcc = (1-allData.allLambda{1,1}.allErr).*100;
load('NAccIncreaseAllButRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end
doubleHist(realAcc,randAcc,'Unit','%','Main','NAcc All but One: Increase','Loc','northeast','xlab','Accuracy (%)');
%% PFC Rest - Increase and Decrease
cd('C:\Users\Pythia\Documents\GreenLab\data\angela\')
load('PFCDecrease.mat')
realAcc = (1-allData.allLambda{1,1}.allErr).*100;
load('PFCDecreaseRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end
doubleHist(realAcc,randAcc,'Unit','%','Main','PFC Rest: Decrease','Loc','northeast','xlab','Accuracy (%)');

load('PFCIncrease.mat')
realAcc = (1-allData.allLambda{1,1}.allErr).*100;
load('PFCIncreaseRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end
doubleHist(realAcc,randAcc,'Unit','%','Main','PFC Rest: Increase','Loc','northeast','xlab','Accuracy (%)');
%% PFC All - Increase and Decrease
cd('C:\Users\Pythia\Documents\GreenLab\data\angela\')
load('PFCDecreaseAll.mat')
realAcc = (1-allData.allLambda{1,1}.allErr)*100;
load('PFCDecreaseAllRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end
doubleHist(realAcc,randAcc,'Unit','%','Main','PFC All: Decrease','Loc','northeast','xlab','Accuracy (%)');

load('PFCIncreaseAll.mat')
realAcc = (1-allData.allLambda{1,1}.allErr).*100;
load('PFCIncreaseAllRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr).*100;
end

doubleHist(realAcc,randAcc,'Unit','%','Main','PFC All: Increase','Loc','northeast','xlab','Accuracy (%)');
%% NAcc Drinking
cd('C:\Users\Pythia\Documents\GreenLab\data\angela\')
load('NAccDrink.mat')
realAcc = allData.allLambda{1,1}.allErr;
load('NAccDrinkRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = allData.allLambda{1,ii}.allErr;
end
doubleHist(realAcc,randAcc,'Unit',' g/kg','Main','NAcc Rest: Baseline Drinking','Loc','northeast','xlab','Mean Absolute Error (g/kg)');

load('NAccDrinkAll.mat')
realAcc = allData.allLambda{1,1}.allErr;
load('NAccDrinkAllRand.mat')
randAcc = zeros(1000,1);
for ii = 1:10
    randAcc(100*ii-99:100*ii,1) = allData.allLambda{1,ii}.allErr;
end
doubleHist(realAcc,randAcc,'Unit',' g/kg','Main','All: Baseline Drinking','Loc','northeast','xlab','Mean Absolute Error (g/kg)');
%% Univariate NAcc
load('angelaAllData.mat')
% Find sets of 5 that have both 1 and 0 out of all possible quintet.
cmbs = nchoosek(1:22,5);
for ii = 1:size(cmbs,1)
    ui(1,ii) = size(unique(nacc(cmbs(ii,:),1)),1);
    ui(2,ii) = size(unique(nacc(cmbs(ii,:),2)),1);
end
cmbs = cmbs(sum(ui,1)==4,:);
% Randomally grab 100 of these quintets
quin = cmbs(randperm(size(cmbs,1),100),:);
inds = 1:22;
%%
for jj = 1:2
    for k = 1:100
        notQuin = inds(~ismember(1:22,quin(k,:)));
        trainX = avgData(notQuin,:);
        trainY = nacc(notQuin,jj);
        testX = avgData(quin(k,:),:);
        testY = nacc(quin(k,:),jj);
        for feat = 1:60
            mdl = fitglm(trainX(:,feat),trainY,'distribution','binomial');
            prob = predict(mdl,testX(:,feat));
            [~,~,~,naccA(jj,k,feat)] = perfcurve(testY,prob,1);
        end
    end
end
%% Univariate PFC
load('angelaAllData.mat')
% Pull out non-existant animals 
avgData = avgData([1:16,19:end],:);
% Find sets of 5 that have both 1 and 0 out of all possible quintet.
cmbs = nchoosek(1:size(pfc,1),5);
for ii = 1:size(cmbs,1)
    ui(1,ii) = size(unique(pfc(cmbs(ii,:),1)),1);
    ui(2,ii) = size(unique(pfc(cmbs(ii,:),2)),1);
end
cmbs = cmbs(sum(ui,1)==4,:);
% Randomally grab 100 of these quintets
quin = cmbs(randperm(size(cmbs,1),100),:);
inds = 1:size(pfc,1);
%%
for jj = 1:2
    for k = 1:100
        notQuin = inds(~ismember(1:size(pfc,1),quin(k,:)));
        trainX = avgData(notQuin,:);
        trainY = pfc(notQuin,jj);
        testX = avgData(quin(k,:),:);
        testY = pfc(quin(k,:),jj);
        for feat = 1:60
            mdl = fitglm(trainX(:,feat),trainY,'distribution','binomial');
            prob = predict(mdl,testX(:,feat));
            [~,~,~,pfcA(jj,k,feat)] = perfcurve(testY,prob,1);
        end
    end
end
%%
load('uniA.mat')
mNaccA = squeeze(mean(naccA,2))';
mPfcA = squeeze(mean(pfcA,2))';
%% Separate NAc responders from mPFC responders (decreased drinking)
data = collateData('C:\Users\Pythia\Documents\GreenLab\data\angela\3second\all\',{'N20';'N24';'N25';'N29';'AH12';'AH16'},{'pow','coh'},'avg','rel');
nacData = [];
for ii = 1:size(data,2)
    for jj = 1:size(data{ii},1)
        nacData = [nacData;data{ii}{jj}];
    end
end
data = collateData('C:\Users\Pythia\Documents\GreenLab\data\angela\3second\all\',{'N33';'AH11';'AH20'},{'pow','coh'},'avg','rel');
pfcData = [];
for ii = 1:size(data,2)
    for jj = 1:size(data{ii},1)
        pfcData = [pfcData;data{ii}{jj}];
    end
end
xData = [nacData;pfcData];
yData = [ones(size(nacData,1),1);zeros(size(pfcData,1),1)];
% save('nacVpfc.mat','xData','yData')
%%
load('C:\Users\Pythia\Documents\GreenLab\data\angela\nacVpfcModel.mat')
rndAcc = [];
for ii = 1:10
    rndAcc = [rndAcc;1-rndData.allLambda{ii}.allErr];
end
doubleHist(1-allData.allLambda{1}.allErr,rndAcc,'xlab','Accuracy','main','NAc vs. mPFC: Optimal Site for Decreased Drinking','loc','ne')
%% Log models - nacVpfc
load('C:\Users\Pythia\Documents\GreenLab\data\angela\nacVpfc.mat')
for ii = 1:60
   mdl = fitglm(xData(:,ii),yData,'distribution','binomial');
   p(ii) = mdl.Coefficients{2,4};
   r(ii) = mdl.Rsquared.Ordinary;
   d(ii) = exp(mdl.Coefficients{2,1});
end
[rSort,sortInd] = sort(r,'descend');
nameVect = names({'SL','SR','PL','PR'},{'d','t','a','b','lg','hg'});
topFeat = nameVect(sortInd)';
sortDir = d(sortInd)';
%% Log models
load('C:\Users\Pythia\Documents\GreenLab\data\angela\angelaAllData.mat')
for ii = 1:60
    % NAc
    mdl = fitglm(avgData(:,ii),nacc(:,1),'distribution','binomial');
    nacR(ii) = mdl.Rsquared.Ordinary;
    nacD(ii) = exp(mdl.Coefficients{2,1});
    % mPFC
    mdl = fitglm(avgData([1:16,19:end],ii),pfc(:,1),'distribution','binomial');
    pfcR(ii) = mdl.Rsquared.Ordinary;
    pfcD(ii) = exp(mdl.Coefficients{2,1});
end
nameVect = names({'SL','SR','PL','PR'},{'d','t','a','b','lg','hg'});
[sNacR,nacRInd] = sort(nacR,'descend');
[sPfcR,pfcRInd] = sort(pfcR,'descend');
sNacR = sNacR';
sPfcR = sPfcR';
nacVar = nameVect(nacRInd)';
pfcVar = nameVect(pfcRInd)';
nacDirSort = nacD(nacRInd)'>1;
pfcDirSort = pfcD(pfcRInd)'>1;
%% Predict response from baseline drinking
% Import the data
[~, ~, raw0_0] = xlsread('C:\Users\Pythia\Documents\GreenLab\data\angela\BAve and Stim Response_For ephys analysis.xlsx','Sheet1','B2:B14');
[~, ~, raw0_1] = xlsread('C:\Users\Pythia\Documents\GreenLab\data\angela\BAve and Stim Response_For ephys analysis.xlsx','Sheet1','D2:D14');
[~, ~, raw0_2] = xlsread('C:\Users\Pythia\Documents\GreenLab\data\angela\BAve and Stim Response_For ephys analysis.xlsx','Sheet1','F2:F14');
raw = [raw0_0,raw0_1,raw0_2];
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));

mdl = fitglm(data(:,1),data(:,2),'distribution','binomial');
mdl2 = fitglm(data(:,1),data(:,3),'distribution','binomial');

figure
plot(data(:,1),data(:,2),'o')
box off
set(gca,'YTick',0:1,'YTickLabel',{'Other','Decrease'})
ylim([-0.1 1.1])
title('NAc: Decrease vs. Other')
xlabel('Baseline Drinking (g/kg)')

figure
plot(data(:,1),data(:,3),'o')
box off
set(gca,'YTick',0:1,'YTickLabel',{'Other','Decrease'})
ylim([-0.1 1.1])
title('mPFC: Decrease vs. Other')
xlabel('Baseline Drinking (g/kg)')
%% PFC Drinking
% cd('C:\Users\Pythia\Documents\GreenLab\data\angela\')
% load('PFCDrinkRest.mat')
% realAcc = allData.allLambda{1,1}.allErr;
% load('PFCDrinkRestRand.mat')
% randAcc = zeros(1000,1);
% for ii = 1:10
%     randAcc(100*ii-99:100*ii,1) = allData.allLambda{1,ii}.allErr;
% end
% doubleHist(realAcc,randAcc,'Unit',' g/kg','Main','PFC Rest: Baseline Drinking','Loc','northeast','xlab','Mean Absolute Error (g/kg)');
% 
% load('PFCDrinkAll.mat')
% realAcc = allData.allLambda{1,1}.allErr;
% load('PFCDrinkAllRand.mat')
% randAcc = zeros(1000,1);
% for ii = 1:10
%     randAcc(100*ii-99:100*ii,1) = allData.allLambda{1,ii}.allErr;
% end
% doubleHist(realAcc,randAcc,'Unit',' g/kg','Main','PFC All: Baseline Drinking','Loc','northeast','xlab','Mean Absolute Error (g/kg)');
%%
% NAcc (probably) not better than chance at predicting response to stim
load('C:\Users\Lucas\Desktop\GreenLab\data\angela\angelaPFCNAcc.mat')
cfg = lassoNetCfg([],'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(NAccX,NAccY,'binomial','class',1,4,1,cfg);
masterRealAcc = 1-allLambda{1,1}.allErr;
% Randomize
cfg = lassoNetCfg([],'y','y','n',100,'1se');
[randAlpha,randLambda,randBeta,randCvFitsArray,randAccArray,randHist] = lassoNet(NAccX,NAccY,'binomial','class',1,4,100,cfg);
masterRandErr = [];
for ii = 1:size(randLambda,2)
    masterRandErr = [masterRandErr;randLambda{1,ii}.allErr];
end
masterRandAcc = 1-masterRandErr;
%%
save('C:\Users\Lucas\Desktop\GreenLab\data\angela\naccData.mat','accArray','allAlpha','allBeta','allLambda','hist','cvFitsArray','masterRandAcc','masterRealAcc','randAccArray','randAlpha','randCvFitsArray','randHist','randLambda')
%% PFC
load('C:\Users\Lucas\Desktop\GreenLab\data\angela\angelaPFCNAcc.mat')
cfg = lassoNetCfg([],'n','y','n',100,'1se');
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(PFCX,PFCY,'binomial','class',1,4,1,cfg);
masterRealAcc = 1-allLambda{1,1}.allErr;
% Randomize
cfg = lassoNetCfg([],'y','y','n',100,'1se');
[randAlpha,randLambda,randBeta,randCvFitsArray,randAccArray,randHist] = lassoNet(PFCX,PFCY,'binomial','class',1,4,100,cfg);
masterRandErr = [];
for ii = 1:size(randLambda,2)
    masterRandErr = [masterRandErr;randLambda{1,ii}.allErr];
end
masterRandAcc = 1-masterRandErr;
%%
save('C:\Users\Lucas\Desktop\GreenLab\data\angela\pfcData.mat','accArray','allAlpha','allBeta','allLambda','hist','cvFitsArray','masterRandAcc','masterRealAcc','randAccArray','randAlpha','randCvFitsArray','randHist','randLambda')
%%
load('C:\Users\Lucas\Desktop\GreenLab\data\angela\pfcData.mat')
realMean = mean(masterRealAcc);
randMean = mean(masterRandAcc);
figure
ecdf(masterRealAcc)
hold on
ecdf(masterRandAcc)
xlabel('Accuracy')
ylabel('Cumulative Density')
title('Stimulation Response Prediction: PFC')
legend({['Real: ',num2str(round((realMean*100),2)),'%'],['Rand: ',num2str(round((randMean*100),2)),'%']},'location','northwest')
%%
bInd = logicFind(0.2,allBeta{1,1}.survBeta,'>');
load('varNames.mat')
vars = names(bInd);
%%
load('C:\Users\Lucas\Desktop\GreenLab\data\angela\naccData.mat')
realMean = mean(masterRealAcc);
randMean = mean(masterRandAcc);
figure
ecdf(masterRealAcc)
hold on
ecdf(masterRandAcc)
xlabel('Accuracy')
ylabel('Cumulative Density')
title('Stimulation Response Prediction: NAcc')
legend({['Real: ',num2str(round((realMean*100),2)),'%'],['Rand: ',num2str(round((randMean*100),2)),'%']},'location','northwest')
bInd = logicFind(0.2,allBeta{1,1}.survBeta,'>');
load('varNames.mat')
vars = names(bInd);