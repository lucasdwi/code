%% In session (notDrink) predicting drinking
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\paper3\drinkNot\'],{'.mat'},{'pow','coh'},'avg','rel');
mInd = [1:18,24,25];
male = cat(1,data{1}{mInd,2});
maleDrink = [0.03;0.00;0.63;0.16;0.14;0.21;0.29;0.16;0.25;0.41;0.05;...
    0.32;0.09;0.04;0.04;0.02;0.11;0.02;0.06;0.02];

% Subset female to have same number of sessions per animal as male
subInd = [20:22,26:40,41:42];
% fInd = [19:23,26:45];
female = cat(1,data{1}{subInd,2});
femaleDrink = [0.39;0.39;0.71;0.08;0.24;0.47;0.89;1.39;0.52;0.90;0.57;...
    0.86;1.13;1.20;0.11;0.73;0.99;0.73;0.83;1.08;0.18;0.40;0.28;0.14;0.11];
femaleDrink = femaleDrink([1:3,6:22]);
save(['C:\Users\Pythia\Documents\GreenLab\data\maleFemale\'...
    'drinkAmountData.mat'],'male','maleDrink','female','femaleDrink')
%%
load('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\drinkAmount.mat')
load(['C:\Users\Pythia\Documents\GreenLab\data\maleFemale\'...
    'drinkAmountData.mat'],'maleDrink','femaleDrink')
doubleHist(maleDat.allLambda.allErr,rndMale.allLambda.allErr,...
    'xlab','Deviance (g/kg)','main','male')
hold on
plot([mean(maleDrink) mean(maleDrink)],[0 0.45],'--k')
plot([median(maleDrink) median(maleDrink)],[0 0.45],'--','color',...
    [0.5 0.5 0.5])

doubleHist(femaleDat.allLambda.allErr,rndFemale.allLambda.allErr,...
    'xlab','Deviance (g/kg)','main','female')
hold on
plot([mean(femaleDrink) mean(femaleDrink)],[0 0.3],'--k')
plot([median(femaleDrink) median(femaleDrink)],[0 0.3],'--','color',...
    [0.5 0.5 0.5])
%% Get features; rank with deviance from univariates using all samples
load(['C:\Users\Pythia\Documents\GreenLab\data\maleFemale\'...
    'drinkAmountData.mat'])
for ii = 1:60
    mMdl = fitglm(male(:,ii),maleDrink);
    mDev(ii) = mMdl.Deviance;
    mR(ii) = mMdl.Rsquared.Ordinary;
    mCoeff(ii) = table2array(mMdl.Coefficients(2,1));
    
    fMdl = fitglm(female(:,ii),femaleDrink);
    fDev(ii) = fMdl.Deviance;
    fR(ii) = fMdl.Rsquared.Ordinary;
    fCoeff(ii) = table2array(fMdl.Coefficients(2,1));
end
% Rank features by r-squared
[mRSort,mRInd] = sort(mR,'descend');
[fRSort,fRInd] = sort(fR,'descend');
% Transpose
mRSort = mRSort';
fRSort = fRSort';
% Sort coeffs
mCoeffSort = mCoeff(mRInd)';
fCoeffSort = fCoeff(fRInd)';
% Get feature names
nameVect = names({'PL','PR','SL','SR'},{'d','t','a','b','lg','hg'});
mFeat = nameVect(mRInd)';
fFeat = nameVect(fRInd)';
%% Baseline 'trait' predicting drinking
[data,~,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\data\'...
    'maleFemale\3SecondRest\'],{'.mat'},{'pow','coh'},'avg','rel');
