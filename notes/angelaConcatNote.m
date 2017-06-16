% Collate data without N20
fNames = fileSearch('C:\Users\Lucas\Desktop\GreenLab\data\angela\processed','.mat','in','N20','ex');
[data,samp] = collateData('C:\Users\Lucas\Desktop\GreenLab\data\angela\processed\',[],{'pow','coh','corr'},'avg',fNames);
%%
% Get names of all files for day1 and day2; exclude N20 (file needing imputing)
fNames1 = fileSearch('C:\Users\Lucas\Desktop\GreenLab\data\angela\processed','base1','in','N20','ex');
fNames2 = fileSearch('C:\Users\Lucas\Desktop\GreenLab\data\angela\processed','base2','in','N20','ex');
%% Load data into corresponding cell arrays
for fI = 1:size(fNames1,2)
    % Load day 1
    load(fNames1{fI})
    coh1{fI} = mean(coh{1,1}.rel,3);
    pow1{fI} = psdTrls{1,1}.avgRelPow;
    corr1{fI} = rVect{1,1};
    % Load day 2
    load(fNames2{fI})
    coh2{fI} = mean(coh{1,1}.rel,3);
    pow2{fI} = psdTrls{1,1}.avgRelPow;
    corr2{fI} = rVect{1,1};
    % Get distributions of change
    powDelta(:,:,fI) = (pow2{fI} - pow1{fI})./pow1{fI};
    cohDelta(:,:,fI) = (coh2{fI} - coh1{fI})./coh1{fI};
    corrDelta(:,:,fI) = (corr2{fI} - corr1{fI})./corr1{fI};
end
%% Get average of delta
powDeltaMu = mean(powDelta,3);
cohDeltaMu = mean(cohDelta,3);
corrDeltaMu = mean(corrDelta,3);
% Get standard devation of delta
powDeltaSigma = std(powDelta,[],3);
cohDeltaSigma = std(cohDelta,[],3);
corrDeltaSigma = std(corrDelta,[],3);
%%
load('N20_Base1_10_26_16_rest.mat')
pow1 = psdTrls{1,1}.avgRelPow;
coh1 = mean(coh{1,1}.rel,3);
corr1 = rVect{1,1};

powImp = (1+powDeltaMu).*pow1;
cohImp = (1+cohDeltaMu).*coh1;
corrImp = (1+corrDeltaMu).*corr1;

load('N20_Base2_11_4_16_rest.mat')
pow2 = psdTrls{1,1}.avgRelPow;
coh2 = mean(coh{1,1}.rel,3);
corr2 = rVect{1,1};
%%
% Inds to impute: 
% Pow: channels 1 and 2
pow2(:,1:2) = powImp(:,1:2);
% Coherence: all combinations including 1 and 2
coh2(1:5,:) = cohImp(1:5,:);
% Correlations: channels 1 and 2
corr2(1:5,:) = corrImp(1:5,:);
% Reshape and concatenate data
thisData = [reshape(pow1,1,20),reshape(coh1,1,30),reshape(corr1,1,30);reshape(pow2,1,20),reshape(coh2,1,30),reshape(corr2,1,30)];
%%
mostDataNacc = cat(1,data{1,1}{:});
mostYnacc = [0,0,0,0,1,1,1,1,1,1,1,1,0,0,1,1,1,1,0,0,0,0];
mostDataPfc = cat(1,data{1,1}{[1:6,9:16,19:22]});
mostYpfc = [1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0];

allDataNacc = [mostDataNacc;thisData];
allDataPfc = [mostDataPfc;thisData];
allYnacc = [mostYnacc,0,0];
allYpfc = [mostYpfc,1,1];
%%
load('Nacc.mat')
realAcc = 1-allData.allLambda{1,1}.allErr;
load('NaccRand.mat')
randAcc = [];
for ii = 1:100
    randAcc = [randAcc;1-allData.allLambda{1,ii}.allErr];
end 
figure
histogram(realAcc,'Normalization','probability','BinWidth',0.02)
hold on
histogram(randAcc,'Normalization','probability','BinWidth',0.02)
legend({'Real','Permuted'},'Location','northeast')
d = distES(realAcc,randAcc);
text(0.7,0.14,['d = ',num2str(d)])
title('Real vs. Permuted Data: NAcc Response')
xlabel('Accuracy')
ylabel('Proportion of CV Models')
%%
load('PFC.mat')
realAcc = 1-allData.allLambda{1,1}.allErr;
load('PFCrand.mat')
randAcc = [];
for ii = 1:100
    randAcc = [randAcc;1-allData.allLambda{1,ii}.allErr];
end 
figure
histogram(realAcc,'Normalization','probability','BinWidth',0.02)
hold on
histogram(randAcc,'Normalization','probability','BinWidth',0.02)
legend({'Real','Permuted'},'Location','northeast')
d = distES(realAcc,randAcc);
text(0.7,0.09,['d = ',num2str(d)])
title('Real vs. Permuted Data: PFC Response')
xlabel('Accuracy')
ylabel('Proportion of CV Models')