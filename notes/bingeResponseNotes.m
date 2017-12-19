load('bingeRespRestOutcome.mat')
shellAcc = (1-allData.allLambda{1,1}.allErr(:,1)).*100;
coreAcc = (1-allData.allLambda{1,1}.allErr(:,2)).*100;

load('bingeRespRestRandOutcome.mat')
randShellAcc = zeros(1000,1);
for ii = 1:10
    randShellAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr(:,1)).*100;
end

randCoreAcc = zeros(1000,1);
for ii = 1:10
    randCoreAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr(:,2)).*100;
end
%%
doubleHist(shellAcc,randShellAcc,'Unit','%','xlab','Accuracy (%)','main','Shell Response: Abs Rest')
doubleHist(coreAcc,randCoreAcc,'Unit','%','xlab','Accuracy (%)','main','Core Response: Abs Rest')
%%
load('bingeRespAllOutcome.mat')
shellAcc = (1-allData.allLambda{1,1}.allErr(:,1)).*100;
coreAcc = (1-allData.allLambda{1,1}.allErr(:,2)).*100;

load('bingeRespAllRandOutcome.mat')
randShellAcc = zeros(1000,1);
for ii = 1:10
    randShellAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr(:,1)).*100;
end

randCoreAcc = zeros(1000,1);
for ii = 1:10
    randCoreAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr(:,2)).*100;
end

doubleHist(shellAcc,randShellAcc,'Unit','%','xlab','Accuracy (%)','main','Shell Response: Abs All')
doubleHist(coreAcc,randCoreAcc,'Unit','%','xlab','Accuracy (%)','main','Core Response: Abs All')
%%
load('bingeRespAllOutcomeRel.mat')
shellAcc = (1-allData.allLambda{1,1}.allErr(:,1)).*100;
coreAcc = (1-allData.allLambda{1,1}.allErr(:,2)).*100;

load('bingeRespAllOutcomeRelRand.mat')
randShellAcc = zeros(1000,1);
for ii = 1:10
    randShellAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr(:,1)).*100;
end

randCoreAcc = zeros(1000,1);
for ii = 1:10
    randCoreAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr(:,2)).*100;
end

doubleHist(shellAcc,randShellAcc,'Unit','%','xlab','Accuracy (%)','main','Shell Response: Rel All')
doubleHist(coreAcc,randCoreAcc,'Unit','%','xlab','Accuracy (%)','main','Core Response: Rel All')
%%
load('bingeRespRestOutcomeRel.mat')
shellAcc = (1-allData.allLambda{1,1}.allErr(:,1)).*100;
coreAcc = (1-allData.allLambda{1,1}.allErr(:,2)).*100;

load('bingeRespRestOutcomeRelRand.mat')
randShellAcc = zeros(1000,1);
for ii = 1:10
    randShellAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr(:,1)).*100;
end

randCoreAcc = zeros(1000,1);
for ii = 1:10
    randCoreAcc(100*ii-99:100*ii,1) = (1-allData.allLambda{1,ii}.allErr(:,2)).*100;
end

doubleHist(shellAcc,randShellAcc,'Unit','%','xlab','Accuracy (%)','main','Shell Response: Rel Rest')
doubleHist(coreAcc,randCoreAcc,'Unit','%','xlab','Accuracy (%)','main','Core Response: Rel Rest')