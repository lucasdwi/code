%% MIA baseline3 kamiak
baseFolder='~/LDcode';
addpath(genpath(baseFolder))
% duals vs. all with 5 second bins as samples
[data,~,files] = collateData('C:\Users\angela.henricks\OneDrive - Washington State University (email.wsu.edu)\Desktop\Research\MatlabAnalysis\MIA\BaselineMIA\processed\',{'Con';'MIA';'AE';'Dual'},{'pow','coh'},'trl','rel');
Con=cat(1,data{1}{:});
MIA=cat(1,data{2}{:});
AE=cat(1,data{3}{:});
Dual=cat(1,data{4}{:});
%% Dual vs ALL by trials
for ii = 1:100
    disp(ii)
    % Generate 15 random indices from dual and 3 from AE, MIA, and Con each 
    aeTrainInd = randperm(size(AE,1),1500);
    conTrainInd = randperm(size(Con,1),1500);
    miaTrainInd = randperm(size(MIA,1),1500);
    dualTrainInd = randperm(size(Dual,1),4500);
    % Use inds to pull out training set
    thisTrainAE = AE(aeTrainInd,:);
    thisTrainCon = Con(conTrainInd,:);
    thisTrainMIA = MIA(miaTrainInd,:);
    thisTrainDual = Dual(dualTrainInd,:);
    % Group together and generate corresponding classifier groups; 1 =
    % Dual, 0 = all other groups
    trainX = [thisTrainDual;thisTrainAE;thisTrainCon;thisTrainMIA];
    trainY = [ones(4500,1);zeros(4500,1)];
    % Use the remaining 2 animals from Dual and one animal from each other
    % group for the test set
    dualTestInd = 1:5627;
    dualTestInd = dualTestInd(~ismember(dualTestInd,dualTrainInd));
    dualTest = Dual(dualTestInd,:);
    
    aeTestInd = 1:4877;
    aeTestInd = aeTestInd(~ismember(aeTestInd,aeTrainInd));
    aeTest = AE(aeTestInd(1),:);
    conTestInd = 1:7680;
    conTestInd = conTestInd(~ismember(conTestInd,conTrainInd));
    conTest = Con(conTestInd(1),:);
    miaTestInd = 1:9872;
    miaTestInd = miaTestInd(~ismember(miaTestInd,miaTrainInd));
    miaTest = MIA(miaTestInd(1),:);

    % Combine and generate classifier group vector
    testX = [aeTest;conTest;miaTest;dualTest];
    testY = [zeros(size([aeTest;conTest;miaTest],1),1);...
        ones(size(dualTest,1),1)];
    % Lasso real
    cfg = lassoNetCfg({testX,testY},[],'n','y','n',100,'1se',[]);
    [~,lambda,beta,fits,acc,hist] = lassoNet(trainX,trainY,'binomial',...
        'class',1,5,1,cfg);
    err(:,ii) = lambda{1}.allErr;
    % accuracy
    a(ii) = acc{1}.acc;
    % Lasso random
    cfg = lassoNetCfg({testX,testY},[],'y','y','n',100,'1se',[]);
    [~,lambdaR,betaR,fitsR,accR,histR] = lassoNet(trainX,trainY,'binomial',...
        'class',1,5,1,cfg);
    errR(:,ii) = lambdaR{1}.allErr;
    aR(ii) = accR{1}.acc;
    % Single feature
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial','binomialsize',30);
        pred = predict(mdl,testX(:,jj));
        [~,~,~,auc(ii,jj)] = perfcurve(testY,pred,1);
    end
end
% plot
% doubleHist(allA,allAR);
pd = fitdist(allA','normal');
pdR = fitdist(allAR','normal');
x = 0:0.01:1;
y = pdf(pd,x);
yR = pdf(pdR,x);
figure
hold on
plot(x,y)
plot(x,yR)
legend({['Real: \mu = ',num2str(round(mean(allA),2)),'\pm',...
    num2str(round(conf(allA,0.95),2))],...
    ['Permuted: \mu = ',num2str(round(mean(allAR),2)),'\pm',...
    num2str(round(conf(allAR,0.95),2))]})
%% names of features
nameVect = names({'ILL','CAL','PLL','NAcL','PLR','CAR','ILR','NAcR'},{'d','t','a','b','lg','hg'});
[allAUCs,sortInd] = sort(mean(allAUC,1),'descend');
sortName = nameVect(sortInd);
DAtiersAUC = tier(allAUC(:,sortInd));