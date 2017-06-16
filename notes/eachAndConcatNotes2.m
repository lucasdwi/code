%% Plot features with significant correlations with changes in voracity
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\featCorr.mat')
nameVect = names({'SL','CL','SR','CR'},{'d','t','a','b','lg','hg'});
% Add worst p-value to list
for ii = 1:length(pInds)
    figure
    % Convert features to gm/ms and %
    scatter(1000.*vorDiff(~isnan(vorDiff)),100.*feats(:,pInds(ii)),'ok','Filled')
    lsline
    thiscorr{ii} = fitlm(1000.*vorDiff(~isnan(vorDiff)),feats(:,pInds(ii)));
    title(nameVect{pInds(ii)})
    text(0,0,['R^2 = ',num2str(round(thiscorr{ii}.Rsquared.Ordinary,2)),char(10),'P = ',num2str(round(thiscorr{ii}.Coefficients.pValue(2),3))])
    xlabel('Change in Voracity (gm/ms)')
    ylabel('Percent Change in Feature')
end
%%
% for ii = 1:length(pInds)
%    tbl(ii,1) = thiscorr{ii}.Rsquared.Ordinary;
%    tbl(ii,2) = thiscorr{ii}.Coefficients.pValue(2);
%    vars{ii} = nameVect{pInds(ii)};
% end
% tbl = array2table(tbl);
% tbl.Properties.VariableNames = {'R2','pValue'};
% tbl.Properties.RowNames{1} = 'SLa';
% tbl.Properties.RowNames{2} = 'SLSRlg';
% tbl.Properties.RowNames{3} = 'SRCRt';
% tbl.Properties.RowNames{4} = 'SLtCL';
% tbl.Properties.RowNames{5} = 'CLtCR';
%% Load all concat files and process
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\500Train_50-50')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\Base1')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\57\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\preBinge\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    ccA(ii) = concatData.auc;
    ccX{ii} = concatData.rocX;
    ccY{ii} = concatData.rocY;
    %
    ceA(ii,:) = eachData.auc;
    ceX(ii,:) = eachData.rocX;
    ceY(ii,:) = eachData.rocY;
end
% Get average and std of 'eachA'
ceAM = mean(ceA,1);
ceS = std(ceA,[],1);
% Get average ROC for 'each'
for ii = 1:12
    ceXM(ii,:) = mean(cat(2,ceX{:,ii}),2);
    ceYM(ii,:) = mean(cat(2,ceY{:,ii}),2);
    ceXS(ii,:) = std(cat(2,ceX{:,ii}),[],2);
    ceYS(ii,:) = std(cat(2,ceY{:,ii}),[],2);
end
% Get average ROC for concat-concat models
ccMX = mean(cat(2,ccX{:}),2);
ccMY = mean(cat(2,ccY{:}),2);
% Get average 'fill' for concat-concat models
[ccXfill,ccYfill] = avgFill(cat(2,ccX{:}),cat(2,ccY{:}),2,1);
%%
% figure
% fill(ccXfill,ccYfill,[.8 .8 .8])
% hold on
% plot(ccMX,ccMY,'-k')
% ylim([0 1])
% xlim([0 1])
%% Load all concatRand files
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\500Train_50-50Rand')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\base1Rand\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\57Rand\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\preBingeRand\')
for ii = 1:20
   load([num2str(ii),'.mat'])
   % Check if roc is 50 line, if so interpolate it; store curve
   if isequal(concatData.rocX,[0;1])
      ccRandX{ii} = (0:1/1500:1)'; 
      ccRandY{ii} = (0:1/1500:1)';
   else
       ccRandX{ii} = concatData.rocX;
       ccRandY{ii} = concatData.rocY;
   end
   % Check if roc is 50 line, if so interpolate it; store curve
   for jj = 1:12
       if isequal(eachData.rocX{jj},[0;1])
           ceRandX{ii}{jj} = (0:1/125:1)';
           ceRandY{ii}{jj} = (0:1/125:1)';
       else
           ceRandX{ii}{jj} = eachData.rocX;
           ceRandY{ii}{jj} = eachData.rocY;
       end
   end
   % Store auc value
    ccRandA(ii) = concatData.auc;
    ceRandA(ii,:) = eachData.auc;
end
% Get average ROC for concat-concat random models
ccRandMX = mean(cat(2,ccRandX{:}),2);
ccRandMY = mean(cat(2,ccRandY{:}),2);
% Get average 'fill' for concat-concat random models
[ccRandXfill,ccRandYfill] = avgFill(cat(2,ccRandX{:}),cat(2,ccRandY{:}),2,1);
%% Load all 'each' files
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\500Trials_50-50\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\base1\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\57\')
for ii = 1:240
    load([num2str(ii),'fix.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    % Get each and self auc and rocs
    eachData.auc(animal) = selfData.auc;
    eeA(animal,:,iter) = eachData.auc;
    eachData.rocX{animal} = selfData.rocX;
    eeX{animal}(:,:,iter) = cat(2,eachData.rocX{:});
    eachData.rocY{animal} = selfData.rocY;
    eeY{animal}(:,:,iter) = cat(2,eachData.rocY{:});
    % Get concat auc and rocs
    ecA(animal,iter) = concatData.auc;
    ecX(animal,:,iter) = concatData.rocX;
    ecY(animal,:,iter) = concatData.rocY;
end
%% Load all 'each' random files 
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\500Trials_50-50Rand\')
% cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\57Rand\')
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    % Add self data to each
    eachData.rocX{animal} = selfData.rocX;
    eachData.rocY{animal} = selfData.rocY;
    % Check if roc is 50 line, if so interpolate it; store curve
    if isequal(concatData.rocX,[0;1])
        ecRandX{animal,iter} = (0:1/1500:1)';
        ecRandY{animal,iter} = (0:1/1500:1)';
    else
        ecRandX{animal,iter} = concatData.rocX;
        ecRandY{animal,iter} = concatData.rocY;
    end
    % Check if roc is 50 line, if so interpolate it; store curve
    if isequal(eachData.rocX{animal},[0,1])
        eeRandX{animal,iter} = (0:1/125:1)';
        eeRandY{animal,iter} = (0:1/125:1)';
    else
        eeRandX{animal,iter} = eachData.rocX;
        eeRandY{animal,iter} = eachData.rocY;
    end
    % Get each and self AUC
    eachData.auc(animal) = selfData.auc;
    eeRandA(animal,:,iter) = eachData.auc;
    % Get concat AUC
    ecRandA(ii) = concatData.auc;
end
%% Plot 'selfs' - diag
cols = distinguishable_colors(12);
figure
hold on
for ii = 1:12
   selfXm(ii,:) = mean(eeX{1,ii}(:,ii,:),3);
   selfYm(ii,:) = mean(eeY{1,ii}(:,ii,:),3);
   plot(mean(eeX{1,ii}(:,ii,:),3),mean(eeY{1,ii}(:,ii,:),3),'Color',cols(ii,:)) 
end
plot(mean(selfXm,1),mean(selfYm,1),'-k','LineWidth',2)
plot([0 1],[0 1],'--k')
title('Each to Self: Diagonal')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
text(0.8,0.2,['AUC: ',num2str(round(mean(diag(mean(eeA,3))),2))])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
%% Plot inds to others - column
cols = distinguishable_colors(12);
mX = squeeze(mean(eeX{1,1},3));
mY = squeeze(mean(eeY{1,1},3));
% Get average curve
otherXM = mean(mX,2);
otherYM = mean(mY,2);
figure
hold on
for ii = 1:12
    plot(mX(:,ii),mY(:,ii),'Color',cols(ii,:))
end
% Plot average
plot(otherXM,otherYM,'-k','LineWidth',2)
plot([0 1],[0 1],'--k')
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
text(0.8,0.2,['AUC: ',num2str(round(mean(mean(eeA(:,1,:),3),1),2))])
title('Each to Each: Column')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
%% Plot all ESs
% Convert AUCs to ESs
for ii = 1:12
    for jj = 1:12
        [eeES(ii,jj),~] = distES(squeeze(eeA(ii,jj,:)),squeeze(eeRandA(ii,jj,:))); 
    end
end
figure
pcolor(padarray(eeES,[1,1],'post'))
set(gca,'XTick',1.5:12.5,'XTickLabel',1:12,'YTick',1.5:12.5,'YTickLabel',1:12)
colormap('viridis')
colorbar
title('Effect Size')
xlabel('Training Set')
ylabel('Test Set')
%% Plot all AUCs
% Average
figure
pcolor(padarray(mean(eeA,3),[1,1],'post'))
set(gca,'XTick',1.5:12.5,'XTickLabel',1:12,'YTick',1.5:12.5,'YTickLabel',1:12)
colormap('viridis')
colorbar
title('AUC: Average')
ylabel('Training Set')
xlabel('Test Set')
% Standard deviation
figure
pcolor(padarray(std(eeA,[],3),[1,1],'post'))
set(gca,'XTick',1.5:12.5,'XTickLabel',1:12,'YTick',1.5:12.5,'YTickLabel',1:12)
colormap('viridis')
colorbar
title('AUC: Standard Deviation')
ylabel('Training Set')
xlabel('Test Set')
%% Plot concat-concat; concat-each
% Get effect size of concat-concat
[d,~] = distES(ccA,ccRandA);
figure
hold on
% Concat-each
for ii = 1:12
   plot(ceXM(ii,:),ceYM(ii,:),'Color',cols(ii,:)) 
end
% Concat-concat fill
fill(ccXfill,ccYfill,[0.8 0.8 0.8])
% Concat-concat average
plot(ccMX,ccMY,'-k','LineWidth',2)
% Concat-concat fill permuted
fill(ccRandXfill,ccRandYfill,[0.8 0.8 0.8])
% Concat-concat average permuted
plot(ccRandMX,ccRandMY,'--k','LineWidth',2)
xlim([0 1])
ylim([0 1])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
title('Concat vs. Concat-Each-Perm')
xlabel('False Positive Rate')
ylabel('True Positive Rate')

text(0.7,0.2,['d = ',num2str(round(d,2))])
text(0.7,0.15,['AUC = ',num2str(round(mean(ccA),2))])
%% Plot example of each-self vs. concat-each
% Get average each-each
for ii = 1:12
   eeXM(:,:,ii) = mean(eeX{ii},3);
   eeYM(:,:,ii) = mean(eeY{ii},3);
end
% Get AUC diffs
diffAUC = diag(mean(eeA,3))-mean(ceA,1)';
% Plot both ROCs and area
figure
hold on
fill([eeXM(:,1,1);flipud(ceXM(1,:)')],[eeYM(:,1,1);flipud(ceYM(1,:)')],[0.8 0.8 0.8]);
plot(eeXM(:,1,1),eeYM(:,1,1),'-k','LineWidth',2)
plot(ceXM(1,:),ceYM(1,:),'--k','LineWidth',2)
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
title('Concat vs. Self: Example')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend({['AUC Diff = ',num2str(round(diffAUC(1),2))],'Self','Concat'})
% Plot distribution of auc differences
figure
scatter(ones(1,12),diffAUC,400,'.k')
set(gca,'XTick',[],'XTickLabel',{})
title('AUC Differences: Self-Concat')
ylabel('AUC Difference')
%% Plot SVM data
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\svm2.mat')
ccSVMa = concatA;
ccSVMx = mean(cat(2,concatX{:}),2);
ccSVMy = mean(cat(2,concatY{:}),2);
[ccSVMxFill,ccSVMyFill] = avgFill(cat(2,concatX{:}),cat(2,concatY{:}),2);
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\svm2Rand.mat')
ccSVMaRand = concatA;
ccSVMxRand = mean(cat(2,concatX{:}),2);
ccSVMyRand = mean(cat(2,concatY{:}),2);
[ccSVMxFillRand,ccSVMyFillRand] = avgFill(cat(2,concatX{:}),cat(2,concatY{:}),2);
% Get effect size
[svmES,~] = distES(ccSVMa,ccSVMaRand);
%%
figure
hold on
% Add concat-concat Lasso
% fill(ccXfill,ccYfill,[0.8 0.8 0.8])
% h1 = plot(ccMX,ccMY,':k','LineWidth',2);
% Then SVM
fill(ccSVMxFill,ccSVMyFill,[0.8 0.8 0.8])
h2 = plot(ccSVMx,ccSVMy,'-k','LineWidth',2);
fill(ccSVMxFillRand,ccSVMyFillRand,[0.8 0.8 0.8])
h3 = plot(ccSVMxRand,ccSVMyRand,'--k','LineWidth',2);
legend([h2,h3],{'SVM','Permuted'},'Location','southeast')
xlim([0 1])
ylim([0 1])
text(0.8,0.3,['d = ',num2str(round(svmES,2))])
text(0.8,0.25,['AUC = ',num2str(round(mean(ccSVMa),2))])
set(gca,'XTick',[0:0.5:1],'YTick',[0:0.5:1])
xlabel('False Postive Rate')
ylabel('True Positive Rate')
title('Concat: SVM')
%% Compare real and permuted auc distributions using Mann-Whitney U
[dC,pC] = distES(ccA,ccRandA);
for ii = 1:12
    [dE(ii),pE(ii)] = distES(ceA(:,ii),ceRandA(:,ii));
end
[~,~,pAdj] = fdr_bh([pC,pE],0.05,'dep');
%% Build and test straight logistic models
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat50-50_500Train.mat')
for ii = 1:20
    mdl{ii} = fitglm(allTrainX{ii},allTrainY{ii},'distribution','binomial');
    prob(:,ii) = predict(mdl{ii},allTestX{ii});
    [concatX(:,ii),concatY(:,ii),~,concatA(ii)] = perfcurve(allTestY{ii},prob(:,ii),1);
end
ccLogA = concatA;
ccLogMX = mean(concatX,2);
ccLogMY = mean(concatY,2);
[ccLogXfill,ccLogYfill] = avgFill(concatX,concatY,2);
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat50-50_500TrainRand.mat')
for ii = 1:20
    mdlRand{ii} = fitglm(allTrainX{ii},allTrainY{ii},'distribution','binomial');
    prob(:,ii) = predict(mdlRand{ii},allTestX{ii});
    [concatRandX(:,ii),concatRandY(:,ii),~,concatRandA(ii)] = perfcurve(allTestY{ii},prob(:,ii),1);
end
ccLogRandA = concatRandA;
ccLogMXrand = mean(concatRandX,2);
ccLogMYrand = mean(concatRandY,2);
[ccLogXrandFill,ccLogYrandFill] = avgFill(concatRandX,concatRandY,2);
[ccES,~] = distES(ccLogA,ccLogRandA);
%%  Plot concat-concat model with permuted data
figure
hold on
% Add concat-concat Lasso
% fill(ccXfill,ccYfill,[0.8 0.8 0.8])
% h1 = plot(ccMX,ccMY,':k','LineWidth',2);
% Then full logisitic
fill(ccLogXfill,ccLogYfill,[0.8 0.8 0.8])
h2 = plot(ccLogMX,ccLogMY,'-k','LineWidth',2);
fill(ccLogXrandFill,ccLogYrandFill,[0.8 0.8 0.8])
h3 = plot(ccLogMXrand,ccLogMYrand,'--k','LineWidth',2);
legend([h2,h3],{'Logistic','Permuted'},'Location','southeast')
xlim([0 1])
ylim([0 1])
xlabel('False Positive Rate')
ylabel('True Positive Rate')
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
title('Concat-Concat: Logistic')
text(0.8,.3,['d = ',num2str(round(ccES,2))])
text(0.8,.25,['AUC = ',num2str(round(mean(concatA),2))])
%% Test models made from concat on other conditions
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\dep24Test.mat')
for ii = 1:20
   load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\500Train_50-50\',num2str(ii),'.mat'])
   testX = [];
   testY = [];
   testX = allTestX{ii};%(:,[1:2,4:34,36:55,57:60]);
   testY = allTestY{ii};
   prob = cvglmnetPredict(concatData.model,zscore(testX),'lambda_1se','response');
   [dep24X{ii},dep24Y{ii},~,dep24A(ii)] = perfcurve(testY,prob,1);
end
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\dep48Test.mat')
for ii = 1:20
   load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\500Train_50-50\',num2str(ii),'.mat'])
   testX = [];
   testY = [];
   testX = allTestX{ii};%(:,[1:2,4:34,36:55,57:60]);
   testY = allTestY{ii};
   prob = cvglmnetPredict(concatData.model,zscore(testX),'lambda_1se','response');
   [dep48X{ii},dep48Y{ii},~,dep48A(ii)] = perfcurve(testY,prob,1);
end
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\chowTest.mat')
for ii = 1:20
   load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\500Train_50-50\',num2str(ii),'.mat'])
   testX = [];
   testY = [];
   testX = allTestX{ii};%(:,[1:2,4:34,36:55,57:60]);
   testY = allTestY{ii};
   prob = cvglmnetPredict(concatData.model,zscore(testX),'lambda_1se','response');
   [chowX{ii},chowY{ii},~,chowA(ii)] = perfcurve(testY,prob,1);
end
%% Prepare for plotting
dep24XM = mean(cat(2,dep24X{:}),2);
dep24YM = mean(cat(2,dep24Y{:}),2);
[dep24Xfill,dep24Yfill] = avgFill(cat(2,dep24X{:}),cat(2,dep24Y{:}),2);
dep48XM = mean(cat(2,dep48X{:}),2);
dep48YM = mean(cat(2,dep48Y{:}),2);
[dep48Xfill,dep48Yfill] = avgFill(cat(2,dep48X{:}),cat(2,dep48Y{:}),2);
chowXM = mean(cat(2,chowX{:}),2);
chowYM = mean(cat(2,chowY{:}),2);
[chowXfill,chowYfill] = avgFill(cat(2,chowX{:}),cat(2,chowY{:}),2);
%% Plot
figure
hold on
% Add concat-concat line
plot(ccMX,ccMY,'-k','LineWidth',2)
% fill(dep24Xfill,dep24Yfill,'b','FaceAlpha',0.5)
plot(dep24XM,dep24YM,'-k','LineWidth',1)
% fill(dep48Xfill,dep48Yfill,'r','FaceAlpha',0.5)
plot(dep48XM,dep48YM,'--k','LineWidth',1.5)
% fill(chowXfill,chowYfill,'y','FaceAlpha',0.5)
plot(chowXM,chowYM,':k','LineWidth',1.5)
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend({['Base: ',num2str(round(mean(ccA),2))],['Dep24: ',num2str(round(mean(dep24A),2))],['Dep48: ',num2str(round(mean(dep48A),2))],['Chow: ',num2str(round(mean(chowA),2))]},'Location','southeast')
title('Baseline Model across Conditions')
ylim([0 1])
xlim([0 1])
%% Test models made from each on other conditions
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\dep24Test.mat')
dep24TestX = eachTestX;
dep24TestY = eachTestY;
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\dep48Test.mat')
dep48TestX = eachTestX;
dep48TestY = eachTestY;
load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\chowTest.mat')
chowTestX = eachTestX;
chowTestY = eachTestY;
for ii = 1:240
   load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\each\500Trials_50-50\',num2str(ii),'fix.mat'])
   animal = ceil(ii/20);
   iter = rem(ii,20);
   if iter == 0
       iter = 20;
   end
   testX = [];
   testY = [];
   testX = dep24TestX{iter,animal};%(:,[1:2,4:34,36:55,57:60]);
   testY = dep24TestY{iter,animal};
   prob = cvglmnetPredict(selfData.model,zscore(testX),'lambda_1se','response');
   [dep24XEach{iter,animal},dep24YEach{iter,animal},~,dep24AEach(iter,animal)] = perfcurve(testY,prob,1);
   if animal <= 9
       testX = [];
       testY = [];
       testX = dep48TestX{iter,animal};%(:,[1:2,4:34,36:55,57:60]);
       testY = dep48TestY{iter,animal};
       prob = cvglmnetPredict(selfData.model,zscore(testX),'lambda_1se','response');
       [dep48XEach{iter,animal},dep48YEach{iter,animal},~,dep48AEach(iter,animal)] = perfcurve(testY,prob,1);
       testX = [];
       testY = [];
       testX = chowTestX{iter,animal};%(:,[1:2,4:34,36:55,57:60]);
       testY = chowTestY{iter,animal};
       prob = cvglmnetPredict(selfData.model,zscore(testX),'lambda_1se','response');
       [chowXEach{iter,animal},chowYEach{iter,animal},~,chowAEach(iter,animal)] = perfcurve(testY,prob,1);
   end
end
%% Prepare for plotting
[dep24xFillEach,dep24yFillEach] = avgFill(cat(2,dep24XEach{:}),cat(2,dep24YEach{:}),2);
dep24XMEach = mean(cat(2,dep24XEach{:}),2);
dep24YMEach = mean(cat(2,dep24YEach{:}),2);
[dep24ES,~] = distES(dep24A,reshape(dep24AEach,1,240));
[dep48xFillEach,dep48yFillEach] = avgFill(cat(2,dep48XEach{:}),cat(2,dep48YEach{:}),2);
dep48XMEach = mean(cat(2,dep48XEach{:}),2);
dep48YMEach = mean(cat(2,dep48YEach{:}),2);
[dep48ES,~] = distES(dep48A,reshape(dep48AEach,1,180));
[chowxFillEach,chowyFillEach] = avgFill(cat(2,chowXEach{:}),cat(2,chowYEach{:}),2);
chowXMEach = mean(cat(2,chowXEach{:}),2);
chowYMEach = mean(cat(2,chowYEach{:}),2);
[chowES,~] = distES(chowA,reshape(chowAEach,1,180));
%% Plot generalized and each together
figure
hold on
fill(dep24xFillEach,dep24yFillEach,[0.8 0.8 0.8])
h1 = plot(dep24XMEach,dep24YMEach,'--k','LineWidth',2);
fill(dep24Xfill,dep24Yfill,[0.8 0.8 0.8])
h2 = plot(dep24XM,dep24YM,'-k','LineWidth',2);
legend([h1,h2],{['Individualized: ',num2str(round(mean(reshape(dep24AEach,1,240)),2))],['Generalized: ',num2str(round(mean(dep24A),2))]},'Location','southeast')
text(0.8,0.2,['d = ',num2str(round(dep24ES,2))])
ylim([0 1])
xlim([0 1])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
title('Dep24: Generalized vs. Individualized')
xlabel('False Positive Rate')
ylabel('True Positive Rate')

figure
hold on
fill(dep48xFillEach,dep48yFillEach,[0.8 0.8 0.8])
h1 = plot(dep48XMEach,dep48YMEach,'--k','LineWidth',2);
fill(dep48Xfill,dep48Yfill,[0.8 0.8 0.8])
h2 = plot(dep48XM,dep48YM,'-k','LineWidth',2);
legend([h1,h2],{['Individualized: ',num2str(round(mean(reshape(dep48AEach,1,180)),2))],['Generalized: ',num2str(round(mean(dep48A),2))]},'Location','southeast')
text(0.8,0.2,['d = ',num2str(round(dep48ES,2))])
ylim([0 1])
xlim([0 1])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
title('Dep48: Generalized vs. Individualized')
xlabel('False Positive Rate')
ylabel('True Positive Rate')

figure
hold on
fill(chowxFillEach,chowyFillEach,[0.8 0.8 0.8])
h1 = plot(chowXMEach,chowYMEach,'--k','LineWidth',2);
fill(chowXfill,chowYfill,[0.8 0.8 0.8])
h2 = plot(chowXM,chowYM,'-k','LineWidth',2);
legend([h1,h2],{['Individualized: ',num2str(round(mean(reshape(chowAEach,1,180)),2))],['Generalized: ',num2str(round(mean(chowA),2))]},'Location','southeast')
text(0.8,0.2,['d = ',num2str(round(chowES,2))])
ylim([0 1])
xlim([0 1])
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
title('Chow: Generalized vs. Individualized')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
%% Build and test univariate logistic from baseline
% Concatenated data
% load('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat50-50_500Train.mat')
for ii = 1:20
    trainX = allTrainX{ii,1};%(:,[1:2,4:34,36:55,57:60]);
    trainY = allTrainY{ii,1};
    testX = allTestX{ii,1};%(:,[1:2,4:34,36:55,57:60]);
    testY = allTestY{ii,1};
    for bi = 1:60
        mdl{ii,bi} = fitglm(trainX(:,bi),trainY,'distribution','binomial');
        prob = predict(mdl{ii,bi},testX(:,bi));
        [~,~,~,catTestA(ii,bi)] = perfcurve(testY,prob,1);
%         for k = 1:11
%             testX = eachTestX{ii
%             prob = predict(mdl{ii,bi},eachTestX{ii,k}(:,bi));
%             [~,~,~,indTestA(ii,bi,k)] = perfcurve(eachTestY{ii,k},prob,1);
%         end
    end
end
%%
mTestA = mean(catTestA,1);
sTestA = std(catTestA,[],1);
[smTestA,inds] = sort(mTestA,'descend');
ssTestA = sTestA(inds);
scatterErr(1:60,smTestA,ssTestA,1)
hold on
% Plot 50 line
plot([0 60],[0.5 0.5],'--k','LineWidth',2)
% Plot line of generalized model at baseline (0.7966)
plot([0 60],[mean(ccA) mean(ccA)],'-k','LineWidth',2)
xlabel('Feature')
ylabel('AUC')
title('Average AUC from Univariate Logistic: Concat')
%% Test best model across times
best = inds(2);
for ii = 1:size(allTestX,1)
   for jj = 1:size(allTestX,2)
       prob = predict(mdl{ii,best},allTestX{ii,jj}(:,best));
       [~,~,~,a(ii,jj)] = perfcurve(allTestY{ii,jj},prob,1);
   end
end
%%
figure
scatterErr(0:60,mean(a,1),std(a,[],1),1)
set(gca,'XTickLabel',2.5:10:62.5)
xlabel('Time before Eating (sec)')
ylabel('AUC')
%% 
load('I6Base_2015-11-24.mat')
disp('Applying 60 Hz filter with filter60.m...')
[LFPTs.data] = filter60(LFPTs,adfreq,'off');
disp('Downsampling signal with dwnSample.m...')
[LFPTs,adfreq] = dwnSample(LFPTs,5,adfreq);
%%
% Calculate spectrogram
[s,f,t] = spectrogram(LFPTs.data(1,:),2048,1024,1:100,400);
% Convert s to dB
pow = 10*log10(abs(s));
% For each window (column of pow) get and normalize high gamma power
hgPowNorm = trapz(pow(70:90,:),1)./trapz(pow(1:100,1));
% Plot normalized high gamma power
figure
hold on
plot(t,hgPowNorm)
% Add binge events
for ii = 1:size(eventTs.t{1,7},1)
    plot([eventTs.t{1,7}(ii) eventTs.t{1,8}(ii)],[1 1],'-k','LineWidth',2)
end
%% logCmbs
for ii = 1:3
    disp(num2str(ii))
    load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\logCmbsFull3\',num2str(ii),'.mat'])
    baseAuc{ii} = A;
    dep24Auc{ii} = dep24A;
    dep48Auc{ii} = dep48A;
    chowAuc{ii} = chowA;
end
%%
baseM = cellfun(@(x) mean(x,1),baseAuc,'UniformOutput',0);
dep24M = cellfun(@(x) mean(x,1),dep24Auc,'UniformOutput',0);
dep48M = cellfun(@(x) mean(x,1),dep48Auc,'UniformOutput',0);
chowM = cellfun(@(x) mean(x,1),chowAuc,'UniformOutput',0);
%% Get top 5 models for each condition and cmbs
for ii  = 1:3
    cmbs{ii} = nchoosek(1:60,ii);
    [sortX,sortInd] = sort(baseM{1,ii},'descend');
    topA(:,ii) = sortX(1:60);
    for jj = 1:60
        topInd{jj,ii} = cmbs{ii}(sortInd(jj),:);
    end
end
nameVect = names({'SL','CL','SR','CR'},{'d','t','a','b','lg','hg'});
for ii = 1:60
    for jj = 1:2
        dyad{ii,jj} = nameVect{topInd{ii,2}(jj)};
    end
    for jj = 1:3
        triad{ii,jj} = nameVect{topInd{ii,3}(jj)};
    end
end
% for ii = 1:3
%    figure
%    hold on
%    scatter(1:length(baseM{ii}),sort(baseM{ii},'descend'),'.')
%    scatter(1:length(dep24M{ii}),sort(dep24M{ii},'descend'),'.')
%    scatter(1:length(dep48M{ii}),sort(dep48M{ii},'descend'),'.')
%    scatter(1:length(chowM{ii}),sort(chowM{ii},'descend'),'.')
% end
%%
[~,ord] = sort(baseM{3},'descend');
figure
subplot(2,2,1)
hold on
scatter(1:length(baseM{3}),sort(baseM{3},'descend'),'.')
scatter(1:length(dep24M{3}),sort(dep24M{3},'descend'),'.')
scatter(1:length(dep48M{3}),sort(dep48M{3},'descend'),'.')
scatter(1:length(chowM{3}),sort(chowM{3},'descend'),'.')

subplot(2,2,2)
hold on
scatter(1:34220,baseM{3}(ord),'.')
scatter(1:34220,dep24M{3}(ord),1,[0.85 0.325 0.098])
scatter(1:34220,sort(dep24M{3},'descend'),'.k')

subplot(2,2,3)
hold on
scatter(1:34220,baseM{3}(ord),'.')
scatter(1:34220,dep48M{3}(ord),1,[0.929 0.694 0.125])
scatter(1:34220,sort(dep48M{3},'descend'),'.k')

subplot(2,2,4)
hold on
scatter(1:34220,baseM{3}(ord),'.')
scatter(1:34220,chowM{3}(ord),1,[0.494 0.184 0.556])
scatter(1:34220,sort(chowM{3},'descend'),'.k')
%%
baseAM = cell2mat(cellfun(@(x) mean(mean(x,1)),baseAuc,'UniformOutput',0));
baseAS = cell2mat(cellfun(@(x) std(mean(x,1),[],2),baseAuc,'UniformOutput',0));
dep24AM = cell2mat(cellfun(@(x) mean(mean(x,1)),dep24Auc,'UniformOutput',0));
dep24AS = cell2mat(cellfun(@(x) std(mean(x,1),[],2),dep24Auc,'UniformOutput',0));
dep48AM = cell2mat(cellfun(@(x) mean(mean(x,1)),dep48Auc,'UniformOutput',0));
dep48AS = cell2mat(cellfun(@(x) std(mean(x,1),[],2),dep48Auc,'UniformOutput',0));
chowAM = cell2mat(cellfun(@(x) mean(mean(x,1)),chowAuc,'UniformOutput',0));
chowAS = cell2mat(cellfun(@(x) std(mean(x,1),[],2),chowAuc,'UniformOutput',0));
%% Prebinge data
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\preBinge2\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    preA(ii,:) = concatData{end}.auc;
end
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\preBinge2Rand\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    preARand(ii,:) = concatData{end}.auc;
end
scatterErr(0:60,mean(preA,1),std(preA,[],1),1)
%% Open preBinge data
for ii = 1:20
    load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\preBinge2\',num2str(ii),'.mat'])
    preA(ii,:) = concatData{61}.auc;
    load(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\preBinge2Rand\',num2str(ii),'.mat'])
    preRandA(ii,:) = concatData{61}.auc;
end
%%
figure
hold on
scatterErr(1:61,mean(preA,1),std(preA,[],1),0)
scatterErr(1:61,mean(preRandA,1),std(preRandA,[],1),0)
set(gca,'XTick',1:10:61,'XTickLabel',-2.5:-10:-62.5,'XDir','Reverse','YAxisLocation','right')
xlim([0 62.5])
xlabel('Time before Feeding (sec)')
ylabel('AUC')
%% Extract all shell left high gamma (ind = 6) from all test sets
preOne = cell(1,61); preZero = cell(1,61);
for ii = 1:size(allTestY,1)
    for jj = 1:size(allTestY,2)
        % Get inds of ones
        oneInd{ii,jj} = logicFind(1,allTestY{ii,jj},'==');
        % Get inds of zeros
        zeroInd{ii,jj} = logicFind(0,allTestY{ii,jj},'==');
        % Use above to extract and concatenate shell left high gamma
        preOne{jj} = [preOne{jj};allTestX{ii,jj}(oneInd{ii,jj},6)];
        preZero{jj} = [preZero{jj};allTestX{ii,jj}(zeroInd{ii,jj},6)];
    end
end
%% Average each time step
preAvg(1,:) = cellfun(@mean,preOne,'UniformOutput',0);
preAvg(2,:) = cellfun(@mean,preZero,'UniformOutput',0);
%%
% Extract binge start time using second column of samp
ind = cell(1);
for k = 1:size(samp{1,1},1)
    starts = [];
    for ii = 1:61
        starts = [starts;samp{1,1}{k,ii}(:,2)+(400*(ii-1))];
    end
    starts = unique(starts);
    for jj = 1:length(starts)
        for ii = 1:61
            thisInd = logicFind(starts(jj)-400*(ii-1),samp{1,1}{k,ii}(:,2),'==');
            if isempty(thisInd)
                ind{k}(jj,ii) = NaN;
            else
                ind{k}(jj,ii) = thisInd;
            end
        end
    end
end
%% Use indices to extract shell left high gamma from pre-trials
for ii = 1:23
    for jj = 1:61
        for k = 1:size(ind{ii},1)
            if isnan(ind{ii}(k,jj))
                hgPow{ii}(k,jj) = NaN;
            else
                hgPow{ii}(k,jj) = data{1,1}{ii,jj}(ind{ii}(k,jj),6);
            end
        end
    end
    hgPowRand{ii} = data{1,1}{ii,end}(:,6);
end
%%
test = diff(samp{1,1}{1,62}(:,2));
contig = [];
contig(1) = 1;
r = 2;
c = 1;
ii = 2;
while ii < size(samp{1,1}{1,62},1)
    if samp{1,1}{1,62}(ii,1) == samp{1,1}{1,62}(ii-1,1)+2000
        c = c+1;
        contig(r-1,c) = ii;
        ii = ii+1;
    else
        c = 1;
        r = r+1;
        contig(r-1,c) = ii;
        ii = ii+1;
    end
end
%%
allPrePow = cat(1,hgPow{:});
allRandPow = 
%%
% Get average ROC for concat-concat models
ccMX = mean(cat(2,ccX{:}),2);
ccMY = mean(cat(2,ccY{:}),2);
% Get average 'fill' for concat-concat models
[ccXfill,ccYfill] = avgFill(cat(2,ccX{:}),cat(2,ccY{:}),2,1);
%% Combine other conditions and plot
for ii = 1:57
   m(ii) = mean([reshape(dep24Auc{1,ii},1,numel(dep24Auc{1,ii})),reshape(dep48Auc{1,ii},1,numel(dep48Auc{1,ii})),reshape(chowAuc{1,ii},1,numel(chowAuc{1,ii}))]); 
   s(ii) = std([reshape(dep24Auc{1,ii},1,numel(dep24Auc{1,ii})),reshape(dep48Auc{1,ii},1,numel(dep48Auc{1,ii})),reshape(chowAuc{1,ii},1,numel(chowAuc{1,ii}))]); 
end
scatterErr(1:57,m,s,1)
%% Plot
figure
subplot(2,2,1)
scatterErr(1:57,baseAM,baseAS,0)
ylim([0.5 0.8])
title('Baseline-Baseline Model')
subplot(2,2,2)
scatterErr(1:57,dep24AM,dep24AS,0)
ylim([0.5 0.8])
title('Baseline-Dep24 Model')
subplot(2,2,3)
scatterErr(1:57,dep48AM,dep48AS,0)
ylim([0.5 0.8])
title('Baseline-Dep48 Model')
subplot(2,2,4)
scatterErr(1:57,chowAM,chowAS,0)
ylim([0.5 0.8])
title('Baseline-Chow Model')

%% Load and analyze given animals/binge epochs
fNames = {'I2BaseDec15.mat','I1BaseNov9.mat','I12BaseNov12.mat','H15BaseSep26.mat'};
trial = [1,35,5,8];
for ii = 1:4
    load(fNames{ii})
    start{ii} = nearest_idx3(eventTs.t{1,7}(trial(ii)),LFPTs.tvec);
    stop{ii} = nearest_idx3(eventTs.t{1,8}(trial(ii)),LFPTs.tvec);
    x1 = LFPTs.data(4,start{ii}-20000:start{ii}+20000);
    y1 = LFPTs.data(3,start{ii}-20000:start{ii}+20000);
    x1 = filter60(x1,adfreq,'off');
    y1 = filter60(y1,adfreq,'off');
    [pow{ii},f{ii},coi{ii}] = cwt(x1,'morse',adfreq);
    [wcoh{ii},~,cohf{ii}] = wcoherence(x1,y1,adfreq);
    %     [s{ii},f{ii},~] = spectrogram(LFPTs.data(4,start{ii}-20000:start{ii}+20000),512,256,1:200,adfreq);
    x2 = LFPTs.data(4,start{ii}-4000:start{ii}+4000);
    y2 = LFPTs.data(3,start{ii}-4000:start{ii}+4000);
    x2 = filter60(x2,adfreq,'off');
    y2 = filter60(y2,adfreq,'off');
    [pow2{ii},f2{ii},coi2{ii}] = cwt(x2,'morse',adfreq);
    [wcoh2{ii},~,cohf2{ii}] = wcoherence(x2,y2,adfreq);
    %     [s2{ii},f2{ii},~] = spectrogram(LFPTs.data(4,start{ii}-4000:start{ii}+4000),256,128,1:200,adfreq);
    t{ii} = (start{ii}-20000:start{ii}+20000)./adfreq;
    t2{ii} = (start{ii}-4000:start{ii}+4000)./adfreq;
end
%% Plot wavelet PSD
c = [0.06,0.08,0.07,0.24];
for ii = 1:4
    figure
    ftest = log2(f{ii});
%     h = pcolor(t{ii},ftest,abs(pow{ii}));
%     h.EdgeColor = 'none';
    imagesc(t{ii},ftest,abs(pow{ii}));
    set(gca,'ydir','normal','ytick',([0:7,7.6439]),'yticklabel',round(2.^([0:7,7.6439])))
    ylim([0 7.6439])
    caxis([0 c(ii)])
    colormap('viridis')
    colorbar
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    hold on
    line([start{ii}/adfreq start{ii}/adfreq],[0 7.6439],'Color','w','LineWidth',2)
    
    title([fNames{ii}])
    
    figure
    ftest = log2(f2{ii});
%     h = pcolor(t2{ii},ftest,abs(pow2{ii}));
%     h.EdgeColor = 'none';
    imagesc(t2{ii},ftest,abs(pow2{ii}));
    set(gca,'ydir','normal','ytick',([0:7,7.6439]),'yticklabel',round(2.^([0:7,7.6439])))
    ylim([0 7.6439])
    caxis([0 c(ii)])
    colormap('viridis')
    colorbar
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    hold on
    line([start{ii}/adfreq start{ii}/adfreq],[0 7.6439],'Color','w','LineWidth',2)
    title([fNames{ii}])
end
%% Plot wavelet COH
for ii = 1:4
    figure
    ftest = log2(cohf{ii});
%     h = pcolor(t{ii},ftest,wcoh{ii});
%     h.EdgeColor = 'none';
    imagesc(t{ii},ftest,wcoh{ii})
    set(gca,'ydir','normal','ytick',([0:7,7.6439]),'yticklabel',round(2.^([0:7,7.6439])))
    ylim([0 7.6439])
    colormap('viridis')
    colorbar
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    hold on
    line([start{ii}/adfreq start{ii}/adfreq],[0 7.6439],'Color','w','LineWidth',2)
    title([fNames{ii}])
    
    figure
    ftest = log2(cohf2{ii});
%     h = pcolor(t2{ii},ftest,wcoh2{ii});
%     h.EdgeColor = 'none';
    imagesc(t2{ii},ftest,wcoh2{ii})
    set(gca,'ydir','normal','ytick',([0:7,7.6439]),'yticklabel',round(2.^([0:7,7.6439])))
    ylim([0 7.6439])
    colormap('viridis')
    colorbar
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    hold on
    line([start{ii}/adfreq start{ii}/adfreq],[0 7.6439],'Color','w','LineWidth',2)
    title([fNames{ii}])
end
%% Subset Indices
% Power indices (channel,band)
powInd = reshape(1:24,6,4)';
% Coherence indices (cmb,band)
cohInd = reshape(25:60,6,6)';
% Delta indices
dInd = [powInd(:,1);cohInd(:,1)]';
% Theta indices
tInd = [powInd(:,2);cohInd(:,2)]';
% Alpha indices
aInd = [powInd(:,3);cohInd(:,3)]';
% Beta indices
bInd = [powInd(:,4);cohInd(:,4)]';
% Low gamma indices
lgInd = [powInd(:,5);cohInd(:,5)]';
% High gamma indices
hgInd = [powInd(:,6);cohInd(:,6)]';
% Left indices
lInd = [powInd([1,3],:);cohInd(2,:)];
% Right indices
rInd = [powInd([2,4],:);cohInd(5,:)];
% Shell indices
sInd = [powInd([1,2],:);cohInd(1,:)];
% Core indices
cInd = [powInd([3,4],:);cohInd(6,:)];
%
inds{1} = powInd;
inds{2} = cohInd;
inds{3} = dInd;
inds{4} = tInd;
inds{5} = aInd;
inds{6} = bInd;
inds{7} = lgInd;
inds{8} = hgInd;
inds{9} = lInd;
inds{10} = rInd;
inds{11} = sInd;
inds{12} = cInd;
%% Load subset files
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\subset\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    for jj = 1:12
       if isequal(concatData.rocX{jj},[0;1])
          concatData.rocX{jj} = (0:1/1500:1)';
          concatData.rocY{jj} = (0:1/1500:1)';
       end
    end
    ccSX(:,:,ii) = cat(2,concatData.rocX{:});
    ccSY(:,:,ii) = cat(2,concatData.rocY{:});
    ccSA(ii,:) = cell2mat(concatData.auc);
end
% Load subsetRand files
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\analyzed\concat\subsetRand\')
for ii = 1:20
    load([num2str(ii),'.mat'])
    for jj = 1:12
       if isequal(concatData.rocX{jj},[0;1])
          concatData.rocX{jj} = (0:1/1500:1)';
          concatData.rocY{jj} = (0:1/1500:1)';
       end
    end
    ccSRandX(:,:,ii) = cat(2,concatData.rocX{:});
    ccSRandY(:,:,ii) = cat(2,concatData.rocY{:});
    ccSRandA(ii,:) = cell2mat(concatData.auc);
end
% Get average fills
for ii = 1:12
    [ccSXfill(ii,:),ccSYfill(ii,:)] = avgFill(ccSX(:,ii,:),ccSY(:,ii,:),3);
    [ccSXrandFill(ii,:),ccSYrandFill(ii,:)] = avgFill(ccSRandX(:,ii,:),ccSRandY(:,ii,:),3);
end
% Get average subset AUC
mccSA = mean(ccSA,1);
sccSA = std(ccSA,[],1);
%% Plot each subset
figure
hold on
for ii = 1:12
   plot(mean(ccSX(:,ii,:),3),mean(ccSY(:,ii,:),3)) 
end
% Add concat-concat line
plot(ccMX,ccMY,'-k','LineWidth',2)
set(gca,'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Subset to Concat ROCs')
%%
figure
for ii = 1:12
    subplot(4,3,ii)
    hold on
    fill(ccXfill(ii,:),ccYfill(ii,:),'b')
    fill(ccXrandFill(ii,:),ccYrandFill(ii,:),'b')
    xlim([0 1])
    ylim([0 1])
end
%% Get effect size for each subset
for ii = 1:12
   [subsetES(ii),~] = distES(ccA(:,ii),ccRandA(:,ii)); 
end
