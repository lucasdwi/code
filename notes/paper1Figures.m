%% Plot real and random histograms
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\cvMeanBWHist.mat')
subtitles = {'Shell All','Shell Strict','Core All','Core Strict'};
figure;
for c = 1:4
    subplot(2,2,c)
    histogram(realData.allAcc(:,c),'FaceColor','k','FaceAlpha',0.5,'Normalization','probability','BinWidth',0.02)
    title(subtitles{c})
    hold on
    histogram(randData.allAcc(:,c),'FaceColor','w','FaceAlpha',0.5,'Normalization','probability','BinWidth',0.02)
    xlim([0.2 1]); ylim([0 .25])
end
% Perform 2 sample KS test
for c = 1:4
    [h(c), p(c)] = kstest2(realData.allErr(:,c),randData.allErr(:,c));
end
%% Plot KS Curves (cdfs)
figure
for r = 1:4
    subplot(2,2,r)
    ecdf(realData.allAcc(:,r))
    hold on
    ecdf(randData.allAcc(:,r))
    h = get(gca,'children');
    % Change real data to black
    set(h(2,1),'Color','k','LineWidth',1);
    % Change rand data to grey
    set(h(1,1),'Color',[0.5 0.5 0.5],'LineWidth',1)
    xlabel('Accuracy'); ylabel('Cumulative Density')
    title(subtitles{r})
    xlim([0.2 1])
end
%% Response site CDFs
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\responseSiteRealRand.mat')
figure
ecdf(realData.acc)
hold on
ecdf(randData.acc)
legend({'Real','Rand'},'location','southeast')
title('Response Site Prediction Accuracy: Real vs. Random Data')
xlabel('Accuracy')
ylabel('Cumulative Density')
%% Response site PCA AUC
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\responseSitePCA.mat')
figure
plot(xMat{1,1},yMat{1,1})
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('AUC of Response Site Prediction from Top PCs')

%% Response Site PCA 3D
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\responseSitePCA.mat')
figure
d1 = scores(:,1); d2 = scores(:,2); d3 = scores(:,3);
inds1 = logicFind(1,y2,'==');
inds0 = logicFind(0,y2,'==');
scatter3(d1(inds1,:),d2(inds1,:),d3(inds1,:),100,'r','MarkerFaceColor','r')
hold on
scatter3(d1(inds0,:),d2(inds0,:),d3(inds0,:),100,'k','MarkerFaceColor','k')
xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3')
%% Response Site Lasso
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\responseSite.mat')
inds1 = logicFind(1,y2,'==');
inds0 = logicFind(0,y2,'==');
figure
plot(x2(inds1,1).*100,x2(inds1,36).*100,'.r','MarkerSize',20)
hold on
plot(x2(inds0,1).*100,x2(inds0,36).*100,'.k','MarkerSize',20)
xlabel('Normalized \theta Power SL')
ylabel('Normalized \theta Coherence SR-CL')
title('\theta Power vs. \theta Coherence: Response Site')
%% Plot PCA Regresssion Results
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\pcaReg3.mat');
%%
figure
% Core
plot(xMat{1,2}{1,2},yMat{1,2}{1,2},'--k','linewidth',1.5)
hold on
% Shell
plot(xMat{1,3}{1,2},yMat{1,3}{1,2},':k','linewidth',1.5)
% All
plot(xMat{1,1}{1,2},yMat{1,1}{1,2},'k','linewidth',1.5)
legend({'Core','Shell','All'},'Location','southeast')
xlabel('False Positive Rate'); ylabel('True Positive Rate')
title('PCA Regression: Core Model')
%% Plot 3D PCA
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\pca3D.mat')
% Scores data obtained by using pcaNotes with ii = 1 to keep data from pca
% on first dataset (all features)
% Plots core responders as black, and non as red
figure
scatter3(d1(inds1),d2(inds1),d3(inds1),[],'k')
hold on
scatter3(d1(inds0),d2(inds0),d3(inds0),[],'r')
xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3')
title('Core Responders Separation by PCA')
%% Run behavioral t-tests
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\5behaviorResponseTable12animal.mat')
behav = table2array(T.Base(:,5:8)); group = table2array(T.Base(:,1:4));
% Go through each behavior (column of behav)
for c = 1:size(behav,2)
    % For each split data into resp vs. notResp
    for r = 1:size(group,2)
        resp = behav(group(:,r)==1,c); notResp = behav(group(:,r)==0,c);
        [~,p(c,r)] = ttest2(resp,notResp);
    end
end
%% Plot percent time at rest Core
basicStripPlot(behav(:,4).*100,group(:,3),'l')
xtick([1 2])
set(gca,'xticklabels',{'R','NR'})
xlabel('Group')
ylabel('% Time Spent at Rest')
sigstar({[1 2]},p(4,3))
title('Time Spent at Rest vs. Core Response')
%% Plot percent time spent at rest Shell
basicStripPlot(behav(:,4).*100,group(:,1),'l')
xtick([1 2])
set(gca,'xticklabels',{'R','NR'})
xlabel('Group')
ylabel('% Time Spent at Rest')
sigstar({[1 2]},p(4,1))
title('Time Spent at Rest vs. Shell Response')
%% Stability
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\paper1data.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\cutoff40Betas.mat')
% Setup day information and perform subtractions
% Days between recordings of H10, H13, H14, H15, I11, I12, I1, I2, I3, I4, I6, I8
day = [18,46,16,12,1,13,14,34,2,71,67,36];
% Get every other row index for subtraction
first = [1:2:size(x,1)];
second = [2:2:size(x,1)];
for s = 1:length(first)
    sub(s,:) = (x(second(s),:) - x(first(s),:))/x(first(s),:);
end
absSub = abs(sub);
% Sum across rows and logicFind to get all betas that are in at least one model
inds = logicFind(0,sum(threshBeta,1),'~=');
% Go through all indices and correlate against time
for b = 1:length(inds)
    thisMd = fitlm(day,sub(:,b));
    ps(b) = thisMd.Coefficients.pValue(2);
end
% Apply FDR
[~,~,adjP] = fdr_bh(ps);
% Get number of significant linear regressions
sigCount = sum(ps <= 0.05);
adjSigCount = sum(adjP <= 0.05);
%% Theta CL-SR Stability
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\paper1data.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\cutoff40Betas.mat')
day = [18,46,16,12,1,13,14,34,2,71,67,36];
% Get every other row index for subtraction
first = [1:2:size(x,1)];
second = [2:2:size(x,1)];
for s = 1:length(first)
    sub(s,:) = (x(second(s),:) - x(first(s),:))./x(first(s),:);
end
md1 = fitlm(day([1:4,6:7,9:end]),sub([1:4,6:7,9:end],36));
% Get average shell and core responder difference
shellInd = [2,4,5,7,8,9,10,12];
shellAvg = mean(sub(shellInd,36));
coreInd = [1,3,11];
coreAvg = mean(sub(coreInd,36));
% Shell-Core/Shell
dif = (shellAvg-coreAvg)/coreAvg;
% Plot stability with responder difference
figure
scatter(day([1:4,6:7,9:end]),sub([1:4,6:7,9:end],36).*100,25,'k','filled')
hold on
% Plot linear regression
lsline
% Plot responder difference on both sides
plot(0:80,repmat(dif,1,81).*100,'-k')
plot(0:80,repmat(-dif,1,81).*100,'-k')
ylim([-100 100])
title('Theta Coherence Stability: Core Left-Shell Right')
xlabel('Days')
ylabel('%\Delta in Coherence')
%% Get significant feature counts
load('cutoff40Betas.mat')
nameVect = names;
% Only keep first 50 features
nameVect = nameVect(1:50);
% First get all count
letters = cell(4,40);
for r = 1:4
    these = logicFind(0,threshBeta(r,:),'~=');
    theseName = nameVect(these);
    C = cellfun(@(s) strsplit(s,'-'), theseName,'UniformOutput',false);
    for ci = 1:length(C)
       letters{r,ci} = C{ci}(2); 
    end
    ci = 1;
    for this = {'t','a','b','lg','hg'}
        count(:,ci) = sum(cellfun(@(sc) strcmpi(sc,this),letters),2);
        ci = ci +1;
    end
end
allcount = count;
all = sum(allcount);
%
tot = sum(allcount,2);
for r = 1:4
    norm(r,:) = allcount(r,:)./tot(r);
end
figure
bar(norm,'stacked'); title('% of Model per Feature')
% Get power
letters = cell(4,40);
for r = 1:4
    these = logicFind(0,threshBeta(r,1:20),'~=');
    theseName = nameVect(these);
    C = cellfun(@(s) strsplit(s,'-'), theseName,'UniformOutput',false);
    for ci = 1:length(C)
       letters{r,ci} = C{ci}(2); 
    end
    ci = 1;
    for this = {'t','a','b','lg','hg'}
        count(:,ci) = sum(cellfun(@(sc) strcmpi(sc,this),letters),2);
        ci = ci +1;
    end
end
powcount = count';
pow = sum(powcount);
% Get coh
letters = cell(4,40);
for r = 1:4
    these = logicFind(0,threshBeta(r,21:50),'~=');
    theseName = nameVect(these);
    C = cellfun(@(s) strsplit(s,'-'), theseName,'UniformOutput',false);
    for ci = 1:length(C)
       letters{r,ci} = C{ci}(2); 
    end
    ci = 1;
    for this = {'t','a','b','lg','hg'}
        count(:,ci) = sum(cellfun(@(sc) strcmpi(sc,this),letters),2);
        ci = ci +1;
    end
end
cohcount = count';
coh = sum(cohcount);

% Interleave pow and coh
allFeat = reshape([powcount(:) cohcount(:)]',2*size(powcount,1),[])';
% Normalize
allFeatNorm = diag(1./sum(allFeat,2))*allFeat;

