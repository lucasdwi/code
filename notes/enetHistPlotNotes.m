%% Plot real and random histograms
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\cvMeanBSHist.mat')
subtitles = {'Shell All','Shell Strict','Core All','Core Strict'};
figure;
for c = 1:4
    subplot(2,2,c)
    histogram(realErr.allAcc(:,c),'FaceColor','k','FaceAlpha',0.5,'Normalization','probability','BinWidth',0.02)
    title(subtitles{c})
    hold on
    histogram(randErr.allAcc(:,c),'FaceColor','w','FaceAlpha',0.5,'Normalization','probability','BinWidth',0.02)
    xlim([0.2 1]); ylim([0 .25])
end
% Perform 2 sample KS test
for c = 1:4
    [h(c), p(c)] = kstest2(realData.allErr(:,c),randData.allErr(:,c));
end
% Plot KS Curves (cdfs)
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
    sub(s,:) = x(second(s),:) - x(first(s),:);
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

