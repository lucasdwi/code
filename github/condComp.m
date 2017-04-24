%% For each coh and power value, get difference from baseline (with and without normalization from rest)
baseFiles = {'H10BaseOct15';'H13BaseSep27';'H14BaseOct15';'H15BaseSep29';'I1BaseNov9';'I2BaseDec15';'I3BaseNov11';'I4BaseSep24';'I6BaseNov24';'I8BaseOct29';'I11BaseOct30';'I12BaseNov12'};
dep24Files = {'H10FoodDep24Sep25';'H13FoodDep24Sep25';'H14FoodDep24Sep28';'H15FoodDep24Sep28';'I1FoodDep24Nov13';'I2FoodDep24Dec16';'I3FoodDep24Nov25';'I4FoodDep24Sep30';'I6FoodDep24Sep30';'I8FoodDep24Sep30';'I11FoodDep24Dec16';'I12FoodDep24Nov13'};
baseFiles2 = {'H10BaseSep27';'H13BaseNov12';'H14BaseSep29';'H15BaseOct11';'I1BaseOct26';'I2BaseNov9';'I3BaseSep24';'I4BaseDec4';;'I6BaseSep18';'I8BaseSep23';'I11BaseOct29';'I12BaseOct30'};
% dep48Files = 
% chowFiles = 
for fi = 1:numel(baseFiles)
   data{fi,1} = baseFiles{fi};
end

%% Just grab last n+1 rest sessions
% For just last trial, use n = 0
n = 0;
for fi = 1:length(baseFiles)
    load(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisRestTest\',baseFiles{fi},'_rest.mat'],'psdTrls','coh')
    lastRestPsd{fi} = psdTrls.event1.Pow(:,end-n:end);
    lastRestCoh{fi} = coh.Cxy(:,:,end-n:end);
end
%%
allLastRestPsd = []; allLastRestCoh = [];
for ii = 1:numel(lastRestCoh)
    %allLastRestPsd(:,:,ii) = lastRestPsd{ii}{1,1};
    allLastRestPsd = cat(3,allLastRestPsd,lastRestPsd{ii}{1,1});
    allLastRestCoh = cat(3,allLastRestCoh,lastRestCoh{ii});
end

%% Grab lots 'o data
for fi = 1:length(baseFiles)
    load([sdir,baseFiles{fi},'_binge.mat'],'psdTrls','coh','relPower')
    basePsd = psdTrls.event1.Overall;
    baseRelPower = relPower;
    baseCoh = coh.mCxy;
    baseRelCoh = coh.rel;
    load([sdir,dep24Files{fi},'_binge.mat'],'psdTrls','coh','relPower')
    dep24Psd = psdTrls.event1.Overall;
    dep24RelPower = relPower;
    dep24Coh = coh.mCxy;
    dep24RelCoh = coh.rel;
    changePsd(:,:,fi) = (dep24Psd-basePsd)./basePsd;
    changeCoh(:,:,fi) = (dep24Coh-baseCoh)./baseCoh;
    changeRelPower(:,:,fi) = (dep24RelPower-baseRelPower)./baseRelPower;
    changeRelCoh(:,:,fi) = (dep24RelCoh-baseRelCoh)./baseRelCoh;
end
%% Plot PSD difference across animals per channel
figure
for ii = 1:4
   subplot(2,2,ii)
   plot(psdTrls.F,sq(changePsd(ii,:,:).*100))
   xlabel('Frequency (Hz)'); ylabel('% Change in Power')
end
mtit('Percent Change in Resting Power from Baseline to Food Dep 24')
%% Plot change in coherence
figure;
for ii = 1:6
   subplot(2,3,ii)
   plot(psdTrls.F(2:end),sq(changeCoh(ii,2:end,:).*100))
   xlabel('Frequency (Hz)'); ylabel('% Change in Coherence')
end
mtit('Percent Change in Resting Coherence from Baseline to Food Dep 24')
%% Plot relative change in power
figure;
for ii = 1:4
   subplot(2,2,ii)
   bar(sq(changeRelPower(:,ii,:)).*100)
   title(['Channel ',num2str(ii)])
   xlabel('Band'); ylabel('Change in Relative Power (%)')
end
mtit('Change in Relative Power during Binge')
%% Plot relative change in coherence
cmb = nchoosek(1:4,2);
figure;
for ii = 1:6
   subplot(2,3,ii)
   bar(sq(changeRelCoh(:,ii,:)).*100)
   title(['Channel ',num2str(cmb(ii,1)),'-',num2str(cmb(ii,2))])
   xlabel('Band'); ylabel('Change in Relative Coherence (%)')
end
mtit('Change in Relative Coherence during Binge','xoff',0.1)
%% Grab relative PSD and Coh from rest and binge of baseline data
for fi = 1:length(baseFiles)
    % Load rest
    load(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisRestTest\',dep24Files{fi},'_rest.mat'],'psdTrls','coh','relPower')
    %restPsd = psdTrls.event1.Overall;
    restRelPower = relPower;
    %restCoh = coh.mCxy;
    restRelCoh = coh.rel;
    % Load binge
    load(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisBingeTest\',dep24Files{fi},'_binge.mat'],'psdTrls','coh','relPower')
    %bingePsd = psdTrls.event1.Overall;
    bingeRelPower = relPower;
    %bingeCoh = coh.mCxy;
    bingeRelCoh = coh.rel;
    % Get difference
    %changePsd(:,:,fi) = (bingePsd-restPsd)./restPsd;
    %changeCoh(:,:,fi) = (bingeCoh-restCoh)./restCoh;
    changeRelPower(:,:,fi) = (bingeRelPower-restRelPower)./restRelPower;
    changeRelCoh(:,:,fi) = (bingeRelCoh-restRelCoh)./restRelCoh;
end
%% Get average across animals
avgPowerChange = mean(changeRelPower,3);
sdPowerChange = std(changeRelPower,[],3);
avgCohChange = mean(changeRelCoh,3);
sdCohChange = std(changeRelCoh,[],3);
% Run t-tests
powT = zeros(size(avgPowerChange));
for ti = 1:numel(avgPowerChange)
    % Get row/col indices
    [r,c] = ind2sub([size(changeRelPower,1),size(changeRelPower,2)],ti);
    [~,p] = ttest(changeRelPower(r,c,:));
    powT(ti) = p;
end
cohT = zeros(size(avgCohChange));
for ti = 1:numel(avgCohChange)
    [r,c] = ind2sub([size(changeRelCoh,1),size(changeRelCoh,2)],ti);
    [~,p] = ttest(changeRelCoh(r,c,:));
    cohT(ti) = p;
end
%% Grab relative PSD and Coh from rest and binge of baseline data
for fi = 1:length(baseFiles)
    % Load rest
    load(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisRestTest\',baseFiles{fi},'_rest.mat'],'psdTrls','coh','relPower')
    %restPsd = psdTrls.event1.Overall;
    restRelPower1 = relPower;
    %restCoh = coh.mCxy;
    restRelCoh1 = coh.rel;
    % Load binge
    load(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisRestTest\',baseFiles2{fi},'_rest.mat'],'psdTrls','coh','relPower')
    %bingePsd = psdTrls.event1.Overall;
    restRelPower2 = relPower;
    %bingeCoh = coh.mCxy;
    restRelCoh2 = coh.rel;
    % Get difference
    %changePsd(:,:,fi) = (bingePsd-restPsd)./restPsd;
    %changeCoh(:,:,fi) = (bingeCoh-restCoh)./restCoh;
    changeRelPower(:,:,fi) = (restRelPower2-restRelPower1)./restRelPower1;
    changeRelCoh(:,:,fi) = (restRelCoh2-restRelCoh1)./restRelCoh1;
end
%% Get average across animals
avgPowerChange = mean(changeRelPower,3);
sdPowerChange = std(changeRelPower,[],3);
avgCohChange = mean(changeRelCoh,3);
sdCohChange = std(changeRelCoh,[],3);
% Run t-tests
powT = zeros(size(avgPowerChange));
for ti = 1:numel(avgPowerChange)
    % Get row/col indices
    [r,c] = ind2sub([size(changeRelPower,1),size(changeRelPower,2)],ti);
    [~,p] = ttest(changeRelPower(r,c,:));
    powT(ti) = p;
end
cohT = zeros(size(avgCohChange));
for ti = 1:numel(avgCohChange)
    [r,c] = ind2sub([size(changeRelCoh,1),size(changeRelCoh,2)],ti);
    [~,p] = ttest(changeRelCoh(r,c,:));
    cohT(ti) = p;
end
figure; barwitherr(sdPowerChange',avgPowerChange')
title('Change in Rest Power: 2 Baselines')

figure; barwitherr(sdCohChange',avgCohChange')
title('Change in Rest Coherence: 2 Baselines')
%% Grab relative PSD and Coh from rest and binge of baseline data
for fi = 1:length(baseFiles)
    % Load rest
    load(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisBingeTest\',baseFiles{fi},'_binge.mat'],'psdTrls','coh','relPower')
    %restPsd = psdTrls.event1.Overall;
    restRelPower1 = relPower;
    %restCoh = coh.mCxy;
    restRelCoh1 = coh.rel;
    % Load binge
    load(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisBingeTest\',baseFiles2{fi},'_binge.mat'],'psdTrls','coh','relPower')
    %bingePsd = psdTrls.event1.Overall;
    restRelPower2 = relPower;
    %bingeCoh = coh.mCxy;
    restRelCoh2 = coh.rel;
    % Get difference
    %changePsd(:,:,fi) = (bingePsd-restPsd)./restPsd;
    %changeCoh(:,:,fi) = (bingeCoh-restCoh)./restCoh;
    changeRelPower(:,:,fi) = (restRelPower2-restRelPower1)./restRelPower1;
    changeRelCoh(:,:,fi) = (restRelCoh2-restRelCoh1)./restRelCoh1;
end
% Get average across animals
avgPowerChange = mean(changeRelPower,3);
sdPowerChange = std(changeRelPower,[],3);
avgCohChange = mean(changeRelCoh,3);
sdCohChange = std(changeRelCoh,[],3);
% Run t-tests
powT = zeros(size(avgPowerChange));
for ti = 1:numel(avgPowerChange)
    % Get row/col indices
    [r,c] = ind2sub([size(changeRelPower,1),size(changeRelPower,2)],ti);
    [~,p] = ttest(changeRelPower(r,c,:));
    powT(ti) = p;
end
cohT = zeros(size(avgCohChange));
for ti = 1:numel(avgCohChange)
    [r,c] = ind2sub([size(changeRelCoh,1),size(changeRelCoh,2)],ti);
    [~,p] = ttest(changeRelCoh(r,c,:));
    cohT(ti) = p;
end
figure; barwitherr(sdPowerChange',avgPowerChange')
title('Change in Binge Power: 2 Baselines')

figure; barwitherr(sdCohChange',avgCohChange')
title('Change in Binge Coherence: 2 Baselines')
%% Grab relative PSD and Coh from rest and binge of baseline data
% Finds percent difference in 2 from 1
for fi = 1:length(baseFiles)
    % Load 1
    load(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisRestTest\',baseFiles{fi},'_rest.mat'],'psdTrls','coh','relPower')
    relPower1 = relPower;
    relCoh1 = coh.rel;
    % Load 2
    load(['C:\Users\Lucas\Desktop\GreenLab\data\fullAnalysisBingeTest\',baseFiles{fi},'_binge.mat'],'psdTrls','coh','relPower')
    relPower2 = relPower;
    relCoh2 = coh.rel;
    % Get difference
    changeRelPower(:,:,fi) = (relPower2-relPower1)./relPower1;
    changeRelCoh(:,:,fi) = (relCoh2-relCoh1)./relCoh1;
end
% Get average across animals
avgPowerChange = mean(changeRelPower,3);
sdPowerChange = std(changeRelPower,[],3);
avgCohChange = mean(changeRelCoh,3);
sdCohChange = std(changeRelCoh,[],3);
% Run t-tests
powT = zeros(size(avgPowerChange));
for ti = 1:numel(avgPowerChange)
    % Get row/col indices
    [r,c] = ind2sub([size(changeRelPower,1),size(changeRelPower,2)],ti);
    [~,p] = ttest(changeRelPower(r,c,:));
    powT(ti) = p;
end
cohT = zeros(size(avgCohChange));
for ti = 1:numel(avgCohChange)
    [r,c] = ind2sub([size(changeRelCoh,1),size(changeRelCoh,2)],ti);
    [~,p] = ttest(changeRelCoh(r,c,:));
    cohT(ti) = p;
end
figure; b1 = barwitherr(sdPowerChange',avgPowerChange');
title('Change in Power: Binge-Rest (Base)')

figure; b2 = barwitherr(sdCohChange',avgCohChange');
title('Change in Coherence: Binge-Rest (Base)')
%%
for ti = 1:numel(changeRelPower(:,:,1))
    [r,c] = ind2sub([size(changeRelPower,1),size(changeRelPower,2)],ti);
    [~,p] = ttest(pow1(r,c,:),pow2(r,c,:));
    powComp(r,c) = p;
end
for ti = 1:numel(changeRelCoh(:,:,1))
    [r,c] = ind2sub([size(changeRelCoh,1),size(changeRelCoh,2)],ti);
    [~,p] = ttest(coh1(r,c,:),coh2(r,c,:));
    cohComp(r,c) = p;
end