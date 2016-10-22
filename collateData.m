function [] = collateData(sdir,searchStr)
%%
% Get file structure for each searchStr
[fileStruct] = fileSearch(sdir,searchStr);
% Go through each file structure and extract data from files therein
for fsi = 1:numel(fileStruct)
    for fi = 1:size(fileStruct{fsi})
        load([sdir{1},fileStruct{fsi}(fi).name],'coh','psdTrls','relPower','powerCorrSort','LFPTs');
        % Stack overall PSD and PSD distributions
        PSD{fsi}(:,:,fi) = psdTrls.event1.Overall;
        relPSD{fsi}(:,:,fi) = relPower;
        powMat{fsi,fi} = [];
        for ii = 1:length(psdTrls.event1.Pow)
            powMat{fsi,fi} = cat(3,powMat{fsi,fi},psdTrls.event1.Pow{1,ii});
        end
        % Stack mean coherence
        meanCoh{fsi}(:,:,fi) = coh.mCxy;
        cohMat{fsi,fi} = coh.Cxy;
        % Stack power correlations
        powCorr{fsi}(:,:,fi) = powerCorrSort.mastercorrTWPlot;
        % Stack number of windows
        n{fsi}(fi) = size(psdTrls.event1.Pow,2);
        % On last iteration grab channel labels
        if fsi == numel(fileStruct) && fi == size(fileStruct{fsi},2)
            chanLabel = LFPTs.label;
        end
        clearvars -except n powCorr meanCoh relPSD PSD cohMat powMat chanLabel fsi fi fileStruct sdir searchStr
    end
end
%% Get weighted means and standard deviations
wAvgPSD = wmean(PSD,n);
wAvgRelPSD = wmean(relPSD,n);
wAvgCoh = wmean(meanCoh,n);
wAvgPowCorr = wmean(powCorr,n);
% Get weighted standard deviations
for wi = 1:numel(n)
   wSdPSD{wi} = std(PSD{wi},n{wi},3);
   wSdRelPSD{wi} = std(relPSD{wi},n{wi},3);
   wSdCoh{wi} = std(meanCoh{wi},n{wi},3);
   wSdPowCorr{wi} = std(powCorr{wi},n{wi},3);
end
%% Plot
col = {'-r','-b','-g','-k'};
figure
for ii = 1:numel(wAvgPSD)
    subplot(2,2,ii)
    hold on;
    for j = 1:size(wAvgPSD{ii},1)
        shadedErrorBar([1:2:100],wAvgPSD{ii}(j,:),wSdPSD{ii}(j,:),col{j},1)
    end
end
%% Plot PSD per channel across conditions
col = {'-r','-b','-g','-k'};
figure
for ii = 1:size(wAvgPSD{1},1) 
   subplot(2,2,ii)
   hold on;
   for j = 1:numel(wAvgPSD)
      h{ii,j} = shadedErrorBar((1:2:100),wAvgPSD{j}(ii,:),wSdPSD{j}(ii,:),col{j},1); 
   end
   xlabel('Frequency (Hz)'); ylabel('Power (dB)');
   title(chanLabel{ii})
   if ii == 1
       legend([h{ii,1}.patch,h{ii,2}.patch,h{ii,3}.patch,h{ii,4}.patch],'Base','Dep24','Dep48','Chow')
   end
    
end
%% Plot Coh per channel across conditions
cmb = nchoosek(1:4,2);
for c = 1:size(cmb,1)
    channelCmb(c,:) = chanLabel(cmb(c,:));
end
figure
for ii = 1:size(wAvgCoh{1},1)
    subplot(2,3,ii)
    hold on;
    for j = 1:numel(wAvgCoh)
        h2{ii,j} = shadedErrorBar((1:2:100),wAvgCoh{j}(ii,:),wSdCoh{j}(ii,:),col{j},1);
    end
    xlabel('Frequency (Hz)'); ylabel('Coherence Index')
    title(strcat(channelCmb(ii,1),'-',channelCmb(ii,2)))
    if ii == 1
        legend([h2{ii,1}.patch,h2{ii,2}.patch,h2{ii,3}.patch,h2{ii,4}.patch],'Base','Dep24','Dep48','Chow')
    end
end
%% Plot PSD distributions
% Take out first column of each
figure
for ii = 1:4
    subplot(2,2,ii)
    bar([wAvgRelPSD{1}(:,ii),wAvgRelPSD{2}(:,ii),wAvgRelPSD{3}(:,ii),wAvgRelPSD{4}(:,ii)])
    title(chanLabel{ii});
    ylabel('Proportion of Power');
    set(gca,'XTick',1:5,'XTickLabel',{'\theta','\alpha','\beta','l\gamma','h\gamma'});
    xlabel('Frequency Band');
    if ii == 1
        legend('Base','Dep24','Dep48','Chow','Location','northwest')
    end
end
%% Plot PowerCorr
titles = {'Base','Dep24','Dep48','Chow'};
for ii = 1:4
   figure
   subplot(1,2,1)
   pcolor(wAvgPowCorr{ii})
   colormap jet
   title(['Average Power Correlation ',titles{ii}])
   subplot(1,2,2)
   pcolor(wSdPowCorr{ii})
   colormap jet
   title(['Power Correlation Standard Deviation ',titles{ii}])
end
