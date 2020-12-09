%% Single feature analysis - colormaps
c = viridis;
sites = {'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC','rNAcC'};
% subInds = [1:12,37:54,79:90,115:126,211:216];
% sites = {'lPFC','rPFC','lNAc','rNAc'};
freq = {'\delta','\theta','\alpha','\beta','l\gamma','h\gamma'};
chan = numel(sites);
thisData = aSalineSingle';
thisDataP = aLSDsingleP';
signA = sign(mean(betaLSD,2));
%%
thisDataM = mean(thisData,1);
[~,p] = ttest2(thisData,thisDataP);
pAdj = p.*numel(thisDataM);
thisDataM(pAdj>=0.05) = NaN;
signA(pAdj>=0.05) = NaN; 
powerN = chan*numel(freq);
cmbs = nchoosek(1:numel(sites),2);

aPowerMat = reshape(thisDataM(1:powerN).*signA(1:powerN)',numel(freq),chan)';
maxN = max([chan,numel(freq)]);
if any(signA(1:powerN)==-1)
    figure
    % padarray for visulization and to make squares
    pcolor(padarray(aPowerMat,[1+maxN-chan 1+maxN-numel(freq)],NaN,'post'))
    % Use lower half of colormap
    colormap(c(1:128,:))
    caxis([-0.8 -0.5])
    set(gca,'xtick',1.5:1:numel(freq)+0.5,'xticklabel',freq,'ytick',1.5:1:chan+0.5,'yticklabel',sites)
    axis square
    hold on
    [row,col] = ind2sub([chan,numel(freq)],logicFind(1,aPowerMat<0,'=='));
    viscircles([col'+0.5,row'+0.5],repmat(0.425,1,numel(row)),'color','none','linewidth',0.1);
end

if any(signA(1:powerN)==1)
    figure
    pcolor(padarray(aPowerMat,[1+maxN-chan 1+maxN-numel(freq)],NaN,'post'))
    colormap(c(129:end,:))
    caxis([0.5 0.8])
    set(gca,'xtick',1.5:1:numel(freq)+0.5,'xticklabel',freq,'ytick',1.5:1:chan+0.5,'yticklabel',sites)
    axis square
    hold on
    [row,col] = ind2sub([chan,numel(freq)],logicFind(1,aPowerMat>0,'=='));
    viscircles([col'+0.5,row'+0.5],repmat(0.425,1,numel(row)),'color','none','linewidth',0.1);
end

aCohMat = thisDataM;
aCohMat(aCohMat<=0.5) = NaN;

siteCmb = cell(1,size(cmbs,1));
for ii = 1:size(cmbs,1)
    siteCmb{ii} = [sites{cmbs(ii,1)},'-',sites{cmbs(ii,2)}];
end
aCohMat = reshape(aCohMat(powerN+1:end).*signA(powerN+1:end)',numel(freq),numel(siteCmb))';
aCohMat = aCohMat([1:3,8,9,14,4:7,10:13,15:end],:);
siteCmb = siteCmb([1:3,8,9,14,4:7,10:13,15:end]);
maxN = max([numel(siteCmb),numel(freq)]);
if any (signA(powerN+1:end)==-1)
    figure
    pcolor(padarray(aCohMat,[1+maxN-numel(siteCmb) 1+maxN-numel(freq)],NaN,'post'))
    colormap(c(1:128,:))
    caxis([-0.8 -0.5])
    set(gca,'xtick',1.5:1:numel(freq)+0.5,'xticklabel',freq,'ytick',1.5:1:numel(siteCmb)+0.5,'yticklabel',siteCmb)
    axis square
    hold on
    [row,col] = ind2sub([numel(siteCmb),numel(freq)],logicFind(1,aCohMat<0,'=='));
    viscircles([col'+0.5,row'+0.5],repmat(0.425,1,numel(row)),'color','none','linewidth',0.1)
end
if any(signA(powerN+1:end)==1)
    figure
    pcolor(padarray(aCohMat,[1+maxN-numel(siteCmb) 1+maxN-numel(freq)],NaN,'post'))
    colormap(c(129:end,:))
    caxis([0.5 0.8])
    set(gca,'xtick',1.5:1:numel(freq)+0.5,'xticklabel',freq,'ytick',1.5:1:numel(siteCmb)+0.5,'yticklabel',siteCmb)
    axis square
    hold on
    [row,col] = ind2sub([numel(siteCmb),numel(freq)],logicFind(1,aCohMat>0,'=='));
    viscircles([col'+0.5,row'+0.5],repmat(0.425,1,numel(row)),'color','none','linewidth',0.1)
end
sum(sum(~isnan(aPowerMat)))/numel(aPowerMat)

sum(sum(~isnan(aCohMat)))/numel(aCohMat)

[~,p,chi2stat,df] = prop_test([sum(sum(~isnan(aPowerMat))),sum(sum(~isnan(aCohMat)))],[numel(aPowerMat),numel(aCohMat)],'true')