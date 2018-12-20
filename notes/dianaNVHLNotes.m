%% Load data - Norm
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\dianaNVHL\processed\'],{'sham';'nvhl'},{'pow','coh'},'avg',...
    'rel');
% [data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
%     'data\dianaNVHL\processed\'],{'sham';'nvhl'},{'pow','coh'},'trl',...
%     'rel');
%% Load data - raw
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\dianaNVHL\processed\'],{'sham';'nvhl'},{'pow','coh'},'avg','');
%% Concatenate data
shamInd = [1:5,7:11,12:16,17:21];
nvhlInd = [1:5,8:12,13:17,19:23];
sham = cat(1,data{1}{shamInd});
nvhl = cat(1,data{2}{nvhlInd});
%% All
sham = cat(1,data{1}{:});
nvhl = cat(1,data{2}{:});
%% Run t-tests
[h,p,pAdj] = bulkT(sham,nvhl,0,'fdr');
%%
nameVect = names({'lPL','lOFC','lCC','lNAc','rCC','rNAc','rPL','rOFC'},...
    {'d','t','a','b','lg','hg'});
%%
feat = nameVect(pAdj<0.05)';
sigFeat = logicFind(0.05,pAdj,'<');
dir = zeros(size(sigFeat,2),1);
for ii = 1:size(sigFeat,2)
   dir(ii,1) = mean(nvhl(:,sigFeat(ii)))>mean(sham(:,sigFeat(ii))); 
end
%%
for ii = 1:size(sigFeat,2)
    figure
    hold on
    plot(zeros(size(sham,1),1),sham(:,sigFeat(ii)),'.b','MarkerSize',15)
    plot(ones(size(nvhl,1),1),nvhl(:,sigFeat(ii)),'.r','MarkerSize',15)
    xlim([-1 2])
    set(gca,'XTickLabel',{'Sham','NVHL'},'XTick',[0,1])
    title(nameVect(sigFeat(ii)))
end
%%
figure
for ii = 1:28
plot(coh{1}.f,mean(coh{1}.Cxy(ii,:,:),3)); ylim([0 1])
pause(3)
end