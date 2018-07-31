%% Grab data
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'male';'di'},{'pow','coh'},'avg','');
male = cat(1,data{1,1}{:});
di = cat(1,data{1,2}{:});
%% Get mean and standard deviation of data
mMale = mean(male,1);
sMale = std(male,[],1);
mDi = mean(di,1);
sDi = std(di,[],1);
%% Reorder for simplicity - channel pair X band
mMale = reshape(mMale(:,25:end),6,6)';
sMale = reshape(sMale(:,25:end),6,6)';
mDi = reshape(mDi(:,25:end),6,6)';
sDi = reshape(sDi(:,25:end),6,6)';
%% Plot coherence; each channel pair gets a panel
titles = {'PL-PR','PL-SL','PL-SR','PR-SL','PR-SR','SL-SR'};
figure
for ii = 1:6
    subplot(2,3,ii)
    barwitherr([sMale(ii,:);sDi(ii,:)]'.*100,[mMale(ii,:);mDi(ii,:)]'.*100)
    title(titles{ii})
    set(gca,'XTickLabel',{'\Delta','\theta','\alpha','\beta','l\gamma','h\gamma'});
    xlim([0.5 6.5])
    box off
end
%%
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\','male','ex','Plots','ex');
mFemale = [];
sFemale = [];
for ii = 1:size(files,2)
    load(files{ii},'coh')
    mFemale = cat(1,mFemale,coh{1}.mCxy(6,:));
%     sFemale = cat(1,sFemale,coh{1}.sdCxy(1,:));
end
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\','male','in');
mMale = [];
sMale = [];
for ii = 1:size(files,2)
    load(files{ii},'coh')
    mMale = cat(1,mMale,coh{1}.mCxy(6,:));
%     sMale = cat(1,sMale,coh{1}.sdCxy(1,:));
end
figure
hold on
shadedErrorBar(coh{1}.freq,mean(mFemale,1),std(mFemale,[],1),'-r',1)
shadedErrorBar(coh{1}.freq,mean(mMale,1),std(mMale,[],1),'-b',1)
xlabel('Frequency (Hz)')
ylabel('Coherence')
title('lNac-rNAc')