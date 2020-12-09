load('IRDM2_singleSite_mPFC_2hr_2019-02-19_PostLSD.mat')
%%
cmbs = nchoosek(1:8,2);
for jj = 1%:28
    figure
    hold on
    plot(coh{1}.f,mean(coh{1}.Cxy(jj,:,:),3))
    plot(f(1:500),squeeze(mean(iCoh(jj,1:500,:),3)))
    for ii = 1:size(trls{1}.trial,3)
        x = trls{1}.trial(cmbs(jj,1),:,ii);
        y = trls{1}.trial(cmbs(jj,2),:,ii);
        [f,iCoh2(jj,:,ii),~,~,~] = cmtm(x,y,1/hist.adfreq,8,0,0,0);
    end
    plot(f(1:500),squeeze(mean(iCoh2(jj,1:500,:),3)))
end