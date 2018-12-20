%%
[data,samp,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\'...
    'data\angelaMIA\processed\'],{'mSaline';'mMIA';'fSaline';'fMIA'},...
    {'pow','coh'},'avg','rel');
mSal = cat(1,data{1,1}{:});
mP = cat(1,data{1,2}{:});
fSal = cat(1,data{1,3}{:});
fP = cat(1,data{1,4}{:});
save(['C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\analyzed\'...
    '2cohortAvgData.mat'],'data','files','samp','mSal','mP','fSal','fP')
%%
for ii = 1:100
    load(['C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\analyzed\'...
        'sexModels\',num2str(ii),'.mat'])
    fA(ii) = female.a;
    fAcc(ii,:) = 1-female.err;
    fARand(ii) = female.aRand;
    fAccRand(ii,:) = 1-female.errRand;
    fX(ii,:) = female.x;
    fY(ii,:) = female.y;
    fB(ii,:) = female.beta{1}.survBeta;
    fBS(ii,:) = female.beta{1}.signBeta;
    
    mA(ii) = male.a;
    mAcc(ii,:) = 1-male.err;
    mARand(ii) = male.aRand;
    mAccRand(ii,:) = 1-male.errRand;
    mX(ii,:) = male.x;
    mY(ii,:) = male.y;
    mB(ii,:) = male.beta{1}.survBeta;
    mBS(ii,:) = male.beta{1}.signBeta;
end
figure
hold on
plot(mean(mX,1),mean(mY,1))
plot(mean(fX,1),mean(fY,1))
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
box off
xlabel('FPR'); ylabel('TPR')
legend({['Male: ',num2str(round(mean(mA),2)),'\pm',...
    num2str(round(std(mA),2))];['Female: ',num2str(round(mean(fA),2)),...
    '\pm',num2str(round(std(fA),2))]},'location','se')
title('Saline vs. Poly I:C')

doubleHist(reshape(mAcc,1,10000),reshape(mAccRand,1,10000),...
    'xlab','Accuracy','main','Male');
doubleHist(reshape(fAcc,1,10000),reshape(fAccRand,1,10000),...
    'xlab','Accuracy','main','Female');
%% male female; logistics
for ii = 1:100
    disp(num2str(ii))
    sInds = randperm(size(mSal,1),3);
    pInds = randperm(size(mP,1),3);
    
    mTrainX = [mSal(~ismember(1:size(mSal,1),sInds),:);...
        mP(~ismember(1:size(mP,1),pInds),:)];
    mTrainY = [zeros(size(mSal,1)-3,1);ones(size(mP,1)-3,1)];
    mTestX = [mSal(sInds,:);mP(pInds,:)];
    mTestY = [0;0;0;1;1;1];
    
    sInds = randperm(size(mSal,1),3);
    pInds = randperm(size(mP,1),3);
    
    fTrainX = [fSal(~ismember(1:size(fSal,1),sInds),:);...
        fP(~ismember(1:size(fP,1),pInds),:)];
    fTrainY = [zeros(size(fSal,1)-3,1);ones(size(fP,1)-3,1)];
    fTestX = [fSal(sInds,:);fP(pInds,:)];
    fTestY = [0;0;0;1;1;1];
    for jj = 1:216
        mdl = fitglm(mTrainX(:,jj),mTrainY,'distribution','binomial');
        prob = predict(mdl,mTestX(:,jj));
        [mLogX{jj}(ii,:),mLogY{jj}(ii,:),~,mLogA(ii,jj)] = perfcurve(...
            mTestY,prob,1,'Tvals',0:.1:1,'UseNearest',0);
        
        mdl = fitglm(fTrainX(:,jj),fTrainY,'distribution','binomial');
        prob = predict(mdl,fTestX(:,jj));
        [fLogX{jj}(ii,:),fLogY{jj}(ii,:),~,fLogA(ii,jj)] = perfcurve(...
            fTestY,prob,1,'Tvals',0:.1:1,'UseNearest',0);
    end
end
fLogAM = mean(fLogA,1); 
mLogAM = mean(mLogA,1);
label = {'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'};
nameVect = names(label,{'d','t','a','b','lg','hg'});
[fFeatV,fInd] = sort(fLogAM','descend');
[mFeatV,mInd] = sort(mLogAM','descend');
fFeat = nameVect(fInd)';
mFeat = nameVect(mInd)';
mSurv = mean(mB,1);
fSurv = mean(fB,1);
mSurvS = mSurv(mInd)';
fSurvS = fSurv(fInd)';
% Get proxy sign of models by comparing means
for ii = 1:216
    mSignLog(ii,1) = mean(mP(:,ii))>mean(mSal(:,ii));
    fSignLog(ii,1) = mean(fP(:,ii))>mean(fSal(:,ii));
end
mSignLog = mSignLog(mInd);
fSignLog = fSignLog(fInd);
% Find tiers for log AUCs
mTier = tier(mLogA(:,mInd));
fTier = tier(fLogA(:,fInd));
%%
load('C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\analyzed\avgData.mat')
nameVect = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},{'d','t','a','b','lg','hg'});
%%
inds = 49:6:216;
s = data(1:15,:);
poly = data(16:36,:);
%%
for ii = 1:size(inds,2)
    [~,p(ii)] = ttest2(s(:,inds(ii)),poly(:,inds(ii)));
    if p(ii)*28 < 0.05
       figure
       plot(ones(1,size(poly,1)),poly(:,inds(ii)),'ob')
       hold on
       plot([0.75 1.25 ],repmat(mean(poly(:,inds(ii))),1,2),'-k')
       
       plot(zeros(1,size(s,1)),s(:,inds(ii)),'or')
       plot([-.25 .25],repmat(mean(s(:,inds(ii))),1,2),'-k')
       
       xlim([-0.5 1.5])
       title(nameVect{inds(ii)})
       set(gca,'XTick',[0 1],'XTickLabel',{'S','Poly'})
       box off
    end
end

