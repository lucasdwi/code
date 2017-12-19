% [dataRel,~,~,trlDataRel] = collateData('C:\Users\Pythia\Documents\GreenLab\data\angela\New folder\',{'.mat','in'},{'pow','coh'},'trl','rel');
[dataAbsAll,~,files,trlDataAbsAll] = collateData('C:\Users\Pythia\Documents\GreenLab\data\angela\all\',{'.mat','in'},{'pow','coh'},'trl','');
[dataAbsRest,~,filesRest,trlDataAbsRest] = collateData('C:\Users\Pythia\Documents\GreenLab\data\angela\rest\',{'rest','in'},{'pow','coh'},'trl','');
%%
load('drinkResp.mat')
load('bingeResp.mat')
%% Standardize pooled data
% xAll = zscore(trlDataAbsAll);
% xRest = zscore(trlDataAbsRest);
xAll = zscore(trlDataRelAll);
xRest = zscore(trlDataRelRest);
%% Do univariate logistics
oneInds = logicFind(1,nacc,'==');
zeroInds = logicFind(0,nacc,'==');
resp = [];
for ii = 1:20
   if ismember(ii,oneInds)
      resp = [resp;ones(size(dataAbsAll{1,1}{ii,1},1),1)];
   else
       resp = [resp;zeros(size(dataAbsAll{1,1}{ii,1},1),1)]; 
   end
end
%%
for ii = 1:20
    trainInd = ismember(1:size(resp,1),randperm(size(resp,1),round(size(resp,1)*.8)));
    testInd = ~trainInd;
    for jj = 1:60
        mdl = fitglm(xAll(trainInd,jj),resp(trainInd),'distribution','binomial');
        yhat = predict(mdl,xAll(testInd,jj));
        [~,~,~,a(ii,jj)] = perfcurve(resp(testInd),yhat,1);
    end
end
%%
oneInds = logicFind(1,nacc,'==');
zeroInds = logicFind(0,nacc,'==');
resp = [];
for ii = 1:20
   if ismember(ii,oneInds)
      resp = [resp;ones(size(dataAbsRest{1,1}{ii,1},1),1)];
   else
       resp = [resp;zeros(size(dataAbsRest{1,1}{ii,1},1),1)]; 
   end
end
%%
for ii = 1:20
    trainInd = ismember(1:size(resp,1),randperm(size(resp,1),round(size(resp,1)*.8)));
    testInd = ~trainInd;
    for jj = 1:60
        mdl = fitglm(xRest(trainInd,jj),resp(trainInd),'distribution','binomial');
        yhat = predict(mdl,xRest(testInd,jj));
        [~,~,~,aRest(ii,jj)] = perfcurve(resp(testInd),yhat,1);
    end
end
%% Get betas
load('bingeRespAllOutcomeRel.mat')
bingeBetaAll = allData.allBeta{1,1}.survBeta;
bingeShellBetas = logicFind(0,bingeBetaAll(1,:),'~=');
bingeCoreBetas = logicFind(0,bingeBetaAll(2,:),'~=');
load('NAccDecreaseAll.mat')
NAccBetas = logicFind(0,allData.allBeta{1,1}.survBeta(1,:),'~=');
load('PFCDecreaseAll.mat')
PFCBetas = logicFind(0,allData.allBeta{1,1}.survBeta(1,:),'~=');
%% Run PCA on xAll
[coeff,score,latent,tsquared,explained,mu] = pca(xAll,'Centered','off');
%% Plot first three PCs
figure
scatter3(score(:,1),score(:,2),score(:,3),'.b');
% Apply loadings (coeff) to new data (rest)
restScore = xRest*coeff';
hold on
scatter3(restScore(:,1),restScore(:,2),restScore(:,3),'.r')
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
legend({'All','Rest'})
% [xiny,yinx,percOver] = nOverlap(score(:,1:3),restScore(:,1:3));
% plot(NaN,NaN,'Marker','none','LineStyle','none')
% legend({['All: ',num2str(round(xiny,2)*100),'%'],['Rest: ',num2str(round(yinx,2)*100),'%'],['Overlap: ',num2str(round(percOver,2)*100),'%']})
title('Drink Rel: All vs. Rest')
%% Get centroids
allCentroid = mean(score(:,1:3),1);
restCentroid = mean(restScore(:,1:3),1);
%% Split into NA shell responders vs. non-responders - All Binge
oneInds = logicFind(1,shell,'==');
zeroInds = logicFind(0,shell,'==');
resp = [];
for ii = 1:24
   if ismember(ii,oneInds)
      resp = [resp;ones(size(allData{ii,1},1),1)];
   else
       resp = [resp;zeros(size(allData{ii,1},1),1)]; 
   end
end

figure
hold on
scatter3(score(resp==1,1),score(resp==1,2),score(resp==1,3),'.g');
scatter3(score(resp==0,1),score(resp==0,2),score(resp==0,3),'.c');
% [xiny,yinx,percOver] = nOverlap(score(resp==1,1:3),score(resp==0,1:3));
% plot(NaN,NaN,'Marker','none','LineStyle','none')
title('Binge All: Shell Response vs. Non-Response')
legend({'Responder','Non-Responder'})
% legend({['Response: ',num2str(round(xiny,2)*100),'%'],['Non-Response: ',num2str(round(yinx,2)*100),'%'],['Overlap: ',num2str(round(percOver,2)*100),'%']})
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
respCentroid = mean(score(resp==1,1:3),1);
nonRespCentroid = mean(score(resp==0,1:3),1);
%% Split into NA shell responders vs. non-responders - Rest Binge
oneInds = logicFind(1,shell,'==');
zeroInds = logicFind(0,shell,'==');
resp = [];
for ii = 1:24
   if ismember(ii,oneInds)
      resp = [resp;ones(size(restData{ii,1},1),1)];
   else
       resp = [resp;zeros(size(restData{ii,1},1),1)]; 
   end
end

figure
hold on
scatter3(restScore(resp==1,1),restScore(resp==1,2),restScore(resp==1,3),'.g');
scatter3(restScore(resp==0,1),restScore(resp==0,2),restScore(resp==0,3),'.c');
% [xiny,yinx,percOver] = nOverlap(restScore(resp==1,1:3),restScore(resp==0,1:3));
% plot(NaN,NaN,'Marker','none','LineStyle','none')
title('Binge Rest: Shell Response vs. Non-Response')
legend({'Responder','Non-Responder'})
% legend({['Response: ',num2str(round(xiny,2)*100),'%'],['Non-Response: ',num2str(round(yinx,2)*100),'%'],['Overlap: ',num2str(round(percOver,2)*100),'%']})
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
%% Split into NA core responders vs. non-responders - All Binge
oneInds = logicFind(1,core,'==');
zeroInds = logicFind(0,core,'==');
resp = [];
for ii = 1:24
   if ismember(ii,oneInds)
      resp = [resp;ones(size(allData{ii,1},1),1)];
   else
       resp = [resp;zeros(size(allData{ii,1},1),1)]; 
   end
end

figure
hold on
scatter3(score(resp==1,1),score(resp==1,2),score(resp==1,3),'.g');
scatter3(score(resp==0,1),score(resp==0,2),score(resp==0,3),'.c');
% [xiny,yinx,percOver] = nOverlap(score(resp==1,1:3),score(resp==0,1:3));
% plot(NaN,NaN,'Marker','none','LineStyle','none')
title('Binge All: Core Response vs. Non-Response')
legend({'Responder','Non-Responder'})
% legend({['Response: ',num2str(round(xiny,2)*100),'%'],['Non-Response: ',num2str(round(yinx,2)*100),'%'],['Overlap: ',num2str(round(percOver,2)*100),'%']})
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
%% Split into NA core responders vs. non-responders - Rest Binge
oneInds = logicFind(1,core,'==');
zeroInds = logicFind(0,core,'==');
resp = [];
for ii = 1:24
   if ismember(ii,oneInds)
      resp = [resp;ones(size(restData{ii,1},1),1)];
   else
       resp = [resp;zeros(size(restData{ii,1},1),1)]; 
   end
end

figure
hold on
scatter3(restScore(resp==1,1),restScore(resp==1,2),restScore(resp==1,3),'.g');
scatter3(restScore(resp==0,1),restScore(resp==0,2),restScore(resp==0,3),'.c');
% [xiny,yinx,percOver] = nOverlap(restScore(resp==1,1:3),restScore(resp==0,1:3));
% plot(NaN,NaN,'Marker','none','LineStyle','none')
title('Binge Rest: Core Response vs. Non-Response')
legend({'Responder','Non-Responder'})
% legend({['Response: ',num2str(round(xiny,2)*100),'%'],['Non-Response: ',num2str(round(yinx,2)*100),'%'],['Overlap: ',num2str(round(percOver,2)*100),'%']})
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
%% Split into any responders vs. non-responders - All Binge
oneInds = [1:8,11:14,17:24];
zeroInds = [9,10,15,16];
resp = [];
for ii = 1:24
   if ismember(ii,oneInds)
      resp = [resp;ones(size(allData{ii,1},1),1)];
   else
       resp = [resp;zeros(size(allData{ii,1},1),1)]; 
   end
end

figure
hold on
scatter3(score(resp==1,1),score(resp==1,2),score(resp==1,3),'.g');
scatter3(score(resp==0,1),score(resp==0,2),score(resp==0,3),'.c');
title('Binge Rest: Any Response vs. Non-Response')
legend({'Responder','Non-Responder'})
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
%% Split into NA shell responders vs. non-responders - All Drink
oneInds = logicFind(1,nacc,'==');
zeroInds = logicFind(0,nacc,'==');
resp = [];
for ii = 1:20
   if ismember(ii,oneInds)
      resp = [resp;ones(size(dataAbsAll{1,1}{ii,1},1),1)];
   else
       resp = [resp;zeros(size(dataAbsAll{1,1}{ii,1},1),1)]; 
   end
end

figure
hold on
scatter3(score(resp==1,1),score(resp==1,2),score(resp==1,3),'.g');
scatter3(score(resp==0,1),score(resp==0,2),score(resp==0,3),'.c');
% [xiny,yinx,percOver] = nOverlap(score(resp==1,1:3),score(resp==0,1:3));
% plot(NaN,NaN,'Marker','none','LineStyle','none')
title('Drink All: NAcc Response vs. Non-Response')
legend({'Responder','Non-Responder'})
% legend({['Response: ',num2str(round(xiny,2)*100),'%'],['Non-Response: ',num2str(round(yinx,2)*100),'%'],['Overlap: ',num2str(round(percOver,2)*100),'%']})
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
%% Split into NA shell responders vs. non-responders - Rest Drink
oneInds = logicFind(1,nacc,'==');
zeroInds = logicFind(0,nacc,'==');
resp = [];
for ii = 1:20
   if ismember(ii,oneInds)
      resp = [resp;ones(size(dataAbsRest{1,1}{ii,1},1),1)];
   else
       resp = [resp;zeros(size(dataAbsRest{1,1}{ii,1},1),1)]; 
   end
end

figure
hold on
scatter3(restScore(resp==1,1),restScore(resp==1,2),restScore(resp==1,3),'.g');
scatter3(restScore(resp==0,1),restScore(resp==0,2),restScore(resp==0,3),'.c');
[xiny,yinx,percOver] = nOverlap(restScore(resp==1,1:3),restScore(resp==0,1:3));
plot(NaN,NaN,'Marker','none','LineStyle','none')
title('Drink Rest: NAcc Response vs. Non-Response')
legend({['Response: ',num2str(round(xiny,2)*100),'%'],['Non-Response: ',num2str(round(yinx,2)*100),'%'],['Overlap: ',num2str(round(percOver,2)*100),'%']})
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
%% Split into PFC responders vs. non-responders - All Drink
oneInds = logicFind(1,pfc,'==');
zeroInds = logicFind(0,pfc,'==');
resp = [];
for ii = [1:14,17:20]
   if ismember(ii,oneInds)
      resp = [resp;ones(size(dataAbsAll{1,1}{ii,1},1),1)];
   else
       resp = [resp;zeros(size(dataAbsAll{1,1}{ii,1},1),1)]; 
   end
end
% Remove missing PFC data
pfcScore = score([1:8513,9700:end],:);
figure
hold on
scatter3(pfcScore(resp==1,1),pfcScore(resp==1,2),pfcScore(resp==1,3),'.g');
scatter3(pfcScore(resp==0,1),pfcScore(resp==0,2),pfcScore(resp==0,3),'.c');
[xiny,yinx,percOver] = nOverlap(score(resp==1,1:3),score(resp==0,1:3));
plot(NaN,NaN,'Marker','none','LineStyle','none')
title('Drink All: PFC Response vs. Non-Response')
legend({['Response: ',num2str(round(xiny,2)*100),'%'],['Non-Response: ',num2str(round(yinx,2)*100),'%'],['Overlap: ',num2str(round(percOver,2)*100),'%']})
xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
%%
[allK,allV] = convhull(pfcScore(resp==1,1),pfcScore(resp==1,2),pfcScore(resp==1,3),'simplify',true);
[restK,restV] = convhull(pfcScore(resp==0,1),pfcScore(resp==0,2),pfcScore(resp==0,3),'simplify',true);
figure
trimesh(allK,pfcScore(resp==1,1),pfcScore(resp==1,2),pfcScore(resp==1,3),'FaceAlpha',0,'EdgeColor','b')
hold on
trimesh(restK,pfcScore(resp==0,1),pfcScore(resp==0,2),pfcScore(resp==0,3),'FaceAlpha',0,'EdgeColor','r')
%% Split into PFC responders vs. non-responders - Rest Drink
oneInds = logicFind(1,nacc,'==');
zeroInds = logicFind(0,nacc,'==');
resp = [];
for ii = [1:14,17:20]
   if ismember(ii,oneInds)
      resp = [resp;ones(size(dataAbsRest{1,1}{ii,1},1),1)];
   else
       resp = [resp;zeros(size(dataAbsRest{1,1}{ii,1},1),1)]; 
   end
end
% Remove missing PFC data
pfcRestScore = restScore([1:988,1208:end],:);
figure
hold on
scatter3(pfcRestScore(resp==1,1),pfcRestScore(resp==1,2),pfcRestScore(resp==1,3),'.g');
scatter3(pfcRestScore(resp==0,1),pfcRestScore(resp==0,2),pfcRestScore(resp==0,3),'.c');
[xiny,yinx,percOver] = nOverlap(restScore(resp==1,1:3),restScore(resp==0,1:3));
plot(NaN,NaN,'Marker','none','LineStyle','none')
title('Drink Rest: PFC Response vs. Non-Response')
legend({['Response: ',num2str(round(xiny,2)*100),'%'],['Non-Response: ',num2str(round(yinx,2)*100),'%'],['Overlap: ',num2str(round(percOver,2)*100),'%']})
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
%%
[allK,allV] = convhull(pfcRestScore(resp==1,1),pfcRestScore(resp==1,2),pfcRestScore(resp==1,3),'simplify',true);
[restK,restV] = convhull(pfcRestScore(resp==0,1),pfcRestScore(resp==0,2),pfcRestScore(resp==0,3),'simplify',true);
figure
trimesh(allK,pfcRestScore(resp==1,1),pfcRestScore(resp==1,2),pfcRestScore(resp==1,3),'FaceAlpha',0,'EdgeColor','b')
hold on
trimesh(restK,pfcRestScore(resp==0,1),pfcRestScore(resp==0,2),pfcRestScore(resp==0,3),'FaceAlpha',0,'EdgeColor','r')
%% Calculate convex hull and volume
[allK,allV] = convhull(score(:,1),score(:,2),score(:,3),'simplify',true);
[restK,restV] = convhull(restScore(:,1),restScore(:,2),restScore(:,3),'simplify',true);
figure
trimesh(allK,score(:,1),score(:,2),score(:,3),'FaceAlpha',0,'EdgeColor','b')
hold on
trimesh(restK,restScore(:,1),restScore(:,2),restScore(:,3),'FaceColor','r','FaceAlpha',1,'EdgeColor','k')
%%
figure
hold on
trimesh(allK,score(:,1),score(:,2),score(:,3),'FaceAlpha',0,'EdgeColor','b')
scatter3(score(:,1),score(:,2),score(:,3),'.b');
trimesh(restK,restScore(:,1),restScore(:,2),restScore(:,3),'FaceAlpha',0,'EdgeColor','r')
scatter3(restScore(:,1),restScore(:,2),restScore(:,3),'.r')
%% Calculate centroids
x1 = mean(score(:,1)); y1 = mean(score(:,2)); z1 = mean(score(:,3));
x2 = mean(restScore(:,1)); y2 = mean(restScore(:,2)); z2 = mean(restScore(:,3));
scatter3(x1,y1,z1,80,'.g')
scatter3(x2,y2,z2,80,'.y')
%%
nacc = [1,1,1,0,0,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0];
%%
one = [1:3,6:10,13,14];
oneAll = cat(1,all{one});%*coeff;
zero = [4,5,11,12,15:20];
zeroAll = cat(1,all{zero});%*coeff;
%%
figure
hold on
scatter3(oneAll(:,1),oneAll(:,2),oneAll(:,3),'.b')
scatter3(zeroAll(:,1),zeroAll(:,2),zeroAll(:,3),'.r')
%%
figure
hold on
scatter3(score(5512:6112,1),score(5512:6112,2),score(5512:6112,3),'.b')
scatter3(restScore(425:585,1),restScore(425:585,2),restScore(425:585,3),'.r')
%%
figure
hold on
scatter3(score(5512:6112,1),score(5512:6112,2),score(5512:6112,3),'.b')
[n,c] = hist3([score(:,1),score(:,2)]);
[~,h] = contour(c{2},c{1},n);
h.ContourZLevel = -10;
[n,c] = hist3([score(:,2),score(:,3)]);
[~,h] = contour(c{2},c{1},n);
h.ContourXLevel = -10;