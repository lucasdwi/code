%% Calculate absolute difference
for i = 1:length(files.Base)
    BaseMeanPSDs{i} = (PSDs.Base.Bands{i,1}+PSDs.Base.Bands{i,2})./2;
    BaseDiffPSDs{i} = (PSDs.Base.Bands{i,1}-PSDs.Base.Bands{i,2});
end
%%
normMat = [1/3,1/3,1/3,1/3;1/5,1/5,1/5,1/5;1/15,1/15,1/15,1/15;1/18,1/18,1/18,1/18;1/30,1/30,1/30,1/30];
for i = 1:length(files.Base)
    normMean{i} = BaseMeanPSDs{i}.*normMat;
    normDiff{i} = BaseDiffPSDs{i}.*normMat;
end
%% Log transformed
catlogx = []; catlogy = [];
for i = 1:length(files.Base)
    LogBaseMeanPSDs{i} = (log2(PSDs.Base.Bands{i,1}.*normMat)+log2(PSDs.Base.Bands{i,2}.*normMat))./2;
    LogBaseDiffPSDs{i} = log2(PSDs.Base.Bands{i,1}.*normMat)-log2(PSDs.Base.Bands{i,2}.*normMat);
    logx = reshape(LogBaseMeanPSDs{i},[20,1]);
    logy = reshape(LogBaseDiffPSDs{i},[20,1]);
    catlogx = vertcat(catlogx,logx); catlogy = vertcat(catlogy,logy);
end
%%
catx = []; caty = [];
for i = 1:length(files.Base)
    x = reshape(normMean{i},[20,1]);
    y = reshape(normDiff{i},[20,1]);
    catx = vertcat(catx,x); caty = vertcat(caty,y);
end
%%
[h,p] = vartest2(catx,caty);
%%
[h,p] = kstest(caty);
figure; histfit(caty);
[h,p] = kstest(catx);
figure; histfit(catx);
%%
figure;
for i = 1:length(files.Base)
    plot(BaseMeanPSDs{1,i},BaseDiffPSDs{1,i},'.');
    hold on;
end
%%
figure;
for i = 1:length(files.Base)
    plot(normMean{1,i},normDiff{1,i},'.');
    hold on;
end
%%
lmod = fitlm(catx,caty);
figure;
plot(catx,caty,'.');
hold on;
% plot mean diff
refline(0,mean(caty));
% plot 2 std
refline(0,mean(caty)+std(caty));
refline(0,mean(caty)-std(caty));
plot(lmod);
    