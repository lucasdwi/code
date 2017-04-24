%% Grab enet results from master files
allMinMSE = zeros(2,size(masterMinLam,2));
for m = 1:size(masterMinLam,2)
   allMinMSE(1,m) = masterMinLam{1,m}(1,3);
   allMinMSE(2,m) = masterMinLam{1,m}(2,3);
end
mMSE(1) = mean(allMinMSE(1,:));
mMSE(2) = mean(allMinMSE(2,:));
sMSE(1) = std(allMinMSE(1,:));
sMSE(2) = std(allMinMSE(2,:));
%%
moi = find(masterSurvBeta{1,1}(2,:)>0); % Feature columns from cpp
% Fit linear models
for m = 1:length(moi)
   models{1,m} = fitlm(T.Base,['CCPPost~',allPredict{moi(m)}]);
   models{2,m} = coefTest(models{1,m});
   thismean=[];
   for c = 1:2:size(T.Base,1)
       thismean = horzcat(thismean,mean(T.Base.(allPredict{moi(m)})(c:c+1)));
   end
   means{m} = thismean;
   models{3,m} = fitlm(means{m},T.Base.CCPPost(2:2:end));
   models{4,m} = coefTest(models{3,m});
%    models{5,m} = fitlm(table2array(T.Base(1:2:end,moi(m)+2)),table2array(T.Base(1:2:end,2)));
%    models{6,m} = coefTest(models{5,m});
%    models{7,m} = fitlm(table2array(T.Base(2:2:end,moi(m)+2)),table2array(T.Base(2:2:end,2)));
%    models{8,m} = coefTest(models{7,m});
end
%%
moi = find(masterSurvBeta{1,1}(1,:)>0);
% Fit linear models
for m = 1:length(moi)
   models{1,m} = fitlm(T.Base,['BingeReduct~',allPredict{moi(m)}]);
   models{2,m} = coefTest(models{1,m});
   thismean=[];
   for c = 1:2:size(T.Base,1)
       thismean = horzcat(thismean,mean(T.Base.(allPredict{moi(m)})(c:c+1)));
   end
   means{m} = thismean;
   models{3,m} = fitlm(means{m},T.Base.BingeReduct(2:2:end));
   models{4,m} = coefTest(models{3,m});
%    models{5,m} = fitlm(table2array(T.Base(1:2:end,moi(m)+1)),table2array(T.Base(1:2:end,2)));
%    models{6,m} = coefTest(models{5,m});
%    models{7,m} = fitlm(table2array(T.Base(2:2:end,moi(m)+1)),table2array(T.Base(2:2:end,2)));
%    models{8,m} = coefTest(models{7,m});
end

%% Use mean values for plotting
% Grab every other CPP or bingereductvalue
cpp = T.Base.CCPPost(2:2:end);
br = T.Base.BingeReduct(2:2:end);
% Get significant indices
modelInds = find(cell2mat(models(4,:))<=0.05);
sigmoi = moi(find(cell2mat(models(4,:))<=0.05));
%% Plot those interactions
for m = 1:length(sigmoi)
   figure;
  %p = plot(T.Base.(allPredict{sigmoi(m)})(1:2:end).*100,br,'.k');
   p = plot(T.Base.(allPredict{sigmoi(m)})(1:2:end).*100,cpp,'.k');
   set(gca,'box','off')
   hold on;
   text(p.Parent.XTick(end-2),p.Parent.YTick(end-2),['p ~= ',num2str(models{4,modelInds(m)},2)]);
   text(p.Parent.XTick(end-2),p.Parent.YTick(end-3),['R^2 ~= ',num2str(models{3,modelInds(m)}.Rsquared.Ordinary,2)]);
   lsline
   %ylabel('Binge Reduction'); xlabel('Percent of Total Power')
   ylabel('CPP'); xlabel('Percent of Total Power')
   thisTitle = strrep(allPredict{sigmoi(m)},'_',' ');
   %thisTitle = strrep(allPredict{sigmoi(m)},'_',' ');
   title(thisTitle)
end

