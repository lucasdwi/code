load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\paper1data.mat')
load('C:\Users\Lucas\Desktop\GreenLab\data\paper1\finalData\indGroups.mat')
%%
y = y(:,[1,3]);
normX = zscore(x);
%%
for ii = 1:size(inds,1)
%     for ri = 1:100
%         for ni = 1:size(normX,1)
%            thisPerm(ni,:) = randperm(size(normX,2));
%         end
% Randomized
%     [coeff,scores,latent,~,explained] = pca(normX(thisPerm(:,inds{ii})),'Economy',false);
% Go through different inds
    [coeff,scores,latent,~,explained] = pca(normX(:,inds{ii}),'Economy',false);
% Use all x
%     [coeff,scores,latent,~,explained] = pca(normX,'Economy',false);
    % Plot number of prinicple components vs. % of variance explained
    percVar = cumsum(explained);
%     figure;
%     stairs(percVar)
%     title('Variance Explained by {\itn} Princinple Components')
%     xlabel('{\itn} Principal Components'); ylabel('% of Variance Explained')
    % Use first n pcs for regression
    ind(ii) = logicFind(80,percVar,'>=','first');
    pcVar(ii) = percVar(ind(ii));
    [AUC(ii,:),xMat{ii},yMat{ii},dev{ii}] = logPredict(scores(:,1:ind(ii)),y);
%     [AUC(ii,:,ri),xMat{ii,ri},yMat{ii,ri},dev{ii,ri}] = logPredict(scores(:,1:ind(ii)),y);
%     for ri = 1:100
%         thisPerm = randperm(size(y,1));
%         yRand = y(thisPerm',:);
%         [randAUC(ii,:,ri),xMatRand{ii,ri},yMatRand{ii,ri},devRand{ii,ri}] = logPredict(scores(:,1:ind(ii)),yRand);
%     end
%     end
end
%% Get random AUCs and ROCs
aucRandMean = mean(randAUC,3);
aucRandSd = std(randAUC,[],3);
for xi = 1:size(xMatRand,1)
    for yi = 1:size(xMatRand,2)
        for zi = 1:2
            if size(xMatRand{xi,yi}{1,zi},1) < 25
                catRandX{zi}(xi,yi,:) = cell2mat([xMatRand{xi,yi}(1,zi);ones(25-size(xMatRand{xi,yi}{1,zi},1),1)]);
                catRandY{zi}(xi,yi,:) = cell2mat([yMatRand{xi,yi}(1,zi);ones(25-size(yMatRand{xi,yi}{1,zi},1),1)]);
            else
                catRandX{zi}(xi,yi,:) = cell2mat(xMatRand{xi,yi}(1,zi));
                catRandY{zi}(xi,yi,:) = cell2mat(yMatRand{xi,yi}(1,zi));
            end
        end
    end
end
xAvgRand = sq(mean(catRandX{2},2));
yAvgRand = sq(mean(catRandY{2},2));
%% Set up real AUC matrix
% aucComp = [AUC(1,:);AUC(2,:);AUC(3,:)];
% devComp = [dev{1,1};dev{1,2};dev{1,3}];
%% Create real AUC table
aucTab = table(AUC(:,1),AUC(:,2),AUC(:,3),AUC(:,4),'RowNames',indName,'VariableNames',{'SA','SS','CA','CS'});