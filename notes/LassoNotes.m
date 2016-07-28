% index for each predictor columns of T
% mnrfit(X,Y) where X = predictors (variables), Y = response (group)
cond = fieldnames(T);
% Run Lasso for each response
alph = [.1:.1:1];
% for c = 1:length(cond)
    
predict = table2array(T.Base(:,8:end))+1;
% Run cross validation 100 times per alpha 
%for k = 1:100
    for j = 1:length(alph)  
        for ii = 2:5
            response = table2array(T.Base(:,ii));
            % Store results of Lasso in structures at ii-1 (1:4)
            [B{j,ii-1,k},FitInfo{j,ii-1,k}] = lassoglm(predict,response,'binomial','CV',10,'Alpha',alph(j));
            % Store min SE and Dev
            SE(j,ii-1,k) = FitInfo{j,ii-1,k}.SE(FitInfo{j,ii-1,k}.Index1SE);
            lambda(j,ii-1,k) = FitInfo{j,ii-1,k}.Lambda1SE;
            lambdaInd(j,ii-1,k) = length(FitInfo{j,ii-1,k}.Intercept) - FitInfo{j,ii-1,k}.Index1SE;
        end
    end
    disp(k)
%end
% end
%% Provides alpha level with least mean SE
for ii = 1:size(SE,1)
    SE(ii,5) = mean(SE(ii,1:4));
end
[min I] = min(SE(:,5));
minalpha = alph(I);
%% Run lassoglm at different numbers of lambda
clear B FitInfo;
%lams = [25];
for ii = 2:5
    response = table2array(T.Base(:,ii));
    [B{ii-1},FitInfo{ii-1}] = lassoglm(predict,response,'binomial','CV',10);
end
    
%% Plot Lasso Traces
for ii = 1:4
    lassoPlot(B{1,ii},FitInfo{1,ii},'PlotType','CV')
    lassoPlot(B{1,ii},FitInfo{1,ii},'PlotType','Lambda','XScale','log')
end

%% Collect B results from Lassos and find all potential variables and
% variables for best lambda for each Lasso

% First extract non-zero rows from each B and cross-reference these with
% RowNames of table

for ii = 1:length(B);
    % Find all variables with any non-zero coefficients across lambda
    nonz{1}(:,ii) = sum(abs(B{ii}),2);
    % Get best fit lambda index
    lambda(ii) = FitInfo{ii}.Index1SE;
    % Find best fit variables with non-zero coefficients
    nonz{2}(:,ii) = B{ii}(:,lambda(ii));
end
%%
predict = table2array(T.Base(:,8:end))+1;
response = table2array(T.Base(:,2));
a = [.1:.1:1];
%%
lambda = []; mse = []; alph = []; minLambda = []; B = {}; FitInfo= {};
% figure;
% xlabel = {'Lambda'}; ylabel = {'MSE'};
%for j = 1:10;
    for ii = 1:length(a)
        [B{j,ii} FitInfo{j,ii}] = lassoglm(predict,response,'binomial','CV',10);
        lambda{j,ii} = FitInfo{j,ii}.Lambda;
        mse{j,ii} = FitInfo{j,ii}.MSE;
    %minLambda(ii) = FitInfo{ii}.Lambda1SE;
%     thislambda = FitInfo.Lambda;
%     lambda = horzcat(lambda,thislambda);
%     thismse = FitInfo.MSE;
%     mse = horzcat(mse,thismse);
%     thisalph(1:length(thislambda)) = a(ii);
%     alph = horzcat(alph,thisalph);
%     thisalph = [];
%     hold on
%     plot(lambda{1,ii},mse{1,ii})
    end
%end
%%
figure; hold on;
for ii = 1:10
    plot(lambda{ii,1},mse{ii,1})
end
