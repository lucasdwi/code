function [AUC,xAx,yAx,dev,tpr,fpr,prob] = logPredict(x,y,testX,testY)
%% Uses data in x to predict y through logisitic regression and ROC
% INPUTS
% x = data matrix; format = observations x variables
% y = response matrix, may be multinomial; format = observations x categories
% testX = optional input if using test data to calculate AUC; format = same
%   as x
% testY = optional input if using test data to calculate AUC; format = same
%   as y, but dimensions should match testX

% OUTPUTS
% AUC = area under the reciever operator curve
% xAx = cell array of x data for plotting
% yAx = cell arary of y data for plotting
%%
% Get number of observations
nObv = size(x,1);
% Cycle through all responses - y columns
% figure 
for yi = 1:size(y,2)
    % Get beta coefficients from logistic regression
   [beta{yi},dev{yi}] = glmfit(x,y(:,yi),'binomial','link','logit');
   % Calulate yhat from betas and data; if using test set, use testX
   if ~isempty(testX)
       % Expand betas to matrix
       betaMat = repmat(beta{yi}',size(testX,1),1);
       yhat = sum([ones(size(testX,1),1),testX].*betaMat,2);
       % Convert yhat to logisitic probability
       prob{yi} = exp(yhat)./(1+exp(yhat));
   % Otherwise use x
   else
       % Expand betas to matrix
       betaMat = repmat(beta{yi}',nObv,1);
       yhat = sum([ones(nObv,1),x].*betaMat,2);
       % Convert yhat to logisitic probability
       prob(:,yi) = exp(yhat)./(1+exp(yhat));
   end
   % Compare probability to actual with varying thresholds; if using testY
   if ~isempty(testY)
       %        [xAx{yi},yAx{yi},~,AUC{yi}] = perfcurve(testY(:,yi),prob,1);
       %        [tpr,fpr] = rocCurve(testY,prob,0:0.1:1);
       AUC = []; xAx = []; yAx = []; dev= []; tpr = []; fpr = [];
   else
       [xAx{yi},yAx{yi},~,AUC{yi}] = perfcurve(y(:,yi),prob(:,yi),1);
       %        [tpr,fpr] = rocCurve(y,prob(:,yi),0:0.1:1);
       tpr = []; fpr = [];
   end
   % Plot ROCs
%    subplot(2,2,yi)
%    hold on
%    plot(xAx{yi},yAx{yi})
end
