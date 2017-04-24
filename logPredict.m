function [AUC,xAx,yAx,dev,tpr,fpr] = logPredict(x,y)
%% Uses data in x to predict y through logisitic regression and ROC
% INPUTS
% x = data matrix; format = observations x variables
% y = response matrix, may be multinomial; format = observations x categories

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
   % Expand betas to matrix
   betaMat = repmat(beta{yi}',nObv,1);
   % Calulate yhat from betas and data
   yhat = sum([ones(nObv,1),x].*betaMat,2);
   % Convert yhat to logisitic probability
   prob = exp(yhat)./(1+exp(yhat));
   % Compare probability to actual with varying thresholds
   [xAx{yi},yAx{yi},~,AUC{yi}] = perfcurve(y(:,yi),prob,1);
   [tpr,fpr] = rocCurve(y,prob,0:0.1:1);
   % Plot ROCs
%    subplot(2,2,yi)
%    hold on
%    plot(xAx{yi},yAx{yi})
end
