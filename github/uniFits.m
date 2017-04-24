function [out] = uniFits(data,resp,reg)
%% Goes through all independent variables (column) and fits regression
% Pull out independent variables, x; assumes dependent variable exists in
% column one

% INPUTS 
% data = data matrix; format: subjects (rows) x observations (columns) with
%   response variable at beginning
% resp = last column index of response variable; format: integer
% reg = type of regression; options: 'binomial','linear'
%%
% Pull out independent variables, x
x = data(:,resp+1:end);
% Pull out dependent variables, y
y = data(:,1:resp);
% Cycle through all x's and perform regression
for xi = 1:size(x,2)
    [~,~,stats] = glmfit(x(:,xi),y,reg);
    % Store beta coefficient; either OR in case of binomial or raw beta for
    % linear
    switch reg
        case 'binomial'
            out(1,xi) = exp(stats.beta(2));
        case 'linear'
            out(1,xi) = stats.b(2);
    end
    % Store p-value
    out(2,xi) = stats.p(2);
end