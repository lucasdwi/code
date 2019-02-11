function [d] = doubleHist(real,perm,varargin)
%% Plots overlapping histograms. Designed for the outputs of machine 
% learning to compare real and permuted datasets. 
%__________________________________________________________________________
% INPUTS
% real = one dimensional data vector corresponding to real data (correct
%   assignment between predictors and outcomes)
% perm = one dimensional data vector corresponding to permuted data
%   (shuffled assignment of predictors and outcomes)
%__________________________________________________________________________
% PROPERTIES
% unit = string tag corresponding to the units along x-axis; string
% nBin = number of bins to fix across x-axis; integer
% main = title of plot; string
% xLab = label for x-axis; string
% yLab = label for y-axis; string 
%__________________________________________________________________________
% OUTPUTS
% d = Cohen's d, calculated by converting the Mann-Whitney U statistic
%__________________________________________________________________________
% DEPENDENCIES
% distES.mat
%__________________________________________________________________________
% LLD 2017
%% Check that real and perm are 1 dimensional
assert(min(size(real))==1,'Real data has more than one dimension');
assert(min(size(perm))==1,'Perm data has more than one dimension');
%% Input parsing
p = inputParser;
p.CaseSensitive = false;
addRequired(p,'real',@isnumeric);
addRequired(p,'perm',@isnumeric);
addParameter(p,'unit','',@isstr);
addParameter(p,'nBin',40,@isscalar);
addParameter(p,'main','',@isstr);
addParameter(p,'xlab','',@isstr);
addParameter(p,'ylab','Percent of Models (%)',@isstr);
addParameter(p,'loc','Best',@isstr);
parse(p,real,perm,varargin{:});
%%
figure
hold on
% Plot real in black with white outline
h(1) = histogram(real,'Normalization','probability','FaceColor','k',...
    'FaceAlpha',1,'EdgeColor','w');
% Plot perm in white with black outline
h(2) = histogram(perm,'Normalization','probability','FaceColor','w',...
    'FaceAlpha',1);
% Plot dummy for effect size in legend
h(3) = plot(NaN,NaN,'Marker','none','LineStyle','none');
% Grab length of x-axis
xax = get(gca,'xlim');
len = xax(2)-xax(1);
% Split into n bins
bin = len/p.Results.nBin;
% Set binWidths to bin
set(h(1),'BinWidth',bin);
set(h(2),'BinWidth',bin);
% Calculate p value from null distribution
z = (mean(real)-mean(perm))/std(perm);
% pValue = 2*normcdf(-abs(z));
% Calculate effect size (Cohen's d)
% d = distES(real,perm);
% Calculate t with 95% two sided 
t = tinv(0.975,numel(real)-1);
% Calulate 95% confidence interval
cR = t*(std(real)/sqrt(numel(real)));
cP = t*(std(perm)/sqrt(numel(perm)));
% Determine rounding
r1 = 0; test = 0;
while test == 0
    r1 = r1+1;
    test = round(cR,r1);
end
r2 = 0; test = 0;
while test == 0
    r2 = r2+1;
    test = round(cR,r2);
end
% Add legend
legend(h,['Real: \mu = ',num2str(round(mean(real),2)),'\pm',...
    num2str(round(cR,r1)),p.Results.unit],...
    ['Permuted: \mu = ',num2str(round(mean(perm),2)),'\pm',...
    num2str(round(cP,r2)),p.Results.unit],...
    ['z = ',num2str(round(z,4))],...
    'Location',p.Results.loc);
% Set figure text (title and axis labels)
title(p.Results.main);
xlabel(p.Results.xlab);
ylabel(p.Results.ylab);
% Convert y-axis into percent (*100)
set(gca,'YTickLabel',cellstr(num2str(get(gca,'YTick')'.*100)))