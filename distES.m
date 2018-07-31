function [d,delt,r] = distES(dist1,dist2)
%% Calculates effect size between two distributions using Wilcoxon Ranksum,
% converting to Mann-Whitney U, then converting it to 1) a r-family
% effect size, and finally a Cohen's d; and 2) Cliff's delta
%
% For Mann-Whitney U, gives the number of times a sample in dist2 precedes
% a sample in dist 2 (i.e, number of times dist2<dist1).
%__________________________________________________________________________
% INPUTS
% dist1 = first distribution; format: vector
% dist2 = second distribution; format: vector
%__________________________________________________________________________
% OUTPUTS
% delt = Cliff's delta; effect size
% d = Cohen' d
% r = correlation coefficient
%__________________________________________________________________________
% LLD 2017
%%
n1 = numel(dist1); n2 = numel(dist2);
% Get Wilcoxon W statistic
[~,~,stats] = ranksum(dist1,dist2);
% Get r-family effect size (Fritz and Morris, 2012)
r = abs(stats.zval/sqrt(n1+n2));
% Get Cohen's d from r (Fritz et al. 2012)
d = (2*r)/sqrt(1-r.^2);
% Convert Wilcoxon W to Mann-Whitney U
u = stats.ranksum-(n1*(n1+1))/2;
% Get Cliff's delta (Cliff 1996)
delt = (2*u)/(n1*n2)-1;