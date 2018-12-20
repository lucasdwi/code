function [h,p,pAdj] = bulkT(data1,data2,pair,mcc)
% Run bulk t-tests comparing dataset 1 and dataset 2, dataset 1 and 0, or
% between columns in dataset 1
% INPUTS:
% data1 = dataset 1, where rows = subjects and columns = variables
% data2 = dataset 2, where rows = subjects and columns = variables, or 0
% pair = if two datasets, whether or not to use a paired ttest, format = 0
%   or 1, false or true
% mcc = applies multiple comparison correction; options: [] = no
%   correction; 'fdr' = false discovery rate Benjamani-Hochberg method;
%   'bc' = Bonferroni correction

% OUTPUTS:
% h = hypothesis test decision, 1 = reject null that means come from the
%   same distribution
% p = p-value of significance testing
% pAdj = adjusted p-value using 'mcc' method
%% Check that number of columns in data1 and data2 are the same
if ~isempty(data2)
    if size(data1,2) ~= size(data2,2) && data2 ~= 0
        error('The number of colums of data to be compared are not equal.') 
    end
end
%% Run t-tests through columns and get power of test
if isempty(data2)
    % Within dataset two-sample t-test
    pairs = nchoosek(1:size(data1,2),2);
    for c = 1:size(pairs,1)
        [h(:,c),p(:,c)] = ttest(data1(:,pairs(c,1)),data1(:,pairs(c,2)));
    end
else
    for c = 1:size(data1,2)
        if data2 == 0
            % One-sample t-test
            [h(:,c),p(:,c)] = ttest(data1(:,c));
        elseif pair == 0
            % Across datasets two-sample t-test
            [h(:,c),p(:,c)] = ttest2(data1(:,c),data2(:,c));
        elseif pair == 1
            % Across datasets paired t-test
            [h(:,c),p(:,c)] = ttest(data1(:,c),data2(:,c));
        end
    end
end
%% Apply multiple comparisons correction
% FDR
if strcmpi(mcc,'fdr')
   pAdj = mafdr(p','BHFDR','true')'; 
end
% Bonferroni
if strcmpi(mcc,'bc')
   pAdj = p.*size(data1,2); 
end
% No correction
if isempty(mcc)
    pAdj = [];
end