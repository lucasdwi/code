function [h,p,pAdj] = bulkT(data1,data2,mcc)
% Run bulk t-tests comparing dataset 1 and dataset 2 or dataset 1 and 0
% INPUTS:
% data1 = dataset 1, where rows = subjects and columns = variables
% data2 = dataset 2, where rows = subjects and columns = variables, or 0
% mcc = applies multiple comparison correction; options: [] = no
%   correction; 'fdr' = false discovery rate Benjamani-Hochberg method;
%   'bc' = Bonferroni correction

% OUTPUTS:
% h = hypothesis test decision, 1 = reject null that means come from the
%   same distribution
% p = p-value of significance testing
% power = power of test
%% Check that number of columns in data1 and data2 are the same
if size(data1,2) ~= size(data2,2) && data2 ~= 0
   error('The number of colums of data to be compared are not equal.') 
end
%% Run t-tests through columns and get power of test
for c = 1:size(data1,2)
    if data2 == 0
        % One-sample t-test
        [h(:,c),p(:,c)] = ttest(data1(:,c));
    else
        % Two-sample t-test
        [h(:,c),p(:,c)] = ttest(data1(:,c),data2(:,c));
    end
    %power(:,c) = sampsizepwr();
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