function [newPow,newCoh] = addNaNs(LFPTs,pow,coh,nChan)
for ii = 1:32
    if ii < 10
        chanNames{ii} = ['FP0',num2str(ii)];
    else
        chanNames{ii} = ['FP',num2str(ii)];
    end
end
chanArray = reshape(chanNames,nChan,4)';
% Find first channel to get which box/row
[r,~] = find(contains(chanArray,LFPTs.label{1}));
% Use row to compare to known chanArray
missing = logicFind(0,contains(chanArray(r,:),LFPTs.label),'==');
% Add in NaNs where missing data are
newPow = zeros(6,nChan,size(pow,3));
% First power; just use channel indices
c = 1;
for ci = 1:nChan
    if any(ci == missing)
        newPow(:,ci,:) = deal(NaN);
    else
        newPow(:,ci,:) = pow(:,c,:);
        c = c+1;
    end
end
% Then coh; first figure out which cmbs are missing
cmbs = nchoosek(1:nChan,2);
newCoh = zeros(size(cmbs,1),6,size(coh,3));
theseCoh = [];
for m = 1:numel(missing)
    theseCoh(:,m) = any(cmbs==missing(m),2);
end
skipCmbs = logicFind(1,any(theseCoh,2),'==');
c = 1;
for ci = 1:size(cmbs,1)
    if any(ci == skipCmbs)
        newCoh(ci,:,:) = deal(NaN);
    else
        newCoh(ci,:,:) = coh(c,:,:);
        c = c+1;
    end
end