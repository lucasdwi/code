function [allCmbs] = splitCmbs(group1,group2)
%%
cmbs1 = [];
if isodd(numel(group1))
    cmbsUp = nchoosek(group1,ceil(numel(group1)/2));
    cmbsDown = nchoosek(group1,floor(numel(group1)/2));
    c = 1;
    for ii = 1:size(cmbsUp,1)
        cmbs1{c,1} = cmbsUp(ii,:);
        c = c+1;
    end
    for ii = 1:size(cmbsDown,1)
        cmbs1{c,1} = cmbsDown(ii,:);
        c = c+1;
    end
else
    cmbs = nchoosek(group1,ceil(numel(group1)/2));
    for ii = 1:size(cmbs,1)
        cmbs1{ii,1} = cmbs(ii,:);
    end
end
%%
cmbs2 = [];
if isodd(numel(group2))
    cmbsUp = nchoosek(group2,ceil(numel(group2)/2));
    cmbsDown = nchoosek(group2,floor(numel(group2)/2));
    c = 1;
    for ii = 1:size(cmbsUp,1)
        cmbs2{c,1} = cmbsUp(ii,:);
        c = c+1;
    end
    for ii = 1:size(cmbsDown,1)
        cmbs2{c,1} = cmbsDown(ii,:);
        c = c+1;
    end
else
    cmbs = nchoosek(group2,ceil(numel(group2)/2));
    for ii = 1:size(cmbs,1)
        cmbs2{ii,1} = cmbs(ii,:);
    end
end
%%
allCmbs = [];
c = 1;
for ii = 1:size(cmbs1,1)
    for jj = 1:size(cmbs2,1)
        allCmbs{c,1} = [cmbs1{ii,:},cmbs2{jj,:}];
        c = c+1;
    end
end