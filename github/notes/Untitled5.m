load('cutoff40Betas.mat')
nameVect = names;
% Only keep first 50 features
nameVect = nameVect(1:50);
%% First get all count
letters = cell(4,40);
for r = 1:4
    these = logicFind(0,threshBeta(r,:),'~=');
    theseName = nameVect(these);
    C = cellfun(@(s) strsplit(s,'-'), theseName,'UniformOutput',false);
    for ci = 1:length(C)
       letters{r,ci} = C{ci}(2); 
    end
    ci = 1;
    for this = {'t','a','b','lg','hg'}
        count(:,ci) = sum(cellfun(@(sc) strcmpi(sc,this),letters),2);
        ci = ci +1;
    end
end
allcount = count;
all = sum(allcount);
%%
tot = sum(allcount,2);
for r = 1:4
    norm(r,:) = allcount(r,:)./tot(r);
end
figure
bar(norm,'stacked'); title('% of Model per Feature')
%% Get power
letters = cell(4,40);
for r = 1:4
    these = logicFind(0,threshBeta(r,1:20),'~=');
    theseName = nameVect(these);
    C = cellfun(@(s) strsplit(s,'-'), theseName,'UniformOutput',false);
    for ci = 1:length(C)
       letters{r,ci} = C{ci}(2); 
    end
    ci = 1;
    for this = {'t','a','b','lg','hg'}
        count(:,ci) = sum(cellfun(@(sc) strcmpi(sc,this),letters),2);
        ci = ci +1;
    end
end
powcount = count';
pow = sum(powcount);
%% Get coh
letters = cell(4,40);
for r = 1:4
    these = logicFind(0,threshBeta(r,21:50),'~=');
    theseName = nameVect(these);
    C = cellfun(@(s) strsplit(s,'-'), theseName,'UniformOutput',false);
    for ci = 1:length(C)
       letters{r,ci} = C{ci}(2); 
    end
    ci = 1;
    for this = {'t','a','b','lg','hg'}
        count(:,ci) = sum(cellfun(@(sc) strcmpi(sc,this),letters),2);
        ci = ci +1;
    end
end
cohcount = count';
coh = sum(cohcount);

%% Interleave pow and coh
allFeat = reshape([powcount(:) cohcount(:)]',2*size(powcount,1),[])';
% Normalize
allFeatNorm = diag(1./sum(allFeat,2))*allFeat;