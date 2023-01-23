%% Build models telling apart two recordings ~14 days apart with increasing number of animals
%% First go through files and find pairs of recordings 14 days apart with 
% minimal difference in AUC
uID = unique(thisData(:,1));
thisDiff = nan(72,122);
for ii = 1:numel(uID)
    disp(ii)
    inds = logicFind(1,cellfun(@(x) contains(x,uID{ii}),thisData(:,1)),'==');
    for jj = 1:numel(inds)
        thisDate = datetime(thisData{inds(jj),2});
        nextDate = thisDate+days(14);
        check = cellfun(@(x) nextDate==datetime(x),thisData(inds,2)) & ...
            logical(cell2mat(cellfun(@(x) ~contains(x,...
            {'MPH','OFC','ATM','GUAN','Shell','Core','GFC','IL'}),...
            thisData(inds,4),'UniformOutput',false)));
        if any(check)
            option{ii,jj} = [inds(jj) inds(check)];
            thisDiff(ii,jj) = thisData{inds(jj),3}-...
                thisData{inds(check),3};
        else
            option{ii,jj} = [];
        end
    end
    minInd = nearest_idx2(0,thisDiff(ii,:));
    if ~isempty(option{ii,minInd})
        best(ii,:) = option{ii,minInd};
    end
end
for ii = 1:size(best,1)
    if best(ii,1) ~= 0
        these(ii,:) = [{strcat(thisData{best(ii,1),1},'_',thisData{best(ii,1),2})},...
        {strcat(thisData{best(ii,2),1},'_',thisData{best(ii,2),2})}];
    end
end
%%  load data
[dataRaw,samps,files] = collateData('F:/irdmRound2/sampleSize/',{'.mat'},...
    {'pow','coh'},'trl','');
[dataRel,~,~] = collateData('F:/irdmRound2/sampleSize/',{'.mat'},...
    {'pow','coh'},'trl','rel');
%% 2
for ii = 1:210
    load(['F:/irdmRound2/sampleSize/sampleSize2/sampleSize2_',num2str(ii),'.mat'])
    raw2(ii) = aRaw{1}.acc;
    rel2(ii) = aRel{1}.acc;
end
%% 3
for ii = 1:1330
    load(['F:/irdmRound2/sampleSize/sampleSize3/sampleSize3_',num2str(ii),'.mat'])
    raw3(ii) = aRaw{1}.acc;
    rel3(ii) = aRel{1}.acc;
end