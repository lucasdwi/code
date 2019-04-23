files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\pulseStim\toProcess\','.mat');
for ii = 1:size(files,2)
    load(files{ii})
    intEnd = [eventTs.t{1}(logicFind(7,round(diff(eventTs.t{1})),...
        '==')+1);eventTs.t{4}(end)+6];
    intStart = [eventTs.t{4}(logicFind(7,round(diff(eventTs.t{4})),'=='))...
        ;eventTs.t{4}(end)]+0.01;
    % Convert stim times to intervals between stim; add in values for last
    % window; start and end flip for stim intervals
    for jj = 1:size(intStart,1)
        eventTs.t = [eventTs.t,{intStart(jj)},{intEnd(jj)}];
        eventTs.label = [eventTs.label,['Stim',num2str(jj),' (Start)'],...
            ['Stim',num2str(jj),' (End)']];
    end
    % Add in baseline
    eventTs.t = [eventTs.t,{LFPTs.tvec(1)},{eventTs.t{5}(1)-1}];
    eventTs.label = [eventTs.label,{'Base (Start)'},{'Base (End)'}];
    save(files{ii},'LFPTs','eventTs','pl2','adfreq')
end
%%
c = 1;
for ii = 1:60
   for jj = 1:9
        this(c,:) = params(:,allCombs(jj,1,ii),allCombs(jj,2,ii));
        c = c+1;
   end
end

% for ii = 1:size(trls{2}.sampleinfo,1)
%     idx = nearest_idx3(trls{2}.sampleinfo(ii,1)/400,hist.eventTs.t{1,5});
%     if abs(trls{2}.sampleinfo(ii,1)/400-hist.eventTs.t{1,5}(idx))<1
%         stimIdx(ii) = idx;
%     end 
% end

% Get pop stats for coherence
mCoh = mean(coh{541}.Cxy,4,'omitnan');
sCoh = std(coh{541}.Cxy,[],4,'omitnan');
% Z-score all stim data
zCoh = cellfun(@(x) (x.Cxy-mCoh)./sCoh,coh,'UniformOutput',0);
% Group stim chunks together by params for averaging
allParams = reshape(params,6,9);
stimCoh = cell(1,9);
for ii = 1:9
    inds = logicFind(6,sum(this == allParams(:,ii)',2),'==');
    for jj = 1:size(inds,2)
        if size(zCoh{inds(jj)},4)<12
           zCoh{inds(jj)} = cat(4,zCoh{inds(jj)},...
               nan(size(zCoh{inds(jj)},1),size(zCoh{inds(jj)},2),...
               size(zCoh{inds(jj)},3),12-size(zCoh{inds(jj)},4))); 
        end
        stimCoh{ii} = cat(3,stimCoh{ii},zCoh{inds(jj)});
%         x = logicFind(inds(jj),stimIdx,'==');
%         if ~isempty(x)
%             stimCoh{ii} = cat(3,stimCoh{ii},zCoh(:,:,x));
%         end
    end
end
% Average across stim param chunks
mStimCoh = cellfun(@(x) mean(x,3,'omitnan'),stimCoh,'uniformoutput',0);
%% Grab data for each animal
files{1} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\pulseStim\processed\','IRDM2');
files{2} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\pulseStim\processed\','IRDM5');
files{3} = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\pulseStim\processed\','IRDM6');
for ii = 1:3
    for jj = 1:size(files{ii},2)
        load(files{ii}{jj},'stimCoh')
        allStimCoh{ii}(jj,:) = stimCoh;
    end
    for jj = 1:size(allStimCoh{1},2)
        catCoh{ii,jj} = cat(3,allStimCoh{ii}{:,jj});
    end
end
mCatCoh = cellfun(@(x) mean(x,3,'omitnan'),catCoh,'uniformoutput',0);

for ii = 1:9
    figure
    for jj = 1:3
        subplot(3,1,jj)
        imagesc(squeeze(mCatCoh{jj,ii}(4,:,:,:)))
        caxis([-1 1])
    end
end

mLow = cellfun(@(x) mean(x(:,:,:,3:6),4),cellfun(@(x) ...
    mean(x(:,6:9,:,:),2),mCatCoh,'uniformoutput',0),'uniformoutput',0);

ci = 1;
for c = [4,5,10,11]
    for ii = 1:3
        for jj = 1:9
            p(ii,jj,ci) = mLow{ii,jj}(c);
        end
    end
    ci = ci+1;
end
%% Get average alpha for each param
for ii = 1:9
    c = 1;
    for jj = [4,5,10,11]
        zParam(ii,c) = mean(mStimCoh{ii}(jj,11:15,1,1));
        c = c+1;
    end
end