stimOff = 10;
stimOn = 50;
base = 15*60;
%% Get timestamps of baselines
stim = eventTs.t{1,1};
baseStart = eventTs.t{1}(logicFind(base,round(diff(stim),-2),'==')-1);
baseStop = eventTs.t{1}(logicFind(base,round(diff(stim),-2),'==')+1);
% Add in first and last which are skipped by above heuristic
baseStart = [LFPTs.tvec(1);baseStart;LFPTs.tvec(nearest_idx3(stim(end),LFPTs.tvec)+1)];
baseStop = [LFPTs.tvec(nearest_idx3(stim(1),LFPTs.tvec)-1);baseStop;LFPTs.tvec(end)];
% Add baselines to struct
struct.t = [];
struct.label = [];
for ii = 1:size(baseStart,1)
    struct.label = [struct.label,{['Base',num2str(ii),' (Start)']},{['Base',num2str(ii),' (End)']}];
    struct.t = [struct.t,{baseStart(ii)},{baseStop(ii)}];
end
%% Get timestamps of data during stim using timestamps of baselines as references
for ii = 1:size(baseStart,1)-1
    thisStartInd = nearest_idx3(LFPTs.tvec(nearest_idx3(baseStop(ii),LFPTs.tvec)-1),stim); 
    thisStopInd = nearest_idx3(LFPTs.tvec(nearest_idx3(baseStart(ii+1),LFPTs.tvec)+1),stim);
    % Add one second to start time to avoid stim noise
    thisInter = stim(thisStartInd:thisStopInd);
    interStart = thisInter(logicFind(stimOff,round(diff(thisInter)),'=='))+1;
    interStop = thisInter(logicFind(stimOff,round(diff(thisInter)),'==')+1);
    % Add last interval (first 10 seconds of next baseline)
    interStart = [interStart;baseStart(ii+1)+1];
    interStop = [interStop;baseStart(ii+1)+11];
    % Combine all starts and stops into structure
    struct.label = [struct.label,{['Inter',num2str(ii),' (Start)']},{['Inter',num2str(ii),' (End)']}];
    struct.t = [struct.t,{interStart},{interStop}];
end
%% Replace old eventTs with new struct
oldEventTs = eventTs;
eventTs = struct;
%%
save('ZA1_dualSiteDelay_3Hz_250A_175A_2018-05-18.mat','LFPTs','eventTs','oldEventTs','adfreq','pl2','-v7.3')
%%
c = 1;
b = 1;
for ii = 3
    zInterPow{c} = (psdTrls{ii}.Pow-mean(psdTrls{c}.Pow(:,:,end-60:end),3))./std(psdTrls{c}.Pow(:,:,end-60:end),[],3);
    zInterCoh{c} = (coh{ii}.Cxy-mean(coh{c}.Cxy(:,:,end-60:end),3))./std(coh{c}.Cxy(:,:,end-60:end),[],3);
%     b = b+1;
       c = c+1;
end
%%
figure
imagesc(squeeze(zInterPow{1}(5,:,:)))
colormap viridis
caxis([-3 3])
title('Post NAc: R-NAc Power')

figure
imagesc(squeeze(zInterPow{2}(5,:,:)))
colormap viridis
caxis([-3 3])
title('Post OFC: R-NAc Power')
%% Loop through all power
for ii = 1:size(zInterPow{1},1)
    figure
    imagesc(squeeze(zInterPow{1}(ii,:,:)))
    colormap viridis
    caxis([-3 3])
    % title('Post NAc: R-NAc Power')
end
%% Loop through all coherence
main = nchoosek({'PR','PL','OR','OL','NR','NL','LR','LL'},2);
for ii = 1:size(zInterCoh{1},1)
    figure
    imagesc(squeeze(zInterCoh{1}(ii,:,:)))
    colormap viridis
    caxis([-3 3])
    title([main{ii,1},'-',main{ii,2}])
end
%%
figure
imagesc(squeeze(zInterPow{1}(5,:,:)))
colormap viridis
caxis([-3 3])
title('Post: NAc R Power')

figure
imagesc(squeeze(zInterPow{2}(5,:,:)))
colormap viridis
caxis([-3 3])
title('Inter: NAc R Power')
%%
figure
imagesc(1:size(zInterCoh{1},3),coh{1}.f,squeeze(zInterCoh{1}(23,:,:)))
colormap viridis
caxis([-3 3])
title('Post NAc: NAc L-R Coherence')

figure
imagesc(1:size(zInterCoh{2},3),coh{1}.f,squeeze(zInterCoh{2}(23,:,:)))
colormap viridis
caxis([-3 3])
title('Post OFC: NAc L-R Coherence')
%%
figure
c = 1;
for ii = [3,1,2]
    subplot(1,3,c)
    imagesc(squeeze(zInterPow{ii}(6,:,:)))
%     imagesc(squeeze(zInterCoh{ii}(4,:,:)))
    colormap viridis
    caxis([-4 4])
    c = c+1;
end
%%
% load('ZA1_dualSiteDelay_230A_3Hz_2018-05-11_Base1.mat')
%%
c = 1;
for ii = 11:19
    zInterPow{c} = (psdTrls{ii}.Pow-mean(psdTrls{c}.Pow(:,:,end-60:end),3))./std(psdTrls{c}.Pow(:,:,end-60:end),[],3);
    zInterCoh{c} = (coh{ii}.Cxy-mean(coh{c}.Cxy(:,:,end-60:end),3))./std(coh{c}.Cxy(:,:,end-60:end),[],3);
    c = c+1;
end
c = 1;
for ii = 11:19
    percInterPow{c} = (psdTrls{ii}.Pow-mean(psdTrls{c}.Pow(:,:,end-60:end),3))./mean(psdTrls{c}.Pow(:,:,end-60:end),3);
    percInterCoh{c} = (coh{ii}.Cxy-mean(coh{c}.Cxy(:,:,end-60:end),3))./mean(coh{c}.Cxy(:,:,end-60:end),3);
    c = c+1;
end
for ii = 2:11
   baseCoh{ii-1} = (coh{ii}.Cxy-mean(coh{ii-1}.Cxy(:,:,end-60:end),3))./mean(coh{ii-1}.Cxy(:,:,end-60:end),3);
end
%%
mains = {'0','+5','+7','+10','+15','-5','-7','-10','-15'};
figure
c = 1;
for ii = [9,8,7,6,1,2,3,4,5]
% for ii = [2,1,3]
    subplot(2,5,c)
%     imagesc((1:size(percInterCoh{ii},3)).*5,coh{1}.f,squeeze(percInterCoh{ii}(11,:,:)))
%     imagesc((1:size(zInterCoh{ii},3)).*5,coh{1}.f,squeeze(zInterCoh{ii}(11,:,:)))
    imagesc(squeeze(zInterPow{ii}(6,:,:)))
    colormap viridis
    caxis([-3 4])
    title(mains{ii})
%     ylim([0 20])
%     xlabel('Time (sec)')
%     ylabel('Frequency (Hz)')
    set(gca,'YDir','normal')
    box off
    c = c+1;
end
%% base
figure
c = 1;
for ii = [9,8,7,6,1,2,3,4,5]
    subplot(2,5,c)
    imagesc((1:size(baseCoh{ii},3)).*5,coh{1}.f,squeeze(baseCoh{ii}(11,:,:)))
    colormap viridis
    caxis([-3 4])
    ylim([0 20])
    set(gca,'YDir','normal')
    title(mains{ii})
    box off
    c = c+1;
end
%%
mDelta = cellfun(@(x) squeeze(mean(mean(x(11,1:20,:)),3)),zInterCoh,'UniformOutput',0);
% mDelta = cellfun(@(x) squeeze(mean(mean(x(4,4:20,:)),3)),percInterCoh,'UniformOutput',0);
figure
plot(1:9,cell2mat(mDelta([9,8,7,6,1,2,3,4,5])),'o')
ylabel('\DeltaCoherence')
xlabel('\Deltat')
set(gca,'XTick',1:9,'XTickLabel',{-15,-10,-7,-5,0,5,7,10,15})
box off
