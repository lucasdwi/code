fNames = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\mike\toProcess\');
%% Get stim times
stimInd = []; stim = []; stimStop = []; stimStart = [];
stimInd = logicFind(1,~cellfun(@isempty,eventTs.t),'==');
stim = sort(cat(1,eventTs.t{stimInd(1)},eventTs.t{stimInd(2)}));
% Get all stops
stimStop = stim(logicFind(0,round(diff(stim)),'~='));
% Add last stop since above algorithm skips
stimStop = [stimStop;stim(end)];
% Get starts using the next index after stops from stim
stimStart = stim(logicFind(0,round(diff(stim)),'~=')+1);
% Add first start since above algorithm skips
stimStart = [stim(1);stimStart];
% Use stim start and stops to determine data start and stops; also include
% first/last data point - assumes that no stim starts or ends a recording
dataStart = [LFPTs.tvec(1);LFPTs.tvec(nearest_idx3(stimStop,LFPTs.tvec)+1)'];
dataStop = [LFPTs.tvec(nearest_idx3(stimStart,LFPTs.tvec)-1)';LFPTs.tvec(end)];
% Store old eventTs
eventTsOld = eventTs;
% Build new eventTs
eventTs = [];
% eventTs2.label = {'Data (Start)','Data (End)','Stim (Start)','Stim (End)'};
eventTs.label = {'Base (Start)','Base (End)','Int1 (Start)','Int1 (Stop)','Int2 (Start)','Int2 (Stop)','Int3 (Start)','Int3 (Stop)','Int4 (Start)','Int4 (Stop)','Int5 (Start)','Int5 (Stop)','Post (Start)','Post (Stop)'};
eventTs.t = [];
for ii = 1:length(dataStart)
    eventTs.t{ii*2-1} = dataStart(ii);
    eventTs.t{ii*2} = dataStop(ii);
end
%% Check by plotting
figure
plot(stim,ones(length(stim),1),'.')
hold on
plot(stimStart,ones(length(stimStart),1),'o')
plot(stimStop,ones(length(stimStop),1),'o')
plot(dataStart,ones(length(dataStart),1),'s')
plot(dataStop,ones(length(dataStop),1),'s')
plot(LFPTs.tvec,LFPTs.data(1,:))
%% Go through all channels checking for clean ones
data.cfg.history.mfun = [];
for ii = 1:16
    data.data = LFPTs.data(ii,:);
    dummy = threshFilt(data,2,0.0125,5,5,adfreq);
    filt(ii,:) = dummy.data;
end
 new = [];
% Grab NAcc shell
% new = LFPTs.data(8,:);
% Find best NAcc core
naccInd = [7,9,10];
[~,minInd] = min(sum(isnan(LFPTs.data(naccInd,:)),2));
new = [new;LFPTs.data(naccInd(minInd),:)];
% Find best mPFC
mpfcInd = [3,4,5,6,11,12];
[~,minInd] = min(sum(isnan(LFPTs.data(mpfcInd,:)),2));
new = [new;LFPTs.data(mpfcInd(minInd),:)];
% Find best OFC
% ofcInd = [2,13,14];
% [~,minInd] = min(sum(isnan(LFPTs.data(ofcInd,:)),2));
% new = [new;LFPTs.data(ofcInd(minInd),:)];
% Find best AIC
% aicInd = [1,15,16];
% [~,minInd] = min(sum(isnan(LFPTs.data(aicInd,:)),2));
% new = [new;LFPTs.data(aicInd(minInd),:)];
%% Check by plotting
figure
hold on
for ii = 1:size(new,1)
    subplot(1,size(new,1),ii)
    plot(new(ii,:)) 
end
%% Rearrange and save
LFPTsOld = LFPTs;
LFPTs.data = new;
LFPTs.label = {'NAccCore','mPFC'}; %,'mPFC','OFC','AIC'};
% save('','adfreq','eventTs','eventTsOld','LFPTs','LFPTsOld','-v7.3')
%%
c = 1;
for ii = 2:3
    zCoh{c} = (coh{1,ii}.Cxy-repmat(coh{1,1}.mCxy,1,1,size(coh{1,ii},3)))./repmat(coh{1,1}.sdCxy,1,1,size(coh{1,ii}.Cxy,3));
    zPow{c} = (psdTrls{1,ii}.Pow-repmat(mean(psdTrls{1,1}.Pow,3),1,1,size(psdTrls{1,ii}.Pow,3)))./repmat(std(psdTrls{1,1}.Pow,[],3),1,1,size(psdTrls{1,ii}.Pow,3));
    zRelPow{c} = (psdTrls{1,ii}.relPow-repmat(psdTrls{1,1}.avgRelPow,1,1,size(psdTrls{1,ii}.relPow,3)))./repmat(psdTrls{1,1}.stdRelPow,1,1,size(psdTrls{1,ii}.relPow,3));
    zRelCoh{c} = (coh{1,ii}.rel-repmat(mean(coh{1,1}.rel,3),1,1,size(coh{1,ii}.rel,3)))./repmat(std(coh{1,1}.rel,[],3),1,1,size(coh{1,ii}.rel,3));
    c = c+1;
end
%%
for jj = 1:2
    figure
    for ii = 1:2
        subplot(1,2,ii)
        h = pcolor(padarray(squeeze(zRelPow{jj}(:,ii,:)),[1,0],'post'));
        set(h,'EdgeColor','none')
        set(gca,'ydir','reverse')
        colormap('viridis')
    end
end
%%
test = squeeze(zRelCoh{1});
% test(abs(test)<=2) = 0;
figure
h = pcolor(padarray(test,[1,0],'post'));
set(h,'EdgeColor','none')
set(gca,'ydir','reverse')
colormap('viridis')
caxis([-4 4])
%%
figure
time = LFPTs.tvec(trls{2}.sampleinfo(:,1))-eventTs.t{3};
subplot(1,2,1)
imagesc(time,1:2:100,squeeze(zPow{1}(2,:,:)))
subplot(1,2,2)
time = LFPTs.tvec(trls{3}.sampleinfo(:,1))-eventTs.t{13};
imagesc(time,1:2:100,squeeze(zPow{2}(2,:,:)))

figure
imagesc(time,1:2:100,squeeze(zCoh{1}(1,:,:)))
figure
imagesc(time,1:2:100,squeeze(zCoh{2}(1,:,:)))
%%
for ii = 1:4
    figure
    for jj = 2:7
        subplot(2,3,jj-1)
        time = LFPTs.tvec(trls{jj}.sampleinfo(:,1))-eventTs.t{jj*2-1};
        imagesc(time,1:2:100,squeeze(zPow{jj-1}(ii,:,:)))
        caxis([-10 10])
    end
end
%%
for ii = 1:6
    figure
    for jj = 2:7
       subplot(2,3,jj-1)
       time = LFPTs.tvec(trls{jj}.sampleinfo(:,1))-eventTs.t{jj*2-1};
       imagesc(time,1:2:100,squeeze(zCoh{jj-1}(ii,:,:)))
       caxis([-10 10])
    end
end
