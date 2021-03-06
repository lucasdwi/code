
files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\lsdStim\toProcess','.mat');
for ii = 2:size(files,2)
    disp(ii)
    load(files{ii})
    % Find first non-empty eventTs.t
    eInd = logicFind(1,cellfun(@(x) ~isempty(x),eventTs.t),'==','first');
    % Grab beginning and end of stim
    stimStart = eventTs.t{eInd}(1);
    stimEnd = eventTs.t{eInd}(end);
    % Find all breaks; add in last one and slight jitter (0.05) to get
    % beyond
    % stim
    stimInds = logicFind(1,round(diff(eventTs.t{eInd}))==10,'==');
    stimIntStart = [eventTs.t{eInd}(stimInds);stimEnd]+0.05;
    % Use the next instance of stim (ind+1) to get ends; add in last end
    stimIntStop = [eventTs.t{eInd}(stimInds+1);stimIntStart(end)+10];
    baseStart = LFPTs.tvec(1);
    baseStop = LFPTs.tvec(nearest_idx3(eventTs.t{eInd}(1),LFPTs.tvec)+1);
    washoutStart = LFPTs.tvec(nearest_idx3(stimIntStop(end),LFPTs.tvec)+1);
    washoutStop = LFPTs.tvec(end);
    eventTs.t = [eventTs.t,{baseStart},{baseStop},{washoutStart},...
        {washoutStop},{stimIntStart},{stimIntStop}];
    eventTs.label = [eventTs.label,'Base (Start)','Base (Stop)',...
        'Washout (Start)','Washout (Stop)','Stim (Start)','Stim (Stop)'];
    save(files{ii},'LFPTs','eventTs','adfreq')
end
%%
[data,files,samps] = collateData('C:\Users\Pythia\Documents\GreenLab\data\lsdStim\processed\',{'IRDM2'},{'pow','coh'},'trl','raw');
%%
for ii = 1:size(data{1},1)
    baseM(ii,:) = mean(data{1,1}{ii,1},1);
    baseS(ii,:) = std(data{1,1}{ii,1},[],1);
    baseZ{ii} = (data{1,1}{ii,1}-baseM(ii,:))./baseS(ii,:);
    stimZ{ii} = (data{1,1}{ii,2}-baseM(ii,:))./baseS(ii,:);
    washZ{ii} = (data{1,1}{ii,3}-baseM(ii,:))./baseS(ii,:);
    baseP(ii,:) = sum(baseZ{ii}>=0,1)/size(baseZ{ii},1);
    stimP(ii,:) = sum(stimZ{ii}>=0,1)/size(stimZ{ii},1);
    washP(ii,:) = sum(washZ{ii}>=0,1)/size(washZ{ii},1);
    for jj = 1:round((size(stimZ{ii},1)-25)/5)
        stimTimeP{ii}(jj,:) = sum(stimZ{ii}(1+(jj-1)*5:25+(jj-1)*5,:)>=0,1)/25;
    end
    for jj = 1:round((size(washZ{ii},1)-25)/5)
        washTimeP{ii}(jj,:) = sum(washZ{ii}(1+(jj-1)*5:25+(jj-1)*5,:)>=0,1)/25;
    end
end
for jj = 1:216
    [~,bb(jj),~] = prop_test([baseP(1,jj)*size(baseZ{1},1),...
        baseP(2,jj)*size(baseZ{2},1)],[size(baseZ{1},1),size(baseZ{2},...
        1)],true);
    [~,bs(jj),~] = prop_test([baseP(1,jj)*size(baseZ{1},1),...
        stimP(1,jj)*size(stimZ{1},1)],[size(baseZ{1},1),size(stimZ{1},...
        1)],true);
    [~,bw(jj),~] = prop_test([baseP(1,jj)*size(baseZ{1},1),...
        washP(1,jj)*size(washZ{1},1)],[size(baseZ{1},1),size(washZ{1},...
        1)],true);
end
%%
[data,samp,files] = collateData('E:\dualSite\processed\all\',{'single';...
    'dual'},{'pow','coh'},'trl','rel');
animals = {'IRDM2','IRDM5','IRDM6','IVSA74','IVSA75','IVSA76'};
for ii = 1:2
    for jj = 1:size(files{ii},2)
        parts = strsplit(files{ii}{jj},'_');
        names{ii,jj} = parts{1};
    end
end
for ii = 1:size(animals,2)
    singleInds = logicFind(1,cellfun(@(x) isequal(x,animals{ii}),...
        names(1,:)),'==');
    dualInds = logicFind(1,cellfun(@(x) isequal(x,animals{ii}),...
        names(2,:)),'==');
    base{ii} = cat(1,data{1,1}{singleInds,1},data{1,2}{dualInds,1});
end
[lsdData,~,lsdFiles] = collateData(['C:\Users\Pythia\Documents\GreenLab'...
    '\data\lsdStim\processed\'],{'IRDM2';'IRDM5';'IRDM6'},{'pow','coh'},...
    'trl','rel');
for ii = 1:3
    thisZ{ii} = cellfun(@(x) (x-mean(base{ii},1))./std(base{ii},[],1),...
        lsdData{ii},'UniformOutput',false);
    over{ii} = cellfun(@(x) sum(x>=0,1),thisZ{ii},'UniformOutput',false);
    n{ii} = cellfun(@(x) size(x,1),thisZ{ii});
    for jj = 1:size(over{ii},1)
        for k = 1:size(over{ii},2)
            pro{ii}(jj,k,:) = over{ii}{jj,k}./n{ii}(jj,k);
        end
    end
end

%% ttest 
for ii = 1:3
    for jj = 1:size(lsdData{ii},1)
        for k = 1:size(lsdData{ii},2)
            for m = 1:216
                [~,p{ii}(jj,k,m)] = ttest2(base{ii}(:,m),...
                    lsdData{ii}{jj,k}(:,m));
            end
        end
    end
end
pAdj = cellfun(@(x) x.*numel(x),p,'uniformoutput',0);
sig = cellfun(@(x) x<=0.05,pAdj,'UniformOutput',0);
% For each animal and each recording count number of significant features

%% Plot mean z-value if p-value <= 0.05
for ii = 1:3
    figure
    this = [mean(thisZ{ii}{1,1},1);mean(thisZ{ii}{1,2},1);...
        mean(thisZ{ii}{1,3})];
    pcolor(this)
end

files = fileSearch('C:\Users\Pythia\Documents\GreenLab\data\lsdStim\toProcess','.mat');
for ii = 2:size(files,2)
    disp(ii)
    load(files{ii})
    % Find first non-empty eventTs.t
    eInd = logicFind(1,cellfun(@(x) ~isempty(x),eventTs.t),'==','first');
    % Grab beginning and end of stim
    stimStart = eventTs.t{eInd}(1);
    stimEnd = eventTs.t{eInd}(end);
    % Find all breaks; add in last one and slight jitter (0.05) to get
    % beyond
    % stim
    stimInds = logicFind(1,round(diff(eventTs.t{eInd}))==10,'==');
    stimIntStart = [eventTs.t{eInd}(stimInds);stimEnd]+0.05;
    % Use the next instance of stim (ind+1) to get ends; add in last end
    stimIntStop = [eventTs.t{eInd}(stimInds+1);stimIntStart(end)+10];
    baseStart = LFPTs.tvec(1);
    baseStop = LFPTs.tvec(nearest_idx3(eventTs.t{eInd}(1),LFPTs.tvec)+1);
    washoutStart = LFPTs.tvec(nearest_idx3(stimIntStop(end),LFPTs.tvec)+1);
    washoutStop = LFPTs.tvec(end);
    eventTs.t = [eventTs.t,{baseStart},{baseStop},{washoutStart},...
        {washoutStop},{stimIntStart},{stimIntStop}];
    eventTs.label = [eventTs.label,'Base (Start)','Base (Stop)',...
        'Washout (Start)','Washout (Stop)','Stim (Start)','Stim (Stop)'];
    save(files{ii},'LFPTs','eventTs','adfreq')
end
%%
[data,files,samps] = collateData('C:\Users\Pythia\Documents\GreenLab\data\lsdStim\processed\',{'IRDM2'},{'pow','coh'},'trl','raw');
%%
for ii = 1:size(data{1},1)
    baseM(ii,:) = mean(data{1,1}{ii,1},1);
    baseS(ii,:) = std(data{1,1}{ii,1},[],1);
    baseZ{ii} = (data{1,1}{ii,1}-baseM(ii,:))./baseS(ii,:);
    stimZ{ii} = (data{1,1}{ii,2}-baseM(ii,:))./baseS(ii,:);
    washZ{ii} = (data{1,1}{ii,3}-baseM(ii,:))./baseS(ii,:);
    baseP(ii,:) = sum(baseZ{ii}>=0,1)/size(baseZ{ii},1);
    stimP(ii,:) = sum(stimZ{ii}>=0,1)/size(stimZ{ii},1);
    washP(ii,:) = sum(washZ{ii}>=0,1)/size(washZ{ii},1);
    for jj = 1:round((size(stimZ{ii},1)-25)/5)
        stimTimeP{ii}(jj,:) = sum(stimZ{ii}(1+(jj-1)*5:25+(jj-1)*5,:)>=0,1)/25;
    end
    for jj = 1:round((size(washZ{ii},1)-25)/5)
        washTimeP{ii}(jj,:) = sum(washZ{ii}(1+(jj-1)*5:25+(jj-1)*5,:)>=0,1)/25;
    end
end
for jj = 1:216
    [~,bb(jj),~] = prop_test([baseP(1,jj)*size(baseZ{1},1),...
        baseP(2,jj)*size(baseZ{2},1)],[size(baseZ{1},1),size(baseZ{2},...
        1)],true);
    [~,bs(jj),~] = prop_test([baseP(1,jj)*size(baseZ{1},1),...
        stimP(1,jj)*size(stimZ{1},1)],[size(baseZ{1},1),size(stimZ{1},...
        1)],true);
    [~,bw(jj),~] = prop_test([baseP(1,jj)*size(baseZ{1},1),...
        washP(1,jj)*size(washZ{1},1)],[size(baseZ{1},1),size(washZ{1},...
        1)],true);
