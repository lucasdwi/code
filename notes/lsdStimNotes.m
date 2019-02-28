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
