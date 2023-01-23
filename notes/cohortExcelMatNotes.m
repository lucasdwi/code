%% Extract data from IRDM cohort spreadsheets (helps automates checking AUCs)
%% Import spreadsheets
files = fileSearch(['H:\Shared drives\Green Lab\',...
    'IRDM Data + Cohort Spreadsheets\Cohort spreadsheets']...
    ,'.xlsx','in','bTBI','ex');
thisData = cell(1,4);
c = 1;
for ii = 1:numel(files)
    [~,sheetName] = xlsfinfo(files{ii});
    for jj = 1:numel(sheetName)
        if ~contains(sheetName{jj},'Summary') && ...
                ~contains(sheetName{jj},'Bias')
            thisSheet = readtable(files{ii},'sheet',sheetName{jj});
            % Get indices of rows that are not: 'no pl2', 'pre-implant',
            % don't have date, 'skip'
            inds = logicFind(1,...
                cellfun(@isempty,strfind(thisSheet.notes,'no pl2')) & ...
                cellfun(@isempty,strfind(thisSheet.notes,'pre-implant'))...
                & ~isnat(thisSheet.Date) & ...
                cellfun(@isempty,strfind(thisSheet.notes,'skip')) & ...
                cellfun(@isempty,strfind(thisSheet.notes,'T3')) & ...
                cellfun(@(x) ~isempty(x),thisSheet.notes),'==');
            for k = inds
                if isdatetime(thisSheet{k,1}) && ~isnan(thisSheet{k,2})
                    thisData(c,:) = {sheetName{jj},...
                        char(datetime(thisSheet{k,1},'format',...
                        'yyyy-MM-dd')),thisSheet.AUC(k),...
                        thisSheet.notes(k)};
                    % Determine delay
                    if contains(thisSheet{k,"notes"},'32')
                        thisDelay{c} = [0,4,8,16,32];
                    elseif contains(thisSheet{k,"notes"},'10')
                        thisDelay{c} = [0,2.5,5,7.5,10];
                    elseif contains(thisSheet{k,"notes"},'12')
                        thisDelay{c} = [0,1,3,6,12];
                    elseif contains(thisSheet{k,"notes"},'60')
                        thisDelay{c} = [0,15,30,45,60];
                    end
                    % Grab percent delay choice
                    percentDelay(c,:) = thisSheet{k,[3,5,7,9,11]}./100;
                    trapzAUCs(c) = trapz(thisDelay{c},percentDelay(c,:));
                    thisFile(c,:) = thisData(c,1:2);
                    c = c+1;
                end
            end
        end
    end
end
% Fill in NaNs in empty spots of thisFile so that it can be searched
thisFile(cellfun(@isempty,thisFile(:,1)),1) = deal({'nemo'});
thisFile(cellfun(@isempty,thisFile(:,2)),2) = deal({'nowhen'});
%% Check effect of intervention
% Number of baselines to compare to
n = 5;
% Interventions to look for
inter = {'IL','OFC','shell','core','MPH3','MPH1','ATM'};
for ii = 1:numel(inter)
    interInd = logicFind(1,cellfun(@(x) contains(x,inter(ii),...
        'IgnoreCase',true),thisData(:,4)),'==');
    for jj = interInd
        disp([inter{ii},': ',num2str(logicFind(jj,interInd,'==')),...
            ' of ',num2str(numel(interInd))])
        thisDate = thisData{jj,2};
        thisID = thisData{jj,1};
        % Find baseline days from this animal before intervention
        baseInds = logicFind(1,cellfun(@(x) contains(x,'base'),...
            thisData(:,4)) & cellfun(@(x) ...
            contains(x,thisID),thisData(:,1)) & cellfun(@(x) ...
            datetime(x)<thisDate,thisData(:,2)),'==');
        mBase = mean(cell2mat(thisData(baseInds(end-n+1:end),3)));
        sBase = std(cell2mat(thisData(baseInds(end-n+1:end),3)));
        if cell2mat(thisData(jj,3))>mBase+2*sBase
            thisData{jj,5} = [1,mBase,sBase];
        elseif cell2mat(thisData(jj,3))<mBase-2*sBase
            thisData{jj,5} = [-1,mBase,sBase];
        else
            thisData{jj,5} = [0,mBase,sBase];
        end
    end
end
%% count number of animals with recordings for each day after intervention
inter = {'sIL_','sCore_','sOFC_','MPH3_'};
uID = unique(thisData(:,1));
postData = cell(numel(inter),7);
procFiles1 = fileSearch('D:\DDT\processed\','.mat');
procFiles2 = fileSearch('E:\processed\','.mat');
procFiles = cell(1,1);
c = 1;
for ii = 1:size(procFiles1,2)
    procFiles{c,1} = ['D:\DDT\processed\',procFiles1{ii}];
    c = c+1;
end
for ii = 1:size(procFiles2,2)
    procFiles{c,1} = ['E:\processed\',procFiles2{ii}];
    c = c+1;
end
files = [];
missing = [];
c = 1;
for ii = 2:numel(inter)
    for jj = 1:numel(uID)
        parts = strsplit(uID{jj},'_');
        interInd = logicFind(1,cellfun(@(x) contains(x,inter(ii),...
            'IgnoreCase',true),thisData(:,4)) & cellfun(@(x) ...
            contains(x,uID(jj)),thisData(:,1)),'==');
        if ~isempty(interInd)
            % Find last intervention
            lastInterDate = datetime(thisData(interInd(end),2));
            loadInd = logicFind(1,cellfun(@(x) contains(x,parts(1)),...
                procFiles) & cellfun(@(x) contains(x,...
                string(lastInterDate,'yyyy-MM-dd')),...
                procFiles),'==');
            if ~isempty(loadInd)
            disp(procFiles{loadInd})
            load(procFiles{loadInd},'psdTrls','coh')
            [b,c,t] = size(psdTrls{1}.bandPow);
            thisPow = reshape(psdTrls{1}.bandPow,b*c,t)';
            [cmb,b,t] =  size(coh{1}.mBandCoh);
            thisCoh = reshape(permute(...
                squeeze(coh{1}.mBandCoh),[2,1,3]),...
                cmb*b,t)';
            lastData{ii,jj} = [thisPow,thisCoh];
            clear thisPow thisCoh

            files = [files;{thisData{interInd(end),1},...
                thisData{interInd(end),2}}];
            % Get baselines
            theseBase = logicFind(1,cellfun(@(x) contains(x,'base'), ...
                thisData(:,4)) & cellfun(@(x) ...
                contains(x,uID(jj)),thisData(:,1)),'==');
            % Find 3 prior baselines
            preInds = theseBase(logicFind(1,...
                datetime(thisData(theseBase,2))<lastInterDate,'=='));
            for k = 0:2
                ind = logicFind(1,contains(procFiles,parts(1)) & ...
                    contains(procFiles,thisData(preInds(end-k),2)),'==');
                if isempty(ind)
                    missing = [missing;parts(1),...
                        thisData(preInds(end-k),2)];
                end
                disp(procFiles{ind})
                load(procFiles{ind},'psdTrls','coh')
                [b,c,t] = size(psdTrls{1}.bandPow);               
                thisPow = reshape(psdTrls{1}.bandPow,b*c,t)';
                [cmb,b,t] =  size(coh{1}.mBandCoh);
                thisCoh = reshape(permute(...
                    squeeze(coh{1}.mBandCoh),[2,1,3]),...
                    cmb*b,t)';
                baseData{ii,jj}{k+1} = [thisPow,thisCoh];
                clear thisPow thisCoh
            end
            files = [files;thisData(preInds(end-2:end),1),...
                thisData(preInds(end-2:end),2)];
            % Find baseline recordings up to 1 day after last intervention
%             postInds = theseBase(logicFind(1,days(datetime( ...
%                 thisData(theseBase,2))-lastInterDate)<=1 & days(datetime( ...
%                 thisData(theseBase,2))-lastInterDate)>0,'=='));
%             if ~isempty(postInds)
%                 ind = logicFind(1,contains(procFiles,parts(1)) & ...
%                     contains(procFiles,thisData(postInds,2)),'==');
%                 if isempty(ind)
%                     missing = [missing;parts(1),thisData(postInds,2)];
%                 end
%                 disp(procFiles{ind})
%                 load(procFiles{ind},'psdTrls','coh')
%                 [b,c,t] = size(psdTrls{1}.bandPow);               
%                 thisPow = reshape(psdTrls{1}.bandPow,b*c,t)';
%                 [cmb,b,t] =  size(coh{1}.mBandCoh);
%                 thisCoh = reshape(permute(...
%                     squeeze(coh{1}.mBandCoh),[2,1,3]),...
%                     cmb*b,t)';
%                 postData{ii,jj}{k+1} = [thisPow,thisCoh];
%                 clear thisPow thisCoh
%             end
%             files = [files;thisData(postInds,1),...
%                 thisData(postInds,2)];
%             if ~isempty(postInds)
%                 for k = 1:7
%                     if any(days(thisData(postInds,2)-lastInterDate) == k)
%                         postData{ii,k} = [postData{ii,k};uID(jj)];
%                     end
%                 end
%             end
            else
                disp([parts{1},string(lastInterDate),': skipped'])
            end
        end
    end
end
%%
inter = {'IL','OFC','shell','core','MPH3','MPH1','ATM'};
for ii = 1:numel(inter)
    interInd = logicFind(1,cellfun(@(x) contains(x,inter(ii),...
        'IgnoreCase',true),thisData(:,4)),'==');
    effect{ii} = cat(1,thisData{interInd,5});
    figure
    histogram(effect{ii}(:,1))
    title(inter{ii})
end
%% Plot trapz AUCs
allDelay = cat(1,thisDelay{:});
u = unique(allDelay,'rows');
figure
for ii = 1:size(u,1)
    this = [];
    for jj = 1:size(allDelay,1)
        if sum(allDelay(jj,:) == u(ii,:)) == 5
            this = [this;percentDelay(jj,:)];
            subplot(1,5,ii)
            hold on
            plot(u(ii,:),percentDelay(jj,:),'o-')
        end
    end
    these(ii,:) = mean(this,1,'omitnan');
end
figure
for ii = 1:5
    subplot(2,3,ii)
    plot(u(ii,:),these(ii,:),'.-')
    set(gca,'xtick',u(ii,:))
    xlim([0 1])
end
%% Plot trapz AUC from just last 3 baseline sessions of each
load('F:\irdmRound2\lastBaseNorm2.mat','delay','percentDelay')
[allDelayAnimal,allPercentAnimal] = deal(cell(1,1));
c = 1;
for ii = 1:size(percentDelay,1)
    if ~isempty(percentDelay{ii,1})
        allDelayAnimal{c} = cat(1,delay{ii,:});
        allPercentAnimal{c} = cat(1,percentDelay{ii,:});
        c = c+1;
    end
end
allDelay = cat(1,allDelayAnimal{:});
allPercent = cat(1,allPercentAnimal{:});
u = unique(allDelay,'rows');
figure
for ii = 1:size(u,1)
    this = [];
    for jj = 1:size(allDelay,1)
        if sum(allDelay(jj,:) == u(ii,:)) == 5
            this = [this;allPercent(jj,:)];
            subplot(1,5,ii)
            hold on
            plot(u(ii,:),allPercent(jj,:),'o-')
        end
    end
    these(ii,:) = mean(this,1,'omitnan');
end
figure
for ii = 1:size(u,1)
    subplot(2,2,ii)
    plot(u(ii,:),these(ii,:),'.-')
    set(gca,'xtick',u(ii,:))
    ylim([0 1])
end
%% Get list of aucs and filenames from processed files
files = fileSearch('F:\irdmRound2\processedContinuous3\','.mat')';
% processedData = cell(1,6);
for ii = 874:numel(files)
    if any(ii == 1:50:numel(files))
        disp(num2str(ii))
    end
    parts = strsplit(files{ii},'_');
    ID = strsplit(parts{1},'-');
    % Find index that matches ID and date
    ind = logicFind(1,contains(thisData(:,1),ID{1}) ...
        & contains(thisData(:,2),parts{end-1}),'==');
    % Check against spreadsheet data - if not in spreadsheet, insert NaN
    if isempty(ind)
        processedData(ii,4) = {NaN};
        processedData(ii,6) = {'not in spreadsheet'};
    else
        processedData{ii,4} = round(thisData{ind,3},4);
        % Check if 'bad events' 'missing events', if so skip running
        % ddtTrials
        if contains(thisData{ind,4},'bad') || ...
                contains(thisData{ind,4},'missing event') || ...
                contains(thisData{ind,4},'seizure') || ...
                contains(thisData{ind,4},'sedated')
            trials = [];
        else
            load(files{ii},'trials')
%             if isempty(trials)
%                 % Re-run ddtTrials to account for fixed delay time rounding
%                 load(files{ii},'LFPTs','hist')
%                 trials = ddtTrials(hist.eventTs,LFPTs,0);
%                 keyboard
%                 save(files{ii},'-append','trials')
%             end
        end
        % Insert notes
        processedData(ii,6) = thisData{ind,4};
        % Swap lever event markers for file with 'flipped lever' but not
        % 'skip'
%         if contains(thisData{ind,4},'flipped') && ...
%                 ~contains(thisData{ind,4},'skip') && ...
%                 processedData{ii,4} ~= round(trials(1).auc,4)
%             load(files{ii},'LFPTs','hist')
%             old = hist.eventTs;
%             hist.eventTs.t{1} = old.t{2};
%             hist.eventTs.t{2} = old.t{1};
%             trials = ddtTrials(hist.eventTs,LFPTs,0);
%             if round(trials(1).auc,3) ~= processedData{ii,4}
%                 keyboard
%             else
%                 % Save over trials
%                 save(files{ii},'-append','trials')
%             end
%         end
    end
    % If trials is empty, check for manual trapz or use NaN
    thisInd = logicFind(1,contains(thisFile(:,1),ID{1}) & ...
        contains(thisFile(:,2),parts{end-1}),'==');
    if isempty(trials)
        if ii <= numel(trapzAUCs)
            if ~isempty(thisInd)
                trapzAUC = trapzAUCs(thisInd);
            end
        else
            auc = NaN;
            trapzAUC = NaN;
        end
    else
        auc = round(trials(1).auc,4);
        trapzAUC = trials(1).trapzAUC;
    end
    % Fill in processedData with ID, date, and auc
    processedData(ii,1:3) = {ID{1},parts{end-1},auc};
    % Count clean trials
    load(files{ii},'psdTrls','hist')
    % Check if stim, if so, account for period before and after stim; also
    % remove data with stim events
    if contains(processedData{ii,6},'sIL') || ...
            contains(processedData{ii,6},'sNAcC')
        if isempty(hist.eventTs.t{9}) || numel(hist.eventTs.t{9})<3
            warning('this files has no stim events')
            load(files{ii},'LFPTs')
            figure
            plot(LFPTs.tvec,LFPTs.data(1,:))
            keyboard
        else
            allStim = sort([hist.eventTs.t{9};hist.eventTs.t{10}],'ascend');
            start = allStim(1);
            stop = allStim(end)+60;
        end
        stimInds = psdTrls{1}.t+1.5 > start & psdTrls{1}.t+1.5 < stop;
        % NaN data that overlaps with stim events First find when stim
        % turns off (longer than 3 second delay between stim pulses)
        offs = [logicFind(1,diff(hist.eventTs.t{9})>3,'=='),...
            numel(hist.eventTs.t{9})];
        % Add one (and first stim) to get on indices
        ons = [1,offs(1:end-1)+1];
        % Convert to time
        offs = hist.eventTs.t{9}(offs);
        ons = hist.eventTs.t{9}(ons);
        % Go through time vector of psdTrls and find any times that fall
        % within a stim interval
        for jj = 1:numel(ons)
            rmv = psdTrls{1}.t-1.5>ons(1) & psdTrls{1}.t+1.5<offs(1);
            psdTrls{1}.Pow(:,:,:,rmv) = NaN;
            psdTrls{1}.Overall(:,:,:,rmv) = NaN;
            psdTrls{1}.OverallStd(:,:,:,rmv) = NaN;
            psdTrls{1}.bandPow(:,:,:,rmv) = NaN;
            psdTrls{1}.totPow(:,:,:,rmv) = NaN;
            psdTrls{1}.relPow(:,:,:,rmv) = NaN;
            psdTrls{1}.avgRelPow(:,:,:,rmv) = NaN;
            psdTrls{1}.stdRelPow(:,:,:,rmv) = NaN;
            
            coh{1}.Cxy(:,:,:,rmv) = NaN;
            coh{1}.mtCxy(:,:,:,rmv) = NaN;
            coh{1}.mBandCoh(:,:,:,rmv) = NaN;
            coh{1}.normBandCoh(:,:,:,rmv) = NaN;
        end
        % N.B. don't save over. the rmv for coh doesn't really work...
        % Save
%         save(files{ii},'psdTrls','coh','-append')
        processedData{ii,7} = sum(~isnan(psdTrls{1}.avgRelPow(1,1,1,...
            stimInds)));

    else
        processedData{ii,7} = sum(~isnan(psdTrls{1}.avgRelPow(1,1,1,:)));
    end
    trials = [];
    processedData{ii,8} = trapzAUC;
    % Fill in delays and percent choices
    processedData{ii,9} = thisDelay{thisInd};
    processedData{ii,10} = percentDelay(thisInd,:);
end
processedData(:,5) = num2cell(cell2mat(processedData(:,3)) == ...
    cell2mat(processedData(:,4)));
mismatch = logicFind(1,~contains(processedData(:,6),{'missing','jam',...
    'seizure','bad','skip','mismatch','spreadsheet','task'}) & ...
    cell2mat(processedData(:,3)) ~= cell2mat(processedData(:,4)),'==');
%save('F:\irdmRound2\ddtFiles.mat','processedData','thisData','mismatch','missing')
%% Plot example animal's trapzAUC over time with interventions, etc
animalID = 'IRDM28';
inds = contains(processedData(:,1),animalID);
these = cell2mat(processedData(inds,8));
theseDay = cell2mat(processedData(inds,2));
these(these==0) = NaN;
[all,allI] = sort(datetime(theseDay));
first = all(1);
nDays = days(all-first);

figure
plot(nDays,these(allI),'.k')
set(gca,'xtick',nDays,'xticklabel',theseDay(allI,:))
% set(gca,'xtick',nDays,'xticklabel',nDays)
xtickangle(45)
box off
ylabel('AUC')
xlabel('date')
%% Plot interventions over time from baseline (z-score)
load('F:\irdmRound2\interData.mat', 'interTable')
%%
c = 1;
inter = {'sIL','sNAcC','mph3','sIL-LSD','sNAcC-LSD'};
% Get each animal ID
uAnimal = unique(interTable.ID);
for ii = 1:numel(uAnimal)
    % Get interventions
    for jj = 1:numel(inter)
        theseInter = find(strcmp(interTable.Rx,inter{jj}) & ...
            strcmp(interTable.ID,uAnimal{ii}));
        if ~isempty(theseInter)
           uAUC = unique(interTable.b_aucM(theseInter));
           if numel(uAUC) > 1
               for aC = 1:numel(uAUC)
                   these = theseInter(interTable.b_aucM(theseInter)==uAUC(aC));
                   thisInterAUCZ{c,jj} = (interTable.auc(these)-...
                    interTable.b_aucM(these))./interTable.b_aucS(these);
                thisInterDays{c,jj} = days(datetime(cell2mat(interTable.Date(...
                    these)))-datetime(cell2mat(interTable.Date(...
                    these(1)))));
               end
           else
                             
                thisInterAUCZ{c,jj} = (interTable.auc(theseInter)-...
                    interTable.b_aucM(theseInter))./interTable.b_aucS(theseInter);
                thisInterDays{c,jj} = days(datetime(cell2mat(interTable.Date(...
                    theseInter)))-datetime(cell2mat(interTable.Date(...
                    theseInter(1)))));
                %             figure
                %             scatter(thisInterDays{c,jj},thisInterAUCZ{c,jj})
                %             title([uAnimal{ii},inter{jj}])
                %             keyboard
                c = c+1;
           end
        end
    end
end
%%
col = {'r','b','k','m','g'};
figure
for ii = 1:numel(inter)
    subplot(3,2,ii)
    theseInd = logicFind(1,~cellfun(@isempty,thisInterAUCZ(:,ii)),'==');
    hold on
    for jj = theseInd
        plot(thisInterDays{jj,ii},thisInterAUCZ{jj,ii},'-o','color',col{ii})
    end
end
%% unpack interAUCZ

%% Use to check old trials and replace
for ii = 1:size(processedData,1)
    if processedData{ii,5} == 0 && ~contains(processedData{ii,6},'bad') ...
            && ~contains(processedData{ii,6},'spread')
        load(['E:\processedMat\',files{ii}(1:end-8),'.mat'],...
            'trials','eventTs')
        hist.eventTs = eventTs;
        trials = ddtTrials(hist.eventTs,LFPTs,0);
        if round(trials(1).auc,3) ~= processedData{ii,4}
            keyboard
        else
            % Save over trials
            save(files{ii},'-append','trials')
        end
    end
end
%% Flip levers
for ii = [392:397,768,769]
    load(files{ii},'LFPTs','hist')
    old = hist.eventTs;
    hist.eventTs.t{1} = old.t{2};
    hist.eventTs.t{2} = old.t{1};
    trials = ddtTrials(hist.eventTs,LFPTs,0);
    if round(trials(1).auc,3) ~= processedData{ii,4}
        keyboard
    else
        % Save over trials
        save(files{ii},'-append','trials')
    end
end
