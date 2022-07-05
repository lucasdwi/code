load('F:\irdmRound2\ddtFiles2.mat')
% Remove files with: 'skip','neural','seizure'
skipInd = contains(processedData(:,6),'skip') | ...
    contains(processedData(:,6),'neural') | ...
    contains(processedData(:,6),'seizure') | ...
    contains(processedData(:,6),'OFC');
processedData(skipInd,:) = [];

% Get list of all processed files
files = fileSearch('F:/irdmRound2/processedContinuous3/','.mat');
% Set n, number of baselines to use
n = 5;
% Set up indices of pow of coh
powInds = reshape(1:48,6,8)';
cohInds = reshape(49:216,6,28)';
% Get unique IDs
uID = unique(processedData(:,1));
% Go through each animal and find the recordings with interventions and
% their prior baselines; also calculate number of clean trials in data
% inter = {'sIL','sNAcC','guan0.3','mph3','mph0.3','mph1','mph0.1',...
%     'sNAcC-LSD','sIL-LSD'};
% inter = {'mph3'};
inter = {'sIL','sNAcC','mph3'};
% Clear table
interTable = array2table(zeros(0,11));
interBaseInd = cell(numel(uID),numel(inter));
% c = 72;
for ii = 1:numel(uID)
    disp(ii)
    % Get logical for animal ID matches
    idLog = contains(processedData(:,1),uID{ii});
    idInd = logicFind(1,idLog,'==');
    % Sort by date within each animal
    [theseDates,sortInd] = sort(datetime(processedData(idLog,2)));
    processedData(idLog,:) = processedData(idInd(sortInd),:);
    % Base inds
    baseInd{ii} = logicFind(1,contains(processedData(:,6),...
        'base') & idLog,'==');
    % Find recordings from each animal (ii) with each intervention (jj)
    for jj = 1:numel(inter)
        % For first two interventions, omit entries with LSD
        if jj <= 2
            intInd{ii,jj} = logicFind(1,contains(processedData(:,6),...
                inter{jj}) & idLog & ...
                ~contains(processedData(:,6),'LSD'),'==');
        else
            intInd{ii,jj} = logicFind(1,contains(processedData(:,6),...
                inter{jj}) & idLog,'==');
        end
    end
    % Figure out which is the first by using indices
    % (sorted by date)
    firstInt = min(cell2mat(cellfun(@min,intInd(ii,:),...
        'UniformOutput',false)));
    % For each intervention, find the n closest baselines before and all
    % baselines after without an effect
    for jj = 1:numel(inter)
        for k = 1:numel(intInd{ii,jj})
            % Grab index of file in processedData
            thisInd = intInd{ii,jj}(k);
            if processedData{thisInd,7} ~= 0
                % Get this date
                thisDate = datetime(processedData(thisInd,2));
                % For loading check if stim
                newx = {'stimIL','stimCore'};
                if any(jj == [1,2])
                    thisInter = newx{jj};
                else
                    thisInter = inter{jj};
                end
                % Find all baselines that are prior
                theseBaseInd = baseInd{ii}(datetime(...
                    processedData(baseInd{ii},2)) < thisDate);
                % Count number days between this intervention and the
                % baselines
                nDays = days(thisDate - ...
                    datetime(processedData(theseBaseInd(end-(n-1):end),2)));
                % Store info, including samples and nDays
                interBase{ii,jj,k} = [processedData(theseBaseInd(end-(n-1):end),:),...
                    num2cell(nDays)];
                interBaseInd{jj,k} = num2str(theseBaseInd(end-(n-1):end));              
            end
        end
    end
    % Load all unique sets of baselines
    uBSets = unique(interBaseInd(~cell2mat(cellfun(@isempty,...
        interBaseInd,'UniformOutput',0))));
    for jj = 1:size(uBSets,1)
        these = str2num(uBSets{jj});
        theseBase = [];
        for k = 1:n
            % Find file using animal and date
            thisFile = files(contains(files,...
                processedData{these(k),1}) & ...
                contains(files,processedData{these(k),2}));
            % Load and remove NaNed data
            load(['F:/irdmRound2/processedContinuous3/',thisFile{1}]);
            thisNotNan = ~isnan(squeeze(psdTrls{1}.bandPow(1,1,1,:))) &...
                ~isnan(squeeze(coh{1}.mBandCoh(1,1,1,:)));
            notNanPow = squeeze(psdTrls{1}.bandPow(:,:,1,thisNotNan));
            notNanCoh = squeeze(coh{1}.mBandCoh(:,:,1,thisNotNan));
            % Fill in NaNs from missing channels if any
            if size(LFPTs.data,1) < 8
                [newNotNanPow,newNotNanCoh] = addNaNs(LFPTs,...
                    notNanPow,notNanCoh,8);
            else
                newNotNanPow = notNanPow;
                newNotNanCoh = notNanCoh;
            end
            % Reshape
            rePow = reshape(newNotNanPow,48,sum(thisNotNan))';
            reCoh = reshape(permute(newNotNanCoh,[2,1,3]),168,...
                sum(thisNotNan))';
            thisBase = [rePow,reCoh];
            % Combine
            theseBase = [theseBase;thisBase];
        end
        uBSets{jj,2} = theseBase;
        uBSets{jj,3} = mean(cell2mat(processedData(these,8)));
        uBSets{jj,4} = std(cell2mat(processedData(these,8)));
    end
    % Then, load each intervention file and normalize by proper baseline
    for jj = 1:numel(inter)
        for k = 1:numel(intInd{ii,jj})
            thisInd = intInd{ii,jj}(k);
            % If the first intervention, set aside baseline data
            if thisInd == firstInt
                lastBase{ii,jj} = uBSets{jj,2};
            end
            % Load file
            thisFile = files(contains(files,processedData{thisInd,1}) & ...
                contains(files,processedData{thisInd,2}));
            load(['F:/irdmRound2/processedContinuous3/',thisFile{1}]);
            if contains(inter{jj},'sIL') || contains(inter{jj},'sNAc')
                % Get first stim time
                first = min([hist.eventTs.t{9}(1),hist.eventTs.t{10}(1)]);
                last = max([hist.eventTs.t{9}(end),...
                    hist.eventTs.t{10}(end)]);
                % Grab any data that comes from between stim (during actual
                % stim should already be removed)
                notNaN = squeeze(~isnan(psdTrls{1}.bandPow(1,1,1,:)))';
                % Use -1.5 to account for window, and +60 to account for
                % washout, and notNaN
                interInd = psdTrls{1}.t > first-1.5 & ...
                    psdTrls{1}.t < last+60 & notNaN;
            else
                % Othewise just take out the NaNs
                interInd = ~isnan(squeeze(psdTrls{1}.bandPow(1,1,1,:)))'...
                    & ~isnan(squeeze(coh{1}.mBandCoh(1,1,1,:)))';
            end
             % Fill missing data in with NaNs
             stimPow = squeeze(psdTrls{1}.bandPow(:,:,1,interInd));
             stimCoh = squeeze(coh{1}.mBandCoh(:,:,1,interInd));
             % Check if empty
             if ~isempty(stimPow)
                 if size(LFPTs.data,1) < 8
                     [newPow,newCoh] = addNaNs(LFPTs,...
                         stimPow,stimCoh,8);
                 else
                     newPow = stimPow;
                     newCoh = stimCoh;
                 end
                 % Reshape
                 rePow = reshape(newPow,48,sum(interInd))';
                 reCoh = reshape(permute(newCoh,[2,1,3]),168,sum(interInd))';
                 thisStim = [rePow,reCoh];
                 % Grab baseline and get pop stats
                 uBSetInd = contains(uBSets(:,1),interBaseInd{jj,k});
                 thisBase = uBSets{uBSetInd,2};
                 normData = (thisStim-mean(thisBase,1))./std(thisBase,[],1);
                 base{ii,jj,k} = thisBase;
                 stim{ii,jj,k} = thisStim;
             else
                 normData = NaN;
             end
             % Determine effect based on mean ± 2 std
             if processedData{thisInd,8} >= uBSets{uBSetInd,3}+...
                     2*uBSets{uBSetInd,4}
                 effect{ii,jj}(k) = 1;
             end
             if processedData{thisInd,8} < uBSets{uBSetInd,3}-...
                     2*uBSets{uBSetInd,4}
                 effect{ii,jj}(k) = -1;
             end
             if processedData{thisInd,8} >= uBSets{uBSetInd,3}-...
                     2*uBSets{uBSetInd,4} && ...
                     processedData{thisInd,8} <  uBSets{uBSetInd,3}+...
                     2*uBSets{uBSetInd,4}
                 effect{ii,jj}(k) = 0;
             end
             % Get sum of base samples
             baseSamp = sum(cell2mat(processedData(...
                 theseBaseInd(end-2:end),7)));
             % Get sum of base samples if using smallest number
             baseMinSamp = min(cell2mat(processedData(...
                 theseBaseInd(end-2:end),7)))*3;
             % Put in table
             thisPosX = cell2table([uID(ii),processedData(thisInd,2),...
                 processedData(thisInd,6),processedData(thisInd,7),...
                 baseSamp,baseMinSamp,processedData(thisInd,8),...
                 uBSets{uBSetInd,3},effect{ii,jj}(k),normData,...
                 uBSets{uBSetInd,4}]);
             interTable = [interTable;thisPosX];
%             interTable(c,:) = thisPosX;
%             c = c+1;
        end
    end
    %
end
%
interTable.Properties.VariableNames = {'ID','Date','Rx','samps',...
    'max_b_samp','min_b_samp','auc','b_aucM','effect','normData','b_aucS'};
% Get number of recordings per intervention
nInt = cellfun(@numel,intInd);
% Get number of samples per intervention
nSamp = cellfun(@(x) sum(cell2mat(processedData(x,7))),intInd);
% Get number of samples per intervention per recording
nSampInt = cellfun(@(x) arrayfun(@(y) cell2mat(processedData(y,7)),x),...
    intInd,'UniformOutput',false);
sampArray = cell(1,9);
for ii = 1:size(effect,1)
    for jj = 1:size(effect,2)
        effectCount{ii,jj} = [sum(effect{ii,jj}==1),...
            sum(effect{ii,jj}==-1),sum(effect{ii,jj}==0)];
        effectSample{ii,jj} = [sum(nSampInt{ii,jj}(effect{ii,jj}==1)),...
            sum(nSampInt{ii,jj}(effect{ii,jj}==-1)),...
            sum(nSampInt{ii,jj}(effect{ii,jj}==0))];
        if any(effectSample{ii,jj} ~= [0,0,0])
            sampArray{jj} = [sampArray{jj};effectSample{ii,jj}];
        end
    end
end
% save('F:/irdmRound2/interData.mat','interTable','effect','effectCount','effectSample')
% save('F:/irdmRound2/lastBase.mat','lastBase')
save('F:/irdmRound2/stimBase.mat','stim','base','inter')
%% Duration of behavioral effect
% Set number of days of baseline to use
n = 5;
% Go through each animal
uID = unique(processedData(:,1));
inter = {'sIL';'sNAcC'};
for ii = 1:numel(uID)
    disp(ii)
    % Get logical for animal ID matches
    idLog = contains(processedData(:,1),uID{ii});
    idInd = logicFind(1,idLog,'==');
    % Sort by date within each animal
    [theseDates,sortInd] = sort(datetime(processedData(idLog,2)));
    processedData(idLog,:) = processedData(idInd(sortInd),:);
    % Base inds
    baseInd{ii} = logicFind(1,contains(processedData(:,6),...
        'base') & idLog,'==');
    % Find recordings from each animal (ii) with each intervention (jj)
    for jj = 1:numel(inter)
        % For first two interventions, omit entries with LSD
        if jj <= 2
            intInd{ii,jj} = logicFind(1,contains(processedData(:,6),...
                inter{jj}) & idLog & ...
                ~contains(processedData(:,6),'LSD'),'==');
        else
            intInd{ii,jj} = logicFind(1,contains(processedData(:,6),...
                inter{jj}) & idLog,'==');
        end
    end
    % Find n recordings before intervention to get baseline
    for jj = 1:numel(inter)
        for k = 1:numel(intInd{ii,jj})
            % Grab index of file in processedData
            thisInd = intInd{ii,jj}(k);
            if processedData{thisInd,7} ~= 0
                % Get this date
                thisDate = datetime(processedData(thisInd,2));
                % For loading check if stim
                newx = {'stimIL','stimCore'};
                if any(jj == [1,2])
                    thisInter = newx{jj};
                else
                    thisInter = inter{jj};
                end
                % Find all baselines that are prior
                theseBaseInd = baseInd{ii}(datetime(...
                    processedData(baseInd{ii},2)) < thisDate);
                % Count number days between this intervention and the
                % baselines
                nDays = days(thisDate - ...
                    datetime(processedData(theseBaseInd(end-(n-1):end),2)));
                % Store info, including samples and nDays
                interBase{ii,jj,k} = [processedData(theseBaseInd(end-(n-1):end),:),...
                    num2cell(nDays)];
                interBaseInd{ii,jj,k} = num2str(theseBaseInd(end-(n-1):end));
                interBaseAUC{ii,jj,k} = processedData(...
                    theseBaseInd(end-(n-1):end),8);
            end
        end
    end
    % Find recordings after an intervention and check for effect
    for jj = 1:numel(inter)
        if ~isempty(intInd{ii,jj})
        % Find breaks in intervention
        breakInd = logicFind(1,diff(intInd{ii,jj}),'~=');
        if isempty(breakInd)
            lastInd = intInd{ii,jj}(end);
        else
            lastInd = [intInd{ii,jj}(breakInd),intInd{ii,jj}(end)];
        end
                               
        for k = 1:numel(lastInd)
            % Get this data
            thisDate = datetime(processedData(lastInd(k),2));
            % Find all baselines that are prior
            theseBaseInd = baseInd{ii}(datetime(...
                processedData(baseInd{ii},2)) < thisDate);
            % Just grab n baselines and count number days between this
            % intervention and the baselines
            nDays = days(thisDate - ...
                datetime(processedData(theseBaseInd(end-(n-1):end),2)));
            baseAUC = cell2mat(processedData(theseBaseInd(end-(n-1):end)...
                ,8));
            
            % Find next baselines
            postBaseInd = baseInd{ii}(datetime(processedData(baseInd{ii}...
                ,2)) > datetime(processedData(lastInd(k),2)));
            postBaseDays{ii,jj,k} = days(datetime(processedData(postBaseInd,2)) - thisDate);
            postAUC = cell2mat(processedData(postBaseInd,8));
            effect{ii,jj,k} = postAUC >= mean(baseAUC) + 2*std(baseAUC) | ...
                postAUC <= mean(baseAUC) - 2*std(baseAUC);
            if ~isempty(effect{ii,jj,k}) && effect{ii,jj,k}(1) == 1
                disp('effect')
                disp(effect{ii,jj,k})
                disp('days')
                disp(postBaseDays{ii,jj,k})
                dummy = 1;
            end
            % Figure out if last intervention day was successful
            lastEffect{ii,jj,k} = cell2mat(processedData(lastInd(k),8)) ...
                >= mean(baseAUC) + 2*std(baseAUC) | ...
                cell2mat(processedData(lastInd(k),8)) <= ...
                mean(baseAUC) - 2*std(baseAUC);
        end
        end
    end
end
%% Cheat hard code
% IRDM35 sIL had seizure right after; remove
effect{16,1,1} = [];
postBaseDays{16,1,1} = [];
% IRDM 35 sNAcC waited too long after stim
effect{16,2,2} = [];
postBaseDays{16,2,2} = [];
% IRDM38 had an early sNAcC
effect{19,2,1} = [];
postBaseDays{19,2,1} = [];
% IRDM41 had bad head detector during stim
effect{28,2,1} = [];
postBaseDays{28,2,1} = [];
%%
dur = cell(1,2);
for ii = 1:size(effect,1)
    for jj = 1:size(effect,2)
        for k = 1:size(effect,3)
            if ~isempty(effect{ii,jj,k})
                if effect{ii,jj,k}(1) == 1
                    dur{jj} = [dur{jj};postBaseDays{ii,jj,k}(1)];
                    if numel(effect{ii,jj,k}) > 1
                        c = 2;
                        while c < numel(effect{ii,jj,k}) & ...
                                effect{ii,jj,k}(c) == 1
                            dur{jj} = [dur{jj};postBaseDays{ii,jj,k}(c)];
                            c = c+1;
                        end
                    end
                end
            end
        end
    end
end
%%
c=1;
dur = cell(1,2);
for jj = 1:size(effect,2)
    c=1;
    for ii = 1:size(effect,1)
        for k = 1:size(effect,3)
            if ~isempty(effect{ii,jj,k})
                n = 1;
                this = effect{ii,jj,k}(n);
                while this == 1 && n+1<numel(effect{ii,jj,k})
                    n = n+1;
                    this = effect{ii,jj,k}(n);
                end
                postBaseDays{ii,jj,k}(n);
                dur{jj}(c) = postBaseDays{ii,jj,k}(n);
%                 dur{jj}(c) = n-1;
                last{jj}(c) = lastEffect{ii,jj,k};
                c = c+1;
            end
        end
    end
end
%%
for ii = 1:2
    c = 1;
    for jj = unique(dur{ii})
        perc{ii}(c) = sum(dur{ii}==jj)/numel(dur{ii});
        c = c+1;
    end
end
figure
hold on
plot(0:2,perc{1},'-o')
plot(0:2,[perc{2},0],'-o')
set(gca,'xtick',0:1:2,'ytick',0:.2:1)
legend({'sIL','sNAcC'})
xlabel('consecutive days with effect')
ylabel('% of interventions')
%%
bar([perc{1};perc{2},0],'stacked')
set(gca,'xticklabel',{'sIL','sNAcC'})
legend({['0 days: ',num2str(round(perc{1}(1),2)),':',...
    num2str(round(perc{2}(1),2))],['1 day: ',...
    num2str(round(perc{1}(2),2)),':',num2str(round(perc{2}(2),2))],...
    ['2 days: ',num2str(round(perc{1}(3),2)),':0']})
%% normalized features figure
inter = {'sIL','sNAcC','mph3'};
for ii = 1:3
    figure
    hold on
    % Plot zeros
    these = logicFind(1,strcmp(inter{ii},interTable.Rx) & interTable.effect == 0,'==');
    allZero = [];
    for jj = these
        allZero = [allZero;interTable.normData{jj}(:,48:6:end)];
    end
    % Plot pos
    these = logicFind(1,strcmp(inter{ii},interTable.Rx) & interTable.effect == 1,'==');
    allPos = [];
    for jj = these
        allPos = [allPos;interTable.normData{jj}(:,48:6:end)];
    end
    % Plot neg
    these = logicFind(1,strcmp(inter{ii},interTable.Rx) & interTable.effect == -1,'==');
    allNeg = [];
    for jj = these
        allNeg = [allNeg;interTable.normData{jj}(:,48:6:end)];
    end
    figure
    hold on
    for jj = 1:28
%         [~,p(jj)] = ttest2(allZero(:,jj),allPos(:,jj));
%         violin({allZero(:,jj),allPos(:,jj)},'x',[(jj-1)*2+1,jj*2])
        plot([1,2],[mean(allZero(:,jj)),mean(allPos(:,jj))])
    end
end
%% Z-score AUC figure
inter = {'sIL','sNAcC','mph3'};
for ii = 1:numel(inter)
    z{ii} = (interTable.auc(contains(interTable.Rx,inter{ii}))-interTable.b_aucM(contains(interTable.Rx,inter{ii})))./interTable.b_aucS(contains(interTable.Rx,inter{ii}));
end

figure
subplot(1,2,1)
hold on
for ii = 1:3
    plot(ones(1,sum(abs(z{ii})<2))*ii,abs(z{ii}(abs(z{ii})<2)),'.k')
    plot(ones(1,sum(abs(z{ii})>=2))*ii,abs(z{ii}(abs(z{ii})>=2)),'.b');
end
set(gca,'xtick',1:3,'xticklabel',inter)
xlim([0.5 3.5])
ylabel('zscore')

subplot(1,2,2)
hold on
for ii = 1:3
    plot(ones(1,numel(z{ii}))*ii,z{ii},'.k')
end
set(gca,'xtick',1:3,'xticklabel',inter)
violin(z,'facealpha',0);
xlim([0.5 3.5])
ylabel('zscore')
%% z-score distributions
figure
hold on
histogram(z{1},'binwidth',2,'normalization','probability',...
    'displaystyle','stairs')
histogram(z{2},'binwidth',2,'normalization','probability',...
    'displaystyle','stairs')
histogram(z{3},'binwidth',2,'normalization','probability',...
    'displaystyle','stairs')
legend(inter)
xlabel('z-score')
ylabel('proportion')
%% Base vs. Stim (i.e., stim efffects)
% cd('F:/irdmRound2/baseVstim/')
for ii = 1:20
%     load(['baseStim',num2str(ii),'.mat'])
    for jj = 1:sum(~cellfun(@(x) isempty(x),acc(1,:)))
        ilA(ii,jj) = acc{1,jj}{1}.acc;
        ilX(ii,jj,:) = acc{1,jj}{1}.x;
        ilY(ii,jj,:) = acc{1,jj}{1}.y;
        ilSA(ii,jj,:) = a{1,jj};
        ilCoeff(ii,jj,:) = coeff{1,jj};
    end
    for jj = 1:sum(~cellfun(@(x) isempty(x),acc(2,:)))
        coreA(ii,jj) = acc{2,jj}{1}.acc;
        coreX(ii,jj,:) = interp1(linspace(0,1,numel(acc{2,jj}{1}.x)),...
            acc{2,jj}{1}.x,linspace(0,1,200));
        coreY(ii,jj,:) = interp1(linspace(0,1,numel(acc{2,jj}{1}.y)),...
            acc{2,jj}{1}.y,linspace(0,1,200));
        coreSA(ii,jj,:) = a{2,jj};
        coreCoeff(ii,jj,:) = coeff{2,jj};
    end
    for jj = 1:sum(~cellfun(@(x) isempty(x),acc(3,:)))
        mphA(ii,jj) = acc{3,jj}{1}.acc;
        mphX(ii,jj,:) = interp1(linspace(0,1,numel(acc{3,jj}{1}.x)),...
            acc{3,jj}{1}.x,linspace(0,1,200));
        mphY(ii,jj,:) = interp1(linspace(0,1,numel(acc{3,jj}{1}.y)),...
            acc{3,jj}{1}.y,linspace(0,1,200));
        mphSA(ii,jj,:) = a{3,jj};
        mphCoeff(ii,jj,:) = coeff{3,jj};
    end
end
%%
figure
hold on
plot(squeeze(mean(ilX,[1,2])),squeeze(mean(ilY,[1,2])))
plot(squeeze(mean(coreX,[1,2])),squeeze(mean(coreY,[1,2])))
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
legend({['il: ',num2str(round(mean(ilA,[1,2]),2)),'\pm',...
    num2str(round(conf(mean(ilA,1),0.95),2))],['core: ',...
    num2str(round(mean(coreA,[1,2]),2)),'\pm',...
    num2str(round(conf(mean(coreA,1),0.95),2))]})
%%
feat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
mILSA = squeeze(mean(ilSA,[1,2]));
mILC = squeeze(mean(ilCoeff,[1,2]));
ilSign = sign(mILC);
mCoreSA = squeeze(mean(coreSA,[1,2]));
mCoreC = squeeze(mean(coreCoeff,[1,2]));
coreSign = sign(mCoreC);

c = distinguishable_colors(6);
figure
hold on
scatter(mILSA(1:48).*ilSign(1:48),mCoreSA(1:48).*coreSign(1:48),[],repmat(c,8,1),'o');
scatter(mILSA(49:end).*ilSign(49:end),mCoreSA(49:end).*coreSign(49:end),[],repmat(c,28,1),'s');
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
xlabel('IL'); ylabel('Core')

[~,ind] = sort(mILSA,'descend');
ilSort = mILSA(ind).*ilSign(ind);
ilFeat = feat(ind)';
[~,ind] = sort(mCoreSA,'descend');
coreSort = mCoreSA(ind).*coreSign(ind);
coreFeat = feat(ind)';
%% 
figure
subplot(1,2,1)
thisPow = padarray(reshape(mILSA(1:48),6,8).*reshape(ilSign(1:48),6,8),...
    [1 1],'post');
thisPow(thisPow<0.6 & thisPow>-0.6) = NaN;
thisPow(thisPow>0.6) = 1;
thisPow(thisPow<-0.6) = -1;
pcolor(thisPow)
set(gca,'xtick',1.5:1:8.5,'xticklabel',...
    {'lPFC','rPFC','lOFC','rOFC','lNAcs','rNAcS','lNAcC','rNAcC'},...
    'ytick',1.5:1:7.5,'yticklabel',...
    {'delta','theta','alpha','beta','lgamma','hgamma'})
colormap viridis
caxis([-0.6 1])
title('IL stim vs. base: pow')
xtickangle(45)
subplot(1,2,2)
thisPow = padarray(reshape(mCoreSA(1:48),6,8).*reshape(coreSign(1:48),6,8),...
    [1 1],'post');
thisPow(thisPow<0.6 & thisPow>-0.6) = NaN;
thisPow(thisPow>0.6) = 1;
thisPow(thisPow<-0.6) = -1;
pcolor(thisPow)
set(gca,'xtick',1.5:1:8.5,'xticklabel',...
    {'lPFC','rPFC','lOFC','rOFC','lNAcs','rNAcS','lNAcC','rNAcC'},...
    'ytick',1.5:1:6.5,'yticklabel',...
    {'delta','theta','alpha','beta','lgamma','hgamma'})
colormap viridis
caxis([-0.6 1])
title('core stim vs. base: pow')
xtickangle(45)

pairs = cell(1,28);
sites = {'lPFC','rPFC','lOFC','rOFC','lNAcs','rNAcS','lNAcC','rNAcC'};
cmbs = nchoosek(1:8,2);
for ii = 1:size(cmbs,1)
    pairs{ii} = [sites{cmbs(ii,1)},'-',sites{cmbs(ii,2)}];
end
figure
subplot(1,2,1)
thisCoh = padarray(reshape(mILSA(49:end),6,28).*reshape(ilSign(49:end),6,28),...
    [1 1],'post');
thisCoh(thisCoh<0.6 & thisCoh>-0.6) = NaN;
thisCoh(thisCoh>0.6) = 1;
thisCoh(thisCoh<-0.6) = -1;
pcolor(thisCoh)
set(gca,'xtick',1.5:1:28.5,'xticklabel',...
   pairs,'ytick',1.5:1:6.5,'yticklabel',...
    {'delta','theta','alpha','beta','lgamma','hgamma'})
colormap viridis
caxis([-0.7 1])
title('IL stim vs. base: coh')
xtickangle(45)
subplot(1,2,2)
thisCoh = padarray(reshape(mCoreSA(49:end),6,28).*reshape(coreSign(49:end),6,28),...
    [1 1],'post');
thisCoh(thisCoh<0.6 & thisCoh>-0.6) = NaN;
thisCoh(thisCoh>0.6) = 1;
thisCoh(thisCoh<-0.6) = -1;
pcolor(thisCoh)
set(gca,'xtick',1.5:1:28.5,'xticklabel',...
    pairs,'ytick',1.5:1:6.5,'yticklabel',...
    {'delta','theta','alpha','beta','lgamma','hgamma'})
colormap viridis
caxis([-0.7 1])
title('core stim vs. base: coh')
xtickangle(45)
%% Load base -> intervention effect models
[coreX,coreY,coreXL2O,coreYL2O] = deal([]);
for ii = 1:100
    disp(ii)
    load(['F:/irdmRound2/baseEffect/baseEffectMean',...
        num2str(ii),'.mat'])
    coreA(ii) = coreMacc{1}.acc;
    coreL2O(ii) = coreL2OMacc{1}.acc;
    x = coreMacc{1}.x;
    y = coreMacc{1}.y;
    coreX(ii,:) = interp1(linspace(0,1,numel(x)),x,linspace(0,1,7));
    coreY(ii,:) = interp1(linspace(0,1,numel(y)),y,linspace(0,1,7));
    
    corePermA(ii) = coreMPermAcc{1}.acc;
    x = coreMPermAcc{1}.x;
    y = coreMPermAcc{1}.y;
    coreXPerm(ii,:) = interp1(linspace(0,1,numel(x)),x,linspace(0,1,7));
    coreYPerm(ii,:) = interp1(linspace(0,1,numel(y)),y,linspace(0,1,7));
    
    x = coreL2OMacc{1}.x;
    y = coreL2OMacc{1}.y;
    coreXL2O(ii,:) = interp1(linspace(0,1,numel(x)),x,linspace(0,1,11));
    coreYL2O(ii,:) = interp1(linspace(0,1,numel(y)),y,linspace(0,1,11));
    ilA(ii) = multiClassAUC(ilMacc{1}.pred{1},ilMhist.cfg.naive.testY);
    ilAPerm(ii) = multiClassAUC(ilMPermAcc{1}.pred{1},...
        ilMPermHist.cfg.naive.testY);
    ilL3O(ii) = multiClassAUC(ilL3OMacc{1}.pred{1},...
        ilL3OMhist.cfg.naive.testY);
    [~,~,~,ilPosA(ii)] = perfcurve(ilPosHist.cfg.naive.testY,...
        ilPosAcc{1}.pred{1},1);
    load(['F:/irdmRound2/baseEffect/mph3/baseEffectMeanMPH',num2str(ii),...
        '.mat'])
    mphA(ii) = mphMacc{1}.acc;
    x = mphMacc{1}.x;
    y = mphMacc{1}.y;
    mphX(ii,:) = interp1(linspace(0,1,numel(x)),x,linspace(0,1,7));
    mphY(ii,:) = interp1(linspace(0,1,numel(y)),y,linspace(0,1,7));
    
    mphPermA(ii) = mphMPermAcc{1}.acc;
    x = mphMPermAcc{1}.x;
    y = mphMPermAcc{1}.y;
    mphXPerm(ii,:) = interp1(linspace(0,1,numel(x)),x,linspace(0,1,7));
    mphYPerm(ii,:) = interp1(linspace(0,1,numel(y)),y,linspace(0,1,7));
    if ii <= 36
        mphL2O(ii) = mphL2OMacc{1}.acc;
        x = mphL2OMacc{1}.x;
        y = mphL2OMacc{1}.y;
        mphXL2O(ii,:) = interp1(linspace(0,1,numel(x)),x,linspace(0,1,11));
        mphYL2O(ii,:) = interp1(linspace(0,1,numel(y)),y,linspace(0,1,11));
    end
%     [~,~,~,ilPosPermA(ii)] = perfcurve(ilMPosPermHist.cfg.naive.testY,ilMPosPermAcc{1}.pred{1},1);
   % Get single features 
   for jj = 1:216
       % IL
       trainY = ilMhist.trainY;
       trainY(trainY==1)=3;
       trainY(trainY==-1)=2;
       trainY(trainY==0)=1;
       beta = mnrfit(ilMhist.trainX(:,jj),trainY);
       pred = mnrval(beta,ilMhist.cfg.naive.testX(:,jj));
       coeffIL(ii,jj,:) = beta(2,:);
       aIL(ii,jj) = multiClassAUC(pred,ilMhist.cfg.naive.testY);
       
       trainY = ilL3OMhist.trainY;
       trainY(trainY==1)=3;
       trainY(trainY==-1)=2;
       trainY(trainY==0)=1;
       beta = mnrfit(ilL3OMhist.trainX(:,jj),trainY);
       pred = mnrval(beta,ilL3OMhist.cfg.naive.testX(:,jj));
       coeffILL3O(ii,jj,:) = beta(2,:);
       aILL3O(ii,jj) = multiClassAUC(pred,ilL3OMhist.cfg.naive.testY);
       
       % Core
       mdl = fitglm(coreMhist.trainX(:,jj),coreMhist.trainY,...
           'distribution','binomial');
       pred = predict(mdl,coreMhist.cfg.naive.testX(:,jj));
       coeffCore(ii,jj) = table2array(mdl.Coefficients(2,1));
       [~,~,~,aSCore(ii,jj)] = perfcurve(coreMhist.cfg.naive.testY,pred,1);
       
       mdl = fitglm(coreL2OMhist.trainX(:,jj),coreL2OMhist.trainY,...
           'distribution','binomial');
       pred = predict(mdl,coreL2OMhist.cfg.naive.testX(:,jj));
       coeffCoreL2O(ii,jj) = table2array(mdl.Coefficients(2,1));
       [~,~,~,aSCoreL2O(ii,jj)] = perfcurve(coreL2OMhist.cfg.naive.testY,pred,1);
   end
end
save('F:\irdmRound2\baseInterModels.mat','coreX','coreY','coreA',...
    'coeffCore','aSCore','coreXL2O','coreYL2O','coreL2O','coeffCoreL2O',...
    'aSCoreL2O','mphX','mphY','mphA','mphXPerm','mphYPerm','mphXL2O',...
    'mphYL2O','mphPermA','mphL2O','coeffIL','aIL','coeffILL3O','aILL3O')
%% Performance figure
figure
hold on
plot(mean(coreX,1),mean(coreY,1))
plot(mean(coreXL2O,1),mean(coreYL2O,1))
legend({['80/20: \mu = ',num2str(round(mean(coreA),2)),'\pm',...
    num2str(round(conf(coreA,0.95),2))],['L2O: \mu = ',...
    num2str(round(mean(coreL2O),2)),'\pm',...
    num2str(round(conf(coreL2O,0.95),2))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
figure
hold on
plot(mean(mphX,1),mean(mphY,1))
plot(mean(mphXPerm,1),mean(mphYPerm,1))
plot(mean(mphXL2O,1),mean(mphYL2O,1))
legend({['80/20: \mu = ',num2str(round(mean(mphA),2)),'\pm',...
    num2str(round(conf(mphA,0.95),2))],['Perm: \mu = ',...
    num2str(round(mean(mphPermA),2)),'\pm',...
    num2str(round(conf(mphPermA,0.95),2))],['L2O: \mu = ',...
    num2str(round(mean(mphL2O),2)),'\pm',...
    num2str(round(conf(mphL2O,0.95),2))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR')
ylabel('TPR')
title('MPH3')
% doubleHist(coreA,coreL2O)
% doubleHist(ilA,ilAPerm)
%% Single features IL vs. NAcC quadrant plot and feature list
fullFeat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
il = mean(aIL,1).*sign(mean(coeffIL(:,:,2),1));
[~,sortInds] = sort(abs(il),'descend');
sortIL = il(sortInds)';
ilFeat = fullFeat(sortInds)';

ilL3O = mean(aILL3O,1).*sign(mean(coeffILL3O(:,:,2),1));
[~,sortInds] = sort(abs(ilL3O),'descend');
sortILL3O = ilL3O(sortInds)';
ilL3OFeat = fullFeat(sortInds)';

core = mean(aCore,1).*sign(mean(coeffCore,1));
[~,sortInds] = sort(abs(core),'descend');
sortCore = core(sortInds)';
coreFeat = fullFeat(sortInds)';

coreL2O = mean(aCoreL2O,1).*sign(mean(coeffCoreL2O,1));
[~,sortInds] = sort(abs(coreL2O),'descend');
sortCoreL2O = coreL2O(sortInds)';
coreL2OFeat = fullFeat(sortInds)';

c = distinguishable_colors(6);
figure
hold on
scatter(coreL2O(1:48),ilL3O(1:48),[],repmat(c,8,1),'o');
scatter(coreL2O(49:end),ilL3O(49:end),[],repmat(c,28,1),'s');
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
xlabel('core L2O'); ylabel('IL L3O')
%% Single features

figure
plot(il,ilLOO,'ok')


figure
plot(core,coreLOO,'ok')

fullFeat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
[sortIL,sortILind] = sort(mean(aILL3O),'descend');
[sortCore,sortCoreInd] = sort(mean(aCoreL2O),'descend');
ilFeat = fullFeat(sortILind)';
coreFeat = fullFeat(sortCoreInd)';
%% Load effect
load('interData.mat','effect')
% Load last baseline data
load('lastBaseNorm2.mat')
% Remove empty data from both
theseInd = logicFind(1,~cellfun(@isempty,lastBase(:,1)),'==');
thisLastBase = lastBase(theseInd,:);
% Round effect to get single value
thisEffect = round(cellfun(@mean,effect));
thisEffect = thisEffect(theseInd,:);
% Average across time for each baseline
mBase = cellfun(@(x) mean(x,'omitnan'),thisLastBase,'UniformOutput',false);
%%
k = 1;
f = 190;
figure
hold on
for ii = 1:23
    if thisEffect(ii,k)==1
        for jj = 1:3
            plot(1,mBase{ii,jj}(f),'ok')
        end
    end
    if thisEffect(ii,k)==0
        for jj = 1:3
            plot(1,mBase{ii,jj}(f),'or')
        end
    end
    if thisEffect(ii,k)==-1
        for jj = 1:3
            plot(1,mBase{ii,jj}(f),'oy')
        end
    end
end
%% Combine intervention data by animal and effect
load('F:/irdmRound2/interData.mat')
inter = {'sIL','sNAcC','guan0.3','mph3','mph0.3','mph1','mph0.1',...
    'sNAcC-LSD','sIL-LSD'};
% inter = {'mph3'};
uID = unique(interTable.ID);
delta = [0,1,-1];
nanInds = zeros(287,1);
nanInds([133,240]) = 1;
for ii = 1:numel(inter)
    for k = 1:numel(delta)
        for jj = 1:numel(uID)
            theseInd = contains(interTable.ID,uID{jj}) & ...
                contains(interTable.Rx,inter{ii}) & ...
                arrayfun(@(x) isequal(x,delta(k)),interTable.effect) & ...
                ~nanInds;
            
            data{ii,k,jj} = cat(1,interTable.normData{theseInd});
        end
    end
end
%% Grab 20 samples of null and pos and/or neg effect
[null,pos,neg] = deal(cell(size(data,1),size(data,3)));
for ii = 1:size(data,1)
    for jj = 1:size(data,3)
        % Check for at least 20 null effect
        if size(data{ii,1,jj},1) >= 20
            null{ii,jj} = data{ii,1,jj};
            % Check for at least 20 pos effect
            if size(data{ii,2,jj},1) >= 20
                pos{ii,jj} = data{ii,2,jj};
            end
            % Check for at least 20 neg effect
            if size(data{ii,3,jj},1) >= 20
                neg{ii,jj} = data{ii,3,jj};
            end
        end
    end
end
% save('interEffectData20.mat','null','neg','pos','inter','uID')
%% Grab animals with at least 20 trials of any response
[null,pos,neg] = deal(cell(size(data,1),size(data,3)));
for ii = 1:size(data,1)
    for jj = 1:size(data,3)
        % Check for at least 20 null effect
        if size(data{ii,1,jj},1) >= 20
            null{ii,jj} = data{ii,1,jj};
        end
        % Check for at least 20 pos effect
        if size(data{ii,2,jj},1) >= 20
            pos{ii,jj} = data{ii,2,jj};
        end
        % Check for at least 20 neg effect
        if size(data{ii,3,jj},1) >= 20
            neg{ii,jj} = data{ii,3,jj};
        end
    end
end
%% Add in MPH data using normalized baseline data as null
load('mph3InterDataBase.mat')
mphBaseNorm = cellfun(@zscore,lastBase(~cellfun(@isempty,lastBase)),...
    'UniformOutput',0);
[mphNull,mphPos,mphNeg] = deal(cell(size(data,1),size(data,3)));
for ii = 1:size(data,1)
    for jj = 1:size(data,3)
        mphNull{ii,jj} = mphBaseNorm{jj};
        % Check for at least 20 pos effect
        if size(data{ii,2,jj},1) >= 20
            mphPos{ii,jj} = data{ii,2,jj};
        end
        % Check for at least 20 neg effect
        if size(data{ii,3,jj},1) >= 20
            mphNeg{ii,jj} = data{ii,3,jj};
        end
    end
end
% Combine with other data
load('F:/irdmRound2/interEffectData20.mat')
% Get mph IDs
mphID = unique(interTable.ID);
for ii = 1:numel(mphID)
    % Find index of mphID in uID
    ind = logicFind(mphID{ii},uID,'==');
    null{4,ind} = mphNull{ii};
    pos{4,ind} = mphPos{ii};
    neg{4,ind} = mphNeg{ii};
end
% save('F:/irdmRound2/interEffectData20.mat','-append','null','pos','neg')
%% Build MPH model
t = load('F:/irdmRound2/interEffectData20.mat','null');
nullBase = t.null;
% Get indices of animals in each group (one: all 1s; zero: all 0s; mix)
posLogical = cellfun(@(x) ~isempty(x),pos);
nullLogical = cellfun(@(x) ~isempty(x),null);
posOnly = logicFind(1,double(posLogical(4,:) & ~nullLogical(4,:)),'==');
nullOnly = logicFind(1,double(~posLogical(4,:) & nullLogical(4,:)),'==');
both = logicFind(1,double(posLogical(4,:) & nullLogical(4,:)),'==');
% Get all combinations of using a both and either a pos or null as test set
cmbs = [];
c = 1;
for ii = 1:numel(both)
    for jj = 1:numel(posOnly)
        cmbs(c,:) = [both(ii),posOnly(jj)];
        c = c+1;
    end
    for jj = 1:numel(nullOnly)
        cmbs(c,:) = [both(ii),nullOnly(jj)];
        c = c+1;
    end
end
%% Go through combinations
% Set samples
n = 100;
for ii = 1:size(cmbs,1)
    disp(ii)
    % Grab 100 samples from each animal
    bothX = [];
    bothY = [];
    % First get the pos and null data from the both rat
    thisBothPos = data{4,2,cmbs(ii,1)};
    thisBothNull = data{4,1,cmbs(ii,1)};
    bothX = [thisBothPos(randperm(size(thisBothPos,1),n),:);...
        thisBothNull(randperm(size(thisBothNull,1),n),:)];
    bothY = [ones(n,1);zeros(n,1)];
    % Then get the other animals contribution (either ones or zeros)
    testX = [];
    testY = [];
    if ismember(cmbs(ii,2),posOnly)
        thisPos = data{4,2,cmbs(ii,2)};
        testX = [testX;thisPos(randperm(size(thisPos,1),n),:)];
        testY = [testY;ones(n,1)];
    end
    if ismember(cmbs(ii,2),nullOnly)
        thisNull = data{4,1,cmbs(ii,2)};
        testX = [testX;thisNull(randperm(size(thisNull,1),n),:)];
        testY = [testY;zeros(n,1)];
    end
    testX = [testX;bothX];
    testY = [testY;bothY];
    % Figure out who is left for training
    otherBoth = both(~ismember(both,cmbs(ii,1)));
    otherPos = posOnly(~ismember(posOnly,cmbs(ii,2)));
    otherNull = nullOnly(~ismember(nullOnly,cmbs(ii,2)));
    % Use 100 samples from the both rat
    trainX = [];
    trainY = [];
    trainX = [];
    thisOtherBothPos = cat(1,data{4,2,otherBoth});
    thisOtherBothNull = cat(1,data{4,1,otherBoth});
    trainX = [trainX;thisOtherBothPos(randperm(size(thisOtherBothPos,1),...
        n),:);thisOtherBothNull(randperm(size(thisOtherBothNull,1),n),:)];
    trainY = [trainY;ones(n,1);zeros(n,1)];
    % Figure out number of samples per pos and null
    posN = round(n*(numel(otherNull)+numel(otherPos))/(numel(otherPos)*2));
    nullN = round(n*(numel(otherNull)+numel(otherPos))/(numel(otherNull)*2));
    % Pos - bring in baseline to serve as null data
    for jj = 1:numel(otherPos)
        thisOtherPos = data{4,2,otherPos(jj)};
        thisOtherBaseNull = nullBase{4,otherPos(jj)};
        trainX = [trainX;thisOtherPos(randperm(size(thisOtherPos,1),posN)...
            ,:);thisOtherBaseNull(randperm(size(thisOtherBaseNull,1),posN),:)];
        trainY = [trainY;ones(posN,1);zeros(posN,1)];
    end
    % Null
    for jj = 1:numel(otherNull)
        thisOtherNull = data{4,1,otherNull(jj)};
        if size(thisOtherNull,1)<=nullN
            trainX = [trainX;thisOtherNull];
            trainY = [trainY;zeros(size(thisOtherNull,1),1)];
        else
            trainX = [trainX;thisOtherNull(randperm(size(thisOtherNull,1),...
                nullN),:)];
            trainY = [trainY;zeros(nullN,1)];
        end
    end
    % Fit model
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [thisX,thisY,~,a(ii)] = perfcurve(testY,prob,1);
    x(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,40));
    y(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,40));
    % Combine train and test; then redistribute using 80:20
    allX = [trainX;testX];
    allY = [trainY;trainY];
    [trainX,trainY,testX,testY] = trainTest(allX,allY,0.2);
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    prob = predict(mdl,testX);
    [thisX,thisY,~,a80(ii)] = perfcurve(testY,prob,1);
    x80(ii,:) = interp1(linspace(0,1,numel(thisX)),thisX,linspace(0,1,40));
    y80(ii,:) = interp1(linspace(0,1,numel(thisY)),thisY,linspace(0,1,40));
    % Single feature model - using 80:20 data
    for jj = 1:216
        mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
        prob = predict(mdl,testX(:,jj));
        [~,~,~,aS(ii,jj)] = perfcurve(testY,prob,1);
        coeff(ii,jj) = table2array(mdl.Coefficients(2,1));
    end
end
save('F:/irdmRound2/mphInter.mat','x','y','a','x80','y80','a80','aS',...
    'coeff')
%%
figure
hold on
plot(mean(x,1),mean(y,1))
plot(mean(x80,1),mean(y80,1))
legend({['LOO: ',num2str(round(mean(a),2)),'\pm',...
    num2str(round(conf(a,0.95),2))],['80:20: ',...
    num2str(round(mean(a80),2)),'\pm',num2str(round(conf(a80,0.95),2))]}...
    ,'location','se')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('MPH: change from baseline')
%% Single features
load('F:/irdmRound2/mphInter.mat')
fullFeat = names({'lmPFC','rmPFc','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
[~,sortInds] = sort(mean(aS,1),'descend');
mphA = mean(aS,1).*sign(mean(coeff,1));
mphAsort = mphA(sortInds)';
mphFeat = fullFeat(sortInds)';
%% Build multinomial model
load('F:/irdmRound2/interEffectData20.mat')
[multiDataX,multiDataY] = deal(cell(1,16));
c = 1;
for ii = 1:23
    if ~isempty(null{1,ii})
        if ~any(any(isnan(null{1,ii}))) && ~any(any(isnan(pos{1,ii})))...
                && ~any(any(isnan(neg{1,ii})))
            these = randi(size(null{1,ii},1),1,20);
            thisNull = null{1,ii}(these,:);
            if ~isempty(pos{1,ii})
                these = randi(size(pos{1,ii},1),1,20);
                thisPos = pos{1,ii}(these,:);
            else
                thisPos = [];
            end
            if ~isempty(neg{1,ii})
                these = randi(size(neg{1,ii},1),1,20);
                thisNeg = neg{1,ii}(these,:);
            else
                thisNeg = [];
            end
            multiDataX{c} = [thisNull;thisPos;thisNeg];
            multiDataY{c} = [ones(size(thisNull,1),1);...
                ones(size(thisPos,1),1)*2;...
                ones(size(thisNeg,1),1)*3];
            c = c+1;
        end
    end
end
%% Generate 100 sets of samples for individual models
load('F:/irdmRound2/interEffectData20.mat')
% sIL, sNAcC, mph3 | inter([1,2,4])
nullN = cellfun(@(x) size(x,1),null);
posN = cellfun(@(x) size(x,1),pos);
negN = cellfun(@(x) size(x,1),neg);
% Which animals have enough samples to be used
posInds = nullN>0 & posN>0;
negInds = nullN>0 & negN>0;
nullInds = nullN>0 & (negN>0 | posN>0);
c = 1;
for interI = [1,2,4]
    % Positive model
    for animI = logicFind(1,posInds(interI,:),'==')
        posSamp{c,animI} = randi(size(pos{interI,animI},1),100,20);
    end
    % Negative model
    for animI = logicFind(1,negInds(interI,:),'==')
        negSamp{c,animI} = randi(size(neg{interI,animI},1),100,20);
    end
    % Null
    for animI = logicFind(1,nullInds(interI,:),'==')
        nullSamp{c,animI} = randi(size(null{interI,animI},1),100,20);
    end
    c = c+1;
end
% save('interEffectSamples20.mat','nullSamp','negSamp','posSamp',...
%     'posInds','negInds','nullInds')
%% Load models
negCounter = 1;
for ii = [1:8,10:21]
    load(['F:/irdmRound2/interEffect2/interEffect',num2str(ii),'.mat'])
    for jj = 1:3
        posAnimal = logicFind(1,~cellfun(@isempty,posAcc(jj,:)),'==');
        if jj == 1
            % Load over the corresponding neg data
            load(['F:/irdmRound2/interEffectNeg/interEffectNeg',...
                num2str(negCounter),'.mat'])
            negAnimal = logicFind(1,~cellfun(@isempty,negAcc(jj,:)),'==');
        end
        for k = 1:23
            if ~isempty(posAcc{jj,k})
                posA{jj,k}(ii) = posAcc{jj,k}{1}.acc;
                posAR{jj,k}(ii) = posAccR{jj,k}{1}.acc;
                if ~isempty(posLOOAcc{jj,k})
                    posLOOA{jj,k}(ii) = posLOOAcc{jj,k}{1}.acc;
                    posLOOX{jj,k,ii} =  interp1(linspace(0,1,...
                            numel(posLOOAcc{jj,k}{1}.x)),...
                            posLOOAcc{jj,k}{1}.x,linspace(0,1,40));
                    posLOOY{jj,k,ii} = interp1(linspace(0,1,...
                            numel(posLOOAcc{jj,k}{1}.y)),...
                            posLOOAcc{jj,k}{1}.y,linspace(0,1,40));
                end
                % Test model on other animals
                thisMdl = posAcc{jj,k}{1}.mdl{1};
                others = posAnimal(logicFind(1,...
                    ~ismember(posAnimal,k),'=='));
                for oi = 1:numel(others)
                    % Skip if missing data
                    if size(posHist{jj,others(oi)}.cfg.naive.testX,2)==size(posHist{jj,k}.cfg.naive.testX,2)
                        p = cvglmnetPredict(thisMdl,posHist{jj,others(oi)}.cfg.naive.testX);
                        % Stored such that oPosA{1,2,3}{5} is the 3rd
                        % iteration, of the 2nd animal, being tested on the
                        % 5th animal, in the 1st intervention
                        [~,~,~,oPosA{jj,k,ii}{oi}] = perfcurve(posHist{jj,...
                            others(oi)}.cfg.naive.testY,p,1);
                    else
                         oPosA{jj,k,ii}{oi} = NaN;
                    end
                end
            end
            % Skip k 22,23 and jj 2,3 since no neg data
            if k <= 21 && jj == 1 
                if ~isempty(negAcc{jj,k})
                    negA{jj,k}(ii) = negAcc{jj,k}{1}.acc;
                    negAR{jj,k}(ii) = negAccR{jj,k}{1}.acc;
                    if jj == 1
                        negLOOA{jj,k}(ii) = negLOOAcc{jj,k}{1}.acc;
                        negLOOX{jj,k,ii} = interp1(linspace(0,1,...
                            numel(negLOOAcc{jj,k}{1}.x)),...
                            negLOOAcc{jj,k}{1}.x,linspace(0,1,40));
                        negLOOY{jj,k,ii} = interp1(linspace(0,1,...
                            numel(negLOOAcc{jj,k}{1}.y)),...
                            negLOOAcc{jj,k}{1}.y,linspace(0,1,40));
                    end
                    % Test model on other animals
                    thisMdl = negAcc{jj,k}{1}.mdl{1};
                    others = negAnimal(logicFind(1,...
                        ~ismember(negAnimal,k),'=='));
                    % Only stimIL has enough animals
                    if jj == 1
                        for oi = 1:numel(others)
                            if size(negHist{jj,others(oi)}.cfg...
                                    .naive.testX,2) == size(...
                                    negHist{jj,k}.cfg.naive.testX,2)
                                p = cvglmnetPredict(thisMdl,negHist{jj,...
                                    others(oi)}.cfg.naive.testX);
                                [~,~,~,oNegA{jj,k,ii}{oi}] = perfcurve(negHist{jj,...
                                    others(oi)}.cfg.naive.testY,p,1);
                            else
                                oNegA{jj,k,ii}{oi} = NaN;
                            end
                        end
                    end
                end
            end
        end
        posAllA(jj,ii) = allPosAcc{jj}{1}.acc;
        posAllX{jj,ii} = allPosAcc{jj}{1}.x;
        posAllY{jj,ii} = allPosAcc{jj}{1}.y;
        % Test on randomized data
        pred = cvglmnetPredict(allPosAcc{jj}{1}.mdl{1},...
            allPosHist{jj}.cfg.naive.testX);
        [~,~,~,posAllAR(jj,ii)] = perfcurve(...
            allPosHist{jj}.cfg.naive.testY(randi(numel(...
            allPosHist{jj}.cfg.naive.testY),1,...
            numel(allPosHist{jj}.cfg.naive.testY))),pred,1);
        if jj == 1
            negAllA(jj,ii) = allNegAcc{jj}{1}.acc;
            negAllX{jj,ii} = allNegAcc{jj}{1}.x;
            negAllY{jj,ii} = allNegAcc{jj}{1}.y;
            
            pred = cvglmnetPredict(allNegAcc{jj}{1}.mdl{1},...
                allNegHist{jj}.cfg.naive.testX);
            [~,~,~,negAllAR(jj,ii)] = perfcurve(...
                allNegHist{jj}.cfg.naive.testY(randi(numel(...
                allNegHist{jj}.cfg.naive.testY),1,...
                numel(allNegHist{jj}.cfg.naive.testY))),pred,1);
        end
    end
    negCounter = negCounter+1;
end
posAM = cellfun(@mean,posA);
posALOOM = cellfun(@mean,posLOOA);
negAM = cellfun(@mean,negA);
negALOOM = cellfun(@mean,negLOOA);
%% Load permuted feature models
freq = [];
site = [];
for ii = 1:5
    load(['F:/irdmRound2/interEffect2/permuted/interEffectPerm',num2str(ii),'.mat'])
    full(ii) = posAcc{1,2}{1}.acc;
    freq(ii,:) = squeeze(cell2mat(permFreqAcc(1,2,:)));
    site(ii,:) = squeeze(cell2mat(permSiteAcc(1,2,:)));
end
% freqChange = (mean(full)-mean(freq))/std(full);
% siteChange = (mean(full)-mean(site))/std(full);
freqChange = mean(full)-mean(freq);
siteChange = mean(full)-mean(site);
for ii = 1:size(freq,2)
    [~,p(ii)] = ttest2(full,freq(:,ii));
end
for ii = 1:size(site,2)
    [~,p(ii)] = ttest2(full,site(:,ii));
end
fullFeat = names({'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
fullFeat = fullFeat(1:216);
sites = {'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'};
%% 80:20 IL performance with permuted features
for ii = 1:20
    load(['F:/irdmRound2/interEffectGroupPerm/interEffectGroupPerm',...
        num2str(ii),'.mat'])
    allAcc(ii) = allPosAcc{1}{1}.acc;
    freq(ii,:) = cell2mat(permFreqAcc(1,22,:));
    site(ii,:) = cell2mat(permSiteAcc(1,22,:));
end
mean(allAcc);
mean(allAcc)-mean(freq);
mean(allAcc)-mean(site);
mA = posAllA(1,[1:8,10:21]);

%%
doubleHist(allAcc,permCorticalAcc)
%% Test models on individual test sets (gen to each animal)
load('interEffectData20.mat')
load('interEffectSamples20.mat')
% Pos
for k = 1:3
    posX = cat(1,pos{k,:});
    nullX = cat(1,null{k,posInds(k,:)});
    x = [posX;nullX];
    sizes = [cellfun(@(x) size(x,1),pos(k,posInds(k,:))),...
        cellfun(@(x) size(x,1),null(k,posInds(k,:)))];
    id = [1:numel(sizes)/2,1:numel(sizes)/2];
    ids = [];
    for ii = 1:numel(id)
        ids = [ids;repmat(id(ii),sizes(ii),1)];
    end
    for ii = [1:8,10:21]
        load(['F:/irdmRound2/interEffect2/interEffect',num2str(ii),'.mat'],...
            'allPosAcc','allPosHist')
        testID = [];
        for jj = 1:size(allPosHist{k}.cfg.naive.testX,1)
            testID(jj) = ids(logicFind(allPosHist{k}.cfg.naive.testX(jj,1),x(:,1),'=='));
        end
        for jj = 1:numel(sizes)/2
            thesePred = allPosAcc{k}{1}.pred{1}(testID==jj);
            if ~isempty(thesePred)
                [~,~,~,thisA(ii,jj,k)] = perfcurve(allPosHist{k}.cfg.naive.testY(testID==jj),thesePred,1);
            end
        end
    end
end

%% Overall stim prediction performance LOO
inter = {'sIL+ ','sNAcC+ ','MPH+ ','sIL- '};
figure
for ii = 1:3
    subplot(2,2,ii)
    hold on
    c = 1;
    thisAnimalX = [];
    thisAnimalY = [];
    for jj = 1:23
        if ~isempty(posLOOX{ii,jj,2})
            thisAnimalX(c,:) = mean(cat(1,posLOOX{ii,jj,:}),1);
            thisAnimalY(c,:) = mean(cat(1,posLOOY{ii,jj,:}),1);
            plot(thisAnimalX(c,:),thisAnimalY(c,:),'color',[0.8 0.8 0.8])
            c = c+1;
        end
    end
    loox = mean(thisAnimalX,1);
    looy = mean(thisAnimalY,1);
    plot(loox,looy,'-k')
    
    plot(mean(cat(2,posAllX{ii,:}),2),mean(cat(2,posAllY{ii,:}),2),'--k')
    set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
    box off
    title([inter{ii},num2str(round(mean(posALOOM(ii,:),'omitnan'),2)),'\pm',...
        num2str(round(conf(posALOOM(ii,~isnan(posALOOM(ii,:))),0.95),2))])
    disp([inter{ii},' permuted = ',num2str(mean(posAllAR(ii,[1:8,10:21]))),'\pm',num2str(conf(posAllAR(ii,[1:8,10:21]),0.95))])
    disp([inter{ii},'all = ',num2str(mean(posAllA(ii,[1:8,10:21]))),'\pm',num2str(conf(posAllA(ii,[1:8,10:21]),0.95))])
end

subplot(2,2,4)
hold on
c = 1;
thisAnimalX = [];
thisAnimalY = [];
for jj = 1:21
    if ~isempty(negLOOX{1,jj,2})
        thisAnimalX(c,:) = mean(cat(1,negLOOX{1,jj,:}),1);
        thisAnimalY(c,:) = mean(cat(1,negLOOY{1,jj,:}),1);
        plot(thisAnimalX(c,:),thisAnimalY(c,:),'color',[0.8 0.8 0.8])
        c = c+1;
    end
end
loox = mean(thisAnimalX,1);
looy = mean(thisAnimalY,1);
plot(loox,looy,'-k')

plot(mean(cat(2,negAllX{:}),2),mean(cat(2,negAllY{:}),2),'--k')
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
box off
title([inter{4},num2str(round(mean(negALOOM(1,:),'omitnan'),2)),'\pm',...
    num2str(round(conf(negALOOM(1,~isnan(negALOOM(1,:))),0.95),2))])
disp(['sIL- permuted = ',num2str(mean(negAllAR)),'\pm',num2str(conf(negAllAR,0.95))])
disp(['sIL- all = ',num2str(mean(negAllA)),'\pm',num2str(conf(negAllA,0.95))])
%%
figure
hold on
plot(mean(x,1),mean(y,1))
plot(mean(xP,1),mean(yP,1),'--')
for ii = 1:3
    plot(squeeze(mean(subX(ii,:,:))),squeeze(mean(subY(ii,:,:))))
end
legend({['all: ',num2str(round(mean(allA),2)),'\pm',...
    num2str(round(conf(allA,0.95),2))],['perm: ',...
    num2str(round(mean(allAP),2)),'\pm',...
    num2str(round(conf(allAP,0.95),2))],['sIL: ',...
    num2str(round(mean(subA(:,1))',2)),'\pm',...
    num2str(round(conf(subA(:,1)',0.95),2))],['sCore: ',...
    num2str(round(mean(subA(:,2))',2)),'\pm',...
    num2str(round(conf(subA(:,2)',0.95),2))],['MPH: ',...
    num2str(round(mean(subA(:,3))',2)),'\pm',...
    num2str(round(conf(subA(:,3)',0.95),2))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
%% Stim effect prediction - cross
load('F:/irdmRound2/lastBaseNorm.mat')
ids = uID(~cellfun(@isempty,lastBase(:,1)));
%% Waffle plot - sIL pos
self = cat(1,posA{1,:});
sILpos = zeros(9,9);
these = [1:3,5,8,15:16,21:22];
this = squeeze(oPosA(1,these,:));
for ii = 1:9
    sILpos(ii,:) = mean(self(ii,:));
    o = logicFind(1,~ismember(1:9,ii),'==');
    sILpos(ii,o) = mean(cell2mat(cat(1,this{ii,:})),1,'omitnan');
end
figure
subplot(2,2,1)
sILpos = sILpos([4,5,1,8,2,9,3,6,7],[4,5,1,8,2,9,3,6,7]);
pcolor(padarray(sILpos,[1,1],'post'))
colormap viridis
xlabel('train')
ylabel('test')
title('sIL +')
set(gca,'xtick',1.5:9.5,'xticklabel',ids(these([4,5,1,8,2,9,3,6,7])),...
    'ytick',1.5:9.5,'yticklabel',ids(these([4,5,1,8,2,9,3,6,7])))
xtickangle(45)
caxis([0 1])
box off

% Get sex; 1 = male
sex = [0,0,0,0,0,1,1,1,1];
% Get left or right; left = 1
leftRight = [1,1,0,0,0,1,1,0,0];
% Set data to use; where dim = test x train
this = sILpos;
% male->female
mf = this(sex==0,sex==1);
% male->male
mm = this(sex==1,sex==1);
mm(1:size(mm,1)+1:end) = NaN;
% female->female
ff = this(sex==0,sex==0);
ff(1:size(ff,1)+1:end) = NaN;
% female->male
fm = this(sex==1,sex==0);

totalMM = sum(sum(~isnan(mm)));
totalFF = sum(sum(~isnan(ff)));
totalMF = sum(sum(~isnan(mf)));
totalFM = sum(sum(~isnan(fm)));

p = 0.7;
nMM = sum(sum(mm>p));
nFF = sum(sum(ff>p));
nMF = sum(sum(mf>p));
nFM = sum(sum(fm>p));

% left->left
ll = this(leftRight==1,leftRight==1);
ll(1:size(ll,1)+1:end) = NaN;
% left->right
lr = this(leftRight==0,leftRight==1);
% right->left 
rl = this(leftRight==1,leftRight==0);
% right->right
rr = this(leftRight==0,leftRight==0);
rr(1:size(rr,1)+1:end) = NaN;

totalLL = sum(sum(~isnan(ll)));
totalRR = sum(sum(~isnan(rr)));
totalLR = sum(sum(~isnan(lr)));
totalRL = sum(sum(~isnan(rl)));

nLL = sum(sum(ll>p));
nRR = sum(sum(rr>p));
nLR = sum(sum(lr>p));
nRL = sum(sum(rl>p));

perc(1,:) = [nMM/totalMM,nMF/totalMF,nFF/totalFF,nFM/totalFM,...
    nLL/totalLL,nLR/totalLR,nRR/totalRR,nRL/totalRL];
n(1,:) = [nMM,nMF,nFF,nFM,nLL,nLR,nRR,nRL];
total(1,:) = [totalMM,totalMF,totalFF,totalFM,totalLL,totalLR,totalRR,...
    totalRL];
% Waffle plot - sIL neg
self = cat(1,negA{1,:});
sILneg = zeros(9,9);
these = [5,7,8,11,15,17,19:21];
this = squeeze(oNegA(1,these,:));
for ii = 1:9
    sILneg(ii,ii) = mean(self(ii,:));
    o = logicFind(1,~ismember(1:9,ii),'==');
    sILneg(ii,o) = mean(cell2mat(cat(1,this{ii,:})),1);
end
subplot(2,2,2)
pcolor(padarray(sILneg,[1,1],'post'))
colormap viridis
xlabel('train')
ylabel('test')
title('sIL -')
set(gca,'xtick',1.5:9.5,'xticklabel',ids(these([1:3,8,4,7,9,6,5])),'ytick',1.5:9.5,...
    'yticklabel',ids(these([1:3,8,4,7,9,6,5])))
xtickangle(45)
caxis([0 1])
box off

% Get sex; 1 = male
sex = [0,0,1,1,1,1,1,1,1];
% Get left or right; left = 1
leftRight = [1,1,1,1,1,1,0,0,0];
% Set data to use; where dim = test x train
this = sILneg;
% male->female
mf = this(sex==0,sex==1);
% male->male
mm = this(sex==1,sex==1);
mm(1:size(mm,1)+1:end) = NaN;
% female->female
ff = this(sex==0,sex==0);
ff(1:size(ff,1)+1:end) = NaN;
% female->male
fm = this(sex==1,sex==0);

totalMM = sum(sum(~isnan(mm)));
totalFF = sum(sum(~isnan(ff)));
totalMF = sum(sum(~isnan(mf)));
totalFM = sum(sum(~isnan(fm)));

nMM = sum(sum(mm>p));
nFF = sum(sum(ff>p));
nMF = sum(sum(mf>p));
nFM = sum(sum(fm>p));

% left->left
ll = this(leftRight==1,leftRight==1);
ll(1:size(ll,1)+1:end) = NaN;
% left->right
lr = this(leftRight==0,leftRight==1);
% right->left 
rl = this(leftRight==1,leftRight==0);
% right->right
rr = this(leftRight==0,leftRight==0);
rr(1:size(rr,1)+1:end) = NaN;

totalLL = sum(sum(~isnan(ll)));
totalRR = sum(sum(~isnan(rr)));
totalLR = sum(sum(~isnan(lr)));
totalRL = sum(sum(~isnan(rl)));

nLL = sum(sum(ll>p));
nRR = sum(sum(rr>p));
nLR = sum(sum(lr>p));
nRL = sum(sum(rl>p));

perc(2,:) = [nMM/totalMM,nMF/totalMF,nFF/totalFF,nFM/totalFM,...
    nLL/totalLL,nLR/totalLR,nRR/totalRR,nRL/totalRL];
n(2,:) = [nMM,nMF,nFF,nFM,nLL,nLR,nRR,nRL];
total(2,:) = [totalMM,totalMF,totalFF,totalFM,totalLL,totalLR,totalRR,...
    totalRL];
% Waffle plot - core +
self = cat(1,posA{2,:});
corePos = zeros(size(self,1),size(self,1));
these = [5,9,12,21,23];
this = squeeze(oPosA(2,these,:));
for ii = 1:size(self,1)
    corePos(ii,ii) = mean(self(ii,:));
    o = logicFind(1,~ismember(1:size(self,1),ii),'==');
    corePos(ii,o) = mean(cell2mat(cat(1,this{ii,:})),1,'omitnan');
end
subplot(2,2,3)
pcolor(padarray(corePos,[1,1],'post'))
colormap viridis
xlabel('train')
ylabel('test')
title('sNAcC +')
set(gca,'xtick',1.5:9.5,'xticklabel',ids(these([1,2,4,3,5])),'ytick',1.5:9.5,...
    'yticklabel',ids(these([1,2,4,3,5])))
xtickangle(45)
caxis([0 1])
box off

% Get sex; 1 = male
sex = [1,1,0,0,0];
% Get left or right; left = 1
leftRight = [1,1,1,0,0];
% Set data to use; where dim = test x train
this = corePos;
% male->female
mf = this(sex==0,sex==1);
% male->male
mm = this(sex==1,sex==1);
mm(1:size(mm,1)+1:end) = NaN;
% female->female
ff = this(sex==0,sex==0);
ff(1:size(ff,1)+1:end) = NaN;
% female->male
fm = this(sex==1,sex==0);

totalMM = sum(sum(~isnan(mm)));
totalFF = sum(sum(~isnan(ff)));
totalMF = sum(sum(~isnan(mf)));
totalFM = sum(sum(~isnan(fm)));

nMM = sum(sum(mm>p));
nFF = sum(sum(ff>p));
nMF = sum(sum(mf>p));
nFM = sum(sum(fm>p));

% left->left
ll = this(leftRight==1,leftRight==1);
ll(1:size(ll,1)+1:end) = NaN;
% left->right
lr = this(leftRight==0,leftRight==1);
% right->left 
rl = this(leftRight==1,leftRight==0);
% right->right
rr = this(leftRight==0,leftRight==0);
rr(1:size(rr,1)+1:end) = NaN;

totalLL = sum(sum(~isnan(ll)));
totalRR = sum(sum(~isnan(rr)));
totalLR = sum(sum(~isnan(lr)));
totalRL = sum(sum(~isnan(rl)));

nLL = sum(sum(ll>p));
nRR = sum(sum(rr>p));
nLR = sum(sum(lr>p));
nRL = sum(sum(rl>p));

perc(3,:) = [nMM/totalMM,nMF/totalMF,nFF/totalFF,nFM/totalFM,...
    nLL/totalLL,nLR/totalLR,nRR/totalRR,nRL/totalRL];
n(3,:) = [nMM,nMF,nFF,nFM,nLL,nLR,nRR,nRL];
total(3,:) = [totalMM,totalMF,totalFF,totalFM,totalLL,totalLR,totalRR,...
    totalRL];
% MPH 3 mg/kg
self = cat(1,posA{3,:});
mphPos = zeros(size(self,1),size(self,1));
these = [18,21];
this = squeeze(oPosA(3,these,:));
for ii = 1:size(self,1)
    mphPos(ii,ii) = mean(self(ii,:));
    o = logicFind(1,~ismember(1:size(self,1),ii),'==');
    mphPos(ii,o) = mean(cell2mat(cat(1,this{ii,:})),1,'omitnan');
end
subplot(2,2,4)
pcolor(padarray(mphPos,[1,1],'post'))
colormap viridis
xlabel('train')
ylabel('test')
title('mph+')
set(gca,'xtick',1.5:3.5,'xticklabel',ids(these),'ytick',1.5:3.5,...
    'yticklabel',ids(these))
xtickangle(45)
caxis([0 1])
box off
colorbar

%% Single feature
negCounter = 1;
for ii = [1:8,10:21]
    load(['F:/irdmRound2/interEffect2/interEffect',num2str(ii),'.mat'],...
        'aPos','aNeg','coeffPos','coeffNeg','allPosHist')
    % Group models
    for k = 1:3
        for jj = 1:size(allPosHist{k}.trainX,2)
            mdl = fitglm(allPosHist{k}.trainX(:,jj),...
                allPosHist{k}.trainY,'distribution','binomial');
            posCGroup{k}(jj,ii) = table2array(mdl.Coefficients(2,1));
            pred = predict(mdl,allPosHist{k}.cfg.naive.testX(:,jj));
            [~,~,~,posGroup{k}(jj,ii)] = perfcurve(...
                allPosHist{k}.cfg.naive.testY,pred,1);
        end
    end
    load(['F:/irdmRound2/interEffectNeg/interEffectNeg',...
        num2str(negCounter),'.mat'],...
        'allNegHist')
    for jj = 1:size(allNegHist{1}.trainX,2)
        mdl = fitglm(allNegHist{1}.trainX(:,jj),allNegHist{1}.trainY,...
            'distribution','binomial');
        negCGroup(jj,negCounter) = table2array(mdl.Coefficients(2,1));
        pred = predict(mdl,allNegHist{1}.cfg.naive.testX(:,jj));
        [~,~,~,negGroup(jj,negCounter)] = perfcurve(...
            allNegHist{1}.cfg.naive.testY,pred,1);
    end
    % Individual models
    posInd(:,:,:,ii) = aPos;
    negInd(:,:,:,negCounter) = aNeg;
    posCInd(:,:,:,ii) = coeffPos;
    negCInd(negCounter,:,:) = coeffNeg;
    negCounter = negCounter+1;
end
%% Plot group model features
% Turn coeffs into sign then add to auc
for ii = 1:3
    groupPosSign = posCGroup{ii}(:,[1:8,10:21])./...
        abs(posCGroup{ii}(:,[1:8,10:21]));
    groupPosAUC(:,ii) = mean(posGroup{ii}(:,[1:8,10:21]).*groupPosSign,2);
end
groupNegSign = negCGroup./abs(negCGroup);
groupNegAUC = mean(negGroup.*groupNegSign,2);

fullFeat = names({'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
%% Comparison plot of all AUC+ models
c = distinguishable_colors(6);
inter = {'sIL+','sNAcC+','MPH+','sIL-'};
cmbs = [1,2;1,3;2,3];
counter = 1;

figure
for k = 1:3
    subplot(2,2,k)
    hold on
    scatter(groupPosAUC(1:48,cmbs(k,1)),groupPosAUC(1:48,cmbs(k,2)),[],...
        repmat(c,8,1),'o')
    scatter(groupPosAUC(49:end,cmbs(k,1)),groupPosAUC(49:end,cmbs(k,2)),...
        [],repmat(c,28,1),'s')
    xlim([-1 1]); ylim([-1 1])
    set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
    plot([-1 -0.5],[-0.5 -0.5],':k')
    plot([1 0.5],[-0.5 -0.5],':k')
    plot([-1 -0.5],[0.5 0.5],':k')
    plot([1 0.5],[0.5 0.5],':k')
    
    plot([-0.5 -0.5],[-1 -0.5],':k')
    plot([-0.5 -0.5],[1 0.5],':k')
    plot([0.5 0.5],[-1 -0.5],':k')
    plot([0.5 0.5],[1 0.5],':k')
    xlabel(inter{cmbs(k,1)})
    ylabel(inter{cmbs(k,2)})
    counter = counter+1;
end

subplot(2,2,4)
hold on
scatter3(groupPosAUC(1:48,1),groupPosAUC(1:48,2),groupPosAUC(1:48,3),[],...
    repmat(c,8,1),'o')
scatter3(groupPosAUC(49:end,1),groupPosAUC(49:end,2),...
    groupPosAUC(49:end,3),[],repmat(c,28,1),'s')
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1,'ztick',-1:0.5:1)
xlabel(inter{1}); ylabel(inter{2}); zlabel(inter{3})
view(-99,29)
% Also table of top features from each
[topPos,topPosInd] = sort(abs(groupPosAUC),'descend');
topPos(topPosInd) = topPos(topPosInd).*sign(groupPosAUC);
for ii = 1:3
    topPosFeat(:,ii) = fullFeat(topPosInd(:,ii));
end
% sIL- model; compare to sIL+
figure
hold on
scatter(groupNegAUC(1:48),groupPosAUC(1:48,1),[],repmat(c,8,1),'o')
scatter(groupNegAUC(49:end),groupPosAUC(49:end,1),[],repmat(c,28,1),'s')
xlim([-1 1]); ylim([-1 1])
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
plot([-1 -0.5],[-0.5 -0.5],':k')
plot([1 0.5],[-0.5 -0.5],':k')
plot([-1 -0.5],[0.5 0.5],':k')
plot([1 0.5],[0.5 0.5],':k')

plot([-0.5 -0.5],[-1 -0.5],':k')
plot([-0.5 -0.5],[1 0.5],':k')
plot([0.5 0.5],[-1 -0.5],':k')
plot([0.5 0.5],[1 0.5],':k')
xlabel('sIL-'); ylabel('sIL+')
% And get top neg features
[topNeg,topNegInd] = sort(abs(groupNegAUC),'descend');
topNeg = topNeg.*sign(groupNegAUC);
topNegFeat = fullFeat(topNegInd)';
%% Plot sorted single features
figure
plot(2:217,abs(topPos(:,1)),'.k')
hold on
c = conf(posGroup{1}(:,[1:8,10:21]),0.95);
for ii = 1:numel(c)
    plot([ii+1 ii+1],[abs(topPos(ii,1))-c(ii) abs(topPos(ii,1))+c(ii)],'-k')
end
plot(1,mean(posAllA(1,[1:8,10:21])))
plot([1 1],[mean(posAllA(1,[1:8,10:21]))-conf(posAllA(1,[1:8,10:21]),0.95)...
    mean(posAllA(1,[1:8,10:21]))+conf(posAllA(1,[1:8,10:21]),0.95)],'-k')
plot([1 217],[mean(posAllAR(1,[1:8,10:21])) mean(posAllAR(1,[1:8,10:21]))],':k')
box off
xlabel('feature')
ylabel('performance (AUC)')
%% Grabbed data with data vis tool (brush)
allPosFeat = cell(1); allNegFeat = cell(1);
for ii = 1:size(allPos,1)
    allPosFeat{ii,1} = fullFeat(logicFind(1,sum(...
        groupPosAUC==allPos(ii,:),2)==3,'=='));
end
for ii = 1:size(allNeg,1)
    allNegFeat{ii,1} =  fullFeat(logicFind(1,sum(...
        groupPosAUC==allNeg(ii,:),2)==3,'=='));
end
%%
this = [groupNegAUC,groupPosAUC(:,1)];
over60 = round(this(sum(abs(this)>0.6,2)==2,:),2);
over60Feat = fullFeat(sum(abs(this)>0.6,2)==2)';
%% Plot individual model features
mPosAll = squeeze(mean(posInd,4));
mNegAll = squeeze(mean(negInd,4));
posSign = sign(mean(posCInd,4));
negSign = squeeze(sign(mean(negCInd,1)));

indsPos = cell(3);
indsNeg = zeros(9,216);
for ii = 1:3
    c1 = 1; c2 = 1;
    for jj = 1:23
        if sum(mPosAll(ii,jj,:)) ~= 0
            [posAsort{ii}(c1,:),indsPos{ii}(c1,:)] = sort(squeeze(...
                mPosAll(ii,jj,:)),...
                'descend');
            c1 = c1+1;
        end
        if ii == 1 && jj <= 21
            if sum(mNegAll(jj,:)) ~= 0
                [negAsort(c2,:),indsNeg(c2,:)] = sort(squeeze(...
                    mNegAll(jj,:)),'descend');
                c2 = c2+1;
            end
        end
    end
end
% Full features
fullFeat = names({'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
% Missing ch 8
sevenChFeat = names({'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS',...
    'lNAcC'},{'d','t','a','b','lg','hg'});
% Missing first two ch
fiveChFeat = names({'lOFC','rOFC','lNAcS','rNAcS','lNAcC','rNAcC'},...
    {'d','t','a','b','lg','hg'});
% Grab top 10 features of each
top5Pos{1} = sevenChFeat(indsPos{1}(1,1:10));
for ii = 2:8
    top5Pos{1}(ii,:) = fullFeat(indsPos{1}(ii,1:10));
end
top5Pos{1}(9,:) = sevenChFeat(indsPos{1}(9,1:10));
for ii = 2:3
    for jj = 1:size(indsPos{ii},1)
        top5Pos{ii}(jj,:) = fullFeat(indsPos{ii}(jj,1:10));
    end
end
top5Neg = fullFeat(indsNeg(:,1:10));
% Plot features
% Add sign back into feature performance
load('F:/irdmRound2/lastBaseNorm.mat')
ids = uID(~cellfun(@isempty,lastBase(:,1)));
[this,sortInd,thisSort] = deal(cell(1,3));

posInter = {'sIL+','sNAcC+','MPH+'};
theseID = {[1:3,5,8,15:16,21:22],[5,9,12,21,23],[18,21]};
for k = 1:3
    uFeat{k} = unique(top5Pos{k});
    this{k} = zeros(size(top5Pos{k},1),numel(uFeat{k}));
    for ii = 1:size(top5Pos{k},1)
        for jj = 1:numel(uFeat{k})
            if any(strcmp(top5Pos{k}(ii,:),uFeat{k}{jj}))
                this{k}(ii,jj) = posAsort{k}(strcmp(top5Pos{k}(ii,:),...
                    uFeat{k}{jj}))*posSign(k,theseID{k}(ii),...
                    indsPos{k}(ii,strcmp(top5Pos{k}(ii,:),uFeat{k}{jj})));
            else
                this{k}(ii,jj) = NaN;
            end
        end
    end
    
    [~,sortInd{k}] = sort(sum(this{k}~=0,1),'descend');
    thisSort{k} = this{k}(:,sortInd{k});
    
    figure
    pcolor(padarray(thisSort{k},[1,1],'post'))
    set(gca,'xtick',1.5:numel(uFeat{k})+1,'xticklabel',...
        uFeat{k}(sortInd{k}),'ytick',1.5:size(thisSort{k},1)+1,...
        'yticklabel',ids(theseID{k}))
    xtickangle(45)
    title(posInter{k})
    caxis([-1 1])
    colormap viridis
    box off
end
% Neg
uFeat = unique(top5Neg);
this = zeros(size(top5Neg,1),numel(uFeat));
theseID = [5,7,8,11,15,17,19:21];
for ii = 1:size(top5Neg,1)
    for jj = 1:numel(uFeat)
        if any(strcmp(top5Neg(ii,:),uFeat{jj}))
            this(ii,jj) = negAsort(strcmp(top5Neg(ii,:),...
                uFeat{jj}))*negSign(theseID(ii),indsNeg(ii,...
                strcmp(top5Neg(ii,:),uFeat{jj})));
        else
            this(ii,jj) = NaN;
        end
    end
end
[~,sortInd] = sort(sum(this~=0,1),'descend');
thisSort = this(:,sortInd);

figure
pcolor(padarray(thisSort,[1,1],'post'))
set(gca,'xtick',1.5:numel(uFeat)+1,'xticklabel',uFeat(sortInd),...
    'ytick',1.5:size(thisSort,1)+1,'yticklabel',ids(theseID))
xtickangle(45)
title('sIL-')
caxis([-1 1])
colorbar
colormap viridis
box off
%% Combined inter effect pos models
load('F:/irdmRound2/interEffectData20.mat')
load('F:/irdmRound2/interEffectSamples20.mat')
[x,y,xP,yP] = deal([]);
count = 1;
for ii = [1:100]
    [allX,allY,allGroup,allAnimal] = deal([]);
    c = 1;
    for interI = [1,2,4]
        for jj = logicFind(1,posInds(interI,:),'==')
            thisPosX = pos{interI,jj}(posSamp{c,jj}(ii,:),:);
            thisPosY = ones(20,1);
            thisNullX = null{interI,jj}(nullSamp{c,jj}(ii,:),:);
            thisNullY = zeros(20,1);
            thisGroup = repmat(interI,40,1);
            if ~any(isnan(thisPosX),[1,2])
                allX = [allX;thisPosX;thisNullX];
                allY = [allY;thisPosY;thisNullY];
                allGroup = [allGroup;thisGroup];
                allAnimal = [allAnimal;repmat(jj,40,1)];
            end
        end
        c = c+1;
    end
    load(['F:/irdmRound2/interEffectPosCmb/interEffectPosCmb',...
        num2str(ii),'.mat'])
    allA(count) = acc{1}.acc;
    x(count,:) = interp1(linspace(0,1,numel(acc{1}.x)),acc{1}.x,...
        linspace(0,1,185));
    y(count,:) = interp1(linspace(0,1,numel(acc{1}.y)),acc{1}.y,...
        linspace(0,1,185));
    allAP(count) = accP{1}.acc;
    xP(count,:) = interp1(linspace(0,1,numel(accP{1}.x)),accP{1}.x,...
        linspace(0,1,185));
    yP(count,:) = interp1(linspace(0,1,numel(accP{1}.y)),accP{1}.y,...
        linspace(0,1,185));
    % Get performance on each group alone
    group = [];
    for jj = 1:size(hist.cfg.naive.testX,1)
        [~,idx] = ismember(hist.cfg.naive.testX(jj,:),allX,'rows');
        group(jj) = allGroup(idx);
    end
    c = 1;
    for interI = [1,2,4]
        theseX = hist.cfg.naive.testX(group==interI,:);
        theseY = hist.cfg.naive.testY(group==interI,:);
        pred = cvglmnetPredict(acc{1}.mdl{1},theseX);
        [theseX,theseY,~,subA(count,c)] = perfcurve(theseY,pred,1);
        subX(c,count,:) = interp1(linspace(0,1,numel(theseX)),...
            theseX,linspace(0,1,50));
        subY(c,count,:) = interp1(linspace(0,1,numel(theseY)),...
            theseY,linspace(0,1,50));
        c = c+1;
    end
    % Combine single A and coeffs
    allSingleA(count,:) = singleA;
    allCoeff(count,:) = singleCoeff;
    count = count+1;
end
fullFeat = names({'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
mSingleA = mean(allSingleA,1);
%% Plot
figure
hold on
plot(mean(x,1),mean(y,1),'-k')
plot(mean(xP,1),mean(yP,1),'--k')
plot(squeeze(mean(subX(1,:,:),2)),squeeze(mean(subY(1,:,:),2)))
plot(squeeze(mean(subX(2,:,:),2)),squeeze(mean(subY(2,:,:),2)))
plot(squeeze(mean(subX(3,:,:),2)),squeeze(mean(subY(3,:,:),2)))
legend({['all: ',num2str(round(mean(allA),2)),'\pm',...
    num2str(round(conf(allA,0.95),2))],['perm: ',...
    num2str(round(mean(allAP),2)),'\pm',...
    num2str(round(conf(allA,0.95),2))],['sIL: ',...
    num2str(round(mean(subA(:,1)),2)),'\pm',...
    num2str(round(conf(subA(:,1)',0.95),2))],['sNAC: ',...
    num2str(round(mean(subA(:,2)),2)),'\pm',...
    num2str(round(conf(subA(:,2)',0.95),2))],['mph: ',...
    num2str(round(mean(subA(:,3)),2)),'\pm',...
    num2str(round(conf(subA(:,3)',0.95),2))]})
%% Inter effect combined groups
for ii = 1:100
    load(['F:\irdmRound2\interEffectPosCmbGroup\interEffectPosCmbGroup',...
        num2str(ii),'.mat'],'accLOO','groupAcc')
    for jj = 1:15
        allA(ii,jj) = accLOO{jj}{1}.acc;
        allX(ii,jj,:) = interp1(linspace(0,1,numel(accLOO{jj}{1}.x)),...
            accLOO{jj}{1}.x,linspace(0,1,50));
        allY(ii,jj,:) = interp1(linspace(0,1,numel(accLOO{jj}{1}.y)),...
            accLOO{jj}{1}.y,linspace(0,1,50));
    end
    for jj = 1:3
        groupA(ii,jj) = groupAcc{ii,jj}{1}.acc;
        groupX(ii,jj,:) = interp1(linspace(0,1,numel(groupAcc{ii,jj}{1}.x)),...
            groupAcc{ii,jj}{1}.x,linspace(0,1,50));
        groupY(ii,jj,:) = interp1(linspace(0,1,numel(groupAcc{ii,jj}{1}.y)),...
            groupAcc{ii,jj}{1}.y,linspace(0,1,50));
        if jj == 1
            prob = cvglmnetPredict(groupAcc{ii,jj}{1}.mdl{1},...
                groupHist{jj}.cfg.naive.testX);
            [x,y,~,permA(ii,jj)] = perfcurve(...
                groupHist{jj}.cfg.naive.testY(randperm(numel(groupHist{jj}.cfg.naive.testY))),prob,1);
            permX(ii,:) = interp1(linspace(0,1,numel(x)),x,...
                linspace(0,1,50));
            permY(ii,:) = interp1(linspace(0,1,numel(y)),y,...
                linspace(0,1,50));
        end
    end
    groupSingleA(ii,:,:) = singleGroupA;
    groupSingleCoeff(ii,:,:) = singleGroupCoeff;
end
%%
figure
scatter3(squeeze(mean(groupSingleA(:,:,1),1)),...
    squeeze(mean(groupSingleA(:,:,2),1)),...
    squeeze(mean(groupSingleA(:,:,3),1)),'.k')
xlabel('sIL+sCore');ylabel('sIL+MPH');zlabel('sCore+MPH')

figure
hold on
plot(squeeze(mean(allX,[1,2])),squeeze(mean(allY,[1,2])),'k-')
plot(squeeze(mean(groupX(:,1,:),1)),squeeze(mean(groupY(:,1,:),1)))
plot(squeeze(mean(groupX(:,2,:),1)),squeeze(mean(groupY(:,2,:),1)))
plot(squeeze(mean(groupX(:,3,:),1)),squeeze(mean(groupY(:,3,:),1)))
legend({'all groups; LOO','sIL+sCore; 80/20','sIL+MPH; 80/20',...
    'sCore+MPH; 80/20'})
%% 
fullFeat = names({'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
[sortStim,sortInd] = sort(squeeze(mean(singleGroupA(:,:,1),1)),'descend');
sortFeat = fullFeat(sortInd);
dir = sign(mean(singleGroupCoeff(:,:,1),1));

figure
hold on
plot(squeeze(mean(groupX(:,1,:),1)),squeeze(mean(groupY(:,1,:),1)),'-k')
plot(mean(permX,1),mean(permY,1),'--k')
legend({['80/20: ',num2str(round(mean(groupA(:,1)),2)),'\pm',...
    num2str(round(conf(groupA(:,1)',0.95),2))],['perm: ',...
    num2str(round(mean(permA),2)),'\pm',...
    num2str(round(conf(permA',0.95),2))]})
box off
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
xlabel('FPR'); ylabel('TPR')
title('sIL or sNAcC')
%% sIL + sCore; male and female
load('F:/irdmRound2/maleFemaleILCore.mat')
for ii = 1:100
    mA(ii) = maleAcc{ii}{1}.acc;
    mX(ii,:) = interp1(linspace(0,1,numel(maleAcc{ii}{1}.x)),...
        maleAcc{ii}{1}.x,linspace(0,1,50));
    mY(ii,:) = interp1(linspace(0,1,numel(maleAcc{ii}{1}.y)),...
        maleAcc{ii}{1}.y,linspace(0,1,50));
    fA(ii) = femaleAcc{ii}{1}.acc;
    fX(ii,:) = interp1(linspace(0,1,numel(femaleAcc{ii}{1}.x)),...
        femaleAcc{ii}{1}.x,linspace(0,1,50));
    fY(ii,:) = interp1(linspace(0,1,numel(femaleAcc{ii}{1}.y)),...
       femaleAcc{ii}{1}.y,linspace(0,1,50));
end
mS = mean(maleA,1);
fM = mean(femaleA,1);
mC = mean(maleCoeff,1);
mCS = sign(mC);
fS = mean(femaleA,1);
fC = mean(femaleCoeff,1);
fCS = sign(fC);
feat = names({'lmPFC','rmPFC','lOFC','rOFC','lNAcS','rNAcS','lNAcC',...
    'rNAcC'},{'d','t','a','b','lg','hg'});
[~,ind] = sort(mS,'descend');
mSort = mS(ind)'.*mCS(ind)';
mFeat = feat(ind)';
[~,ind] = sort(fS,'descend');
fSort = fS(ind)'.*fCS(ind)';
fFeat = feat(ind)';
%%
figure
hold on
plot(mean(mX,1),mean(mY,1))
plot(mean(fX,1),mean(fY,1))
legend({['male: ',num2str(round(mean(mA),2)),'\pm',...
    num2str(round(conf(mA,0.95),2))],['female: ',...
    num2str(round(mean(fA),2)),'\pm',num2str(round(conf(fA,0.95),2))]})
set(gca,'xtick',0:0.5:1,'ytick',0:0.5:1)
c = distinguishable_colors(6);
figure
hold on
scatter(mS(1:48).*mCS(1:48),fS(1:48).*fCS(1:48),[],repmat(c,8,1),'o');
scatter(mS(49:end).*mCS(49:end),fS(49:end).*fCS(49:end),[],repmat(c,28,1),'s');
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
xlabel('male'); ylabel('female')
%% Get total samples for each effect for each intervention for each animal
effects = [1,-1,0];
for ii = 1:numel(inter)
    interInds = logicFind(1,contains(interTable.Rx,inter{ii}),'==');
    for jj = 1:3
        effectInds = interInds(interTable.effect(interInds) == effects(jj));
        uID = unique(interTable.ID(effectInds));
        for k = 1:numel(uID)
            idInds = effectInds(contains(interTable.ID(effectInds),uID{k}));
            samples{ii,jj}{k} = interTable.samps(idInds);
        end
    end
end
%% Plot auc histograms separated by effect
figure
histogram(interTable.auc)
histogram(interTable.auc(interTable.effect==0))
hold on
histogram(interTable.auc(interTable.effect==1))
histogram(interTable.auc(interTable.effect==-1))
figure
histogram(interTable.auc(interTable.effect==0)-interTable.b_auc(interTable.effect==0))
hold on
histogram(interTable.auc(interTable.effect==1)-interTable.b_auc(interTable.effect==1))
histogram(interTable.auc(interTable.effect==-1)-interTable.b_auc(interTable.effect==-1))
%% Figures
nAnimal = size(nSampInt,1);
% plain stims
figure
for jj = 1:2
    subplot(2,1,jj)
    hold on
    for ii = 1:nAnimal
        if ~isempty(nSampInt{ii,jj})
            plot(ones(1,numel(nSampInt{ii,jj})).*ii,nSampInt{ii,jj},'.',...
                'color',[0.5 0.5 0.5],'markersize',15)
        end
        plot(ii,mean(nSampInt{ii,jj}),'rs','markerfacecolor','r')
    end
    set(gca,'xtick',1:nAnimal,'xticklabel',uID,'xlim',[0 nAnimal+1])
    xtickangle(90)
    title(inter{jj})
end
% guan
figure
for jj = 3
    hold on
    for ii = 1:nAnimal
        if ~isempty(nSampInt{ii,jj})
            plot(ones(1,numel(nSampInt{ii,jj})).*ii,nSampInt{ii,jj},'.',...
                'color',[0.5 0.5 0.5],'markersize',15)
        end
        plot(ii,mean(nSampInt{ii,jj}),'rs','markerfacecolor','r')
    end
    set(gca,'xtick',1:nAnimal,'xticklabel',uID,'xlim',[0 nAnimal+1])
    xtickangle(90)
    title(inter{jj})
end
% mph
figure
c = 1;
for jj = 4:7
    subplot(2,2,c)
    hold on
    for ii = 1:nAnimal
        if ~isempty(nSampInt{ii,jj})
            plot(ones(1,numel(nSampInt{ii,jj})).*ii,nSampInt{ii,jj},'.',...
                'color',[0.5 0.5 0.5],'markersize',15)
        end
        plot(ii,mean(nSampInt{ii,jj}),'rs','markerfacecolor','r')
    end
    set(gca,'xtick',1:nAnimal,'xticklabel',uID,'xlim',[0 nAnimal+1])
    xtickangle(90)
    title(inter{jj})
    c = c+1;
end
% stim-LSD
figure
c = 1;
for jj = 8:9
    subplot(1,2,c)
    hold on
    for ii = 1:nAnimal
        if ~isempty(nSampInt{ii,jj})
            plot(ones(1,numel(nSampInt{ii,jj})).*ii,nSampInt{ii,jj},'.',...
                'color',[0.5 0.5 0.5],'markersize',15)
        end
        plot(ii,mean(nSampInt{ii,jj}),'rs','markerfacecolor','r')
    end
    set(gca,'xtick',1:nAnimal,'xticklabel',uID,'xlim',[0 nAnimal+1])
    xtickangle(90)
    title(inter{jj})
    c = c+1;
end
%% Plot histograms of samples depending on effect
figure
hold on
histogram(interTable.samps(interTable.effect==1))
histogram(interTable.samps(interTable.effect==-1))
histogram(interTable.samps(interTable.effect==0))
%% Unpack effectSample
posSamp = cell(1,9);
negSamp = cell(1,9);
nullSamp = cell(1,9);
for ii = 1:size(effectSample,1)
    for jj = 1:size(effectSample,2)
        % Skip if no null effect
        if effectSample{ii,jj}(3) ~= 0
            % Pull out positive effect
            if effectSample{ii,jj}(1) ~= 0
                posSamp{jj} = [posSamp{jj};effectSample{ii,jj}([1,3])];
            end
            % Pull out negative effect
            if effectSample{ii,jj}(2) ~= 0
                negSamp{jj} = [negSamp{jj};effectSample{ii,jj}([2,3])];
            end
        end
    end
end
req = 10:1000;
for jj = 1:9
    c = 1;
    for ii = req
        % Account for needing null samples as well by checking both
        pos(jj,c) = sum(sum(posSamp{jj}>=ii,2)==2);
        neg(jj,c) = sum(sum(negSamp{jj}>=ii,2)==2);
        c = c+1;
    end
end
% Get degree of imbalance in data
for jj = 1:9
    %     if
end
% Plot
figure
for ii = 1:2
    subplot(2,1,ii)
    plot(req,pos(ii,:))
    hold on
    plot(req,neg(ii,:))
    title(inter{ii})
end

figure
plot(req,pos(3,:))
hold on
plot(req,neg(3,:))
title(inter{3})

figure
c = 1;
for ii = 4:7
    subplot(2,2,c)
    plot(req,pos(ii,:))
    hold on
    plot(req,neg(ii,:))
    title(inter{ii})
    c = c+1;
end

figure
c = 1;
for ii = 8:9
    subplot(2,1,c)
    plot(req,pos(ii,:))
    hold on
    plot(req,neg(ii,:))
    title(inter{ii})
    c = c+1;
end
%% Build models from sIL and sNAc
inter = {'sIL','sNAc'};
n = 3;
for ii = 1:numel(inter)
    for jj = 1:numel(uID)
        for k = [0,1,-1]
            % Group files by, intervention, ID, and effect of intervention
            int = contains(interTable.Rx,inter{ii});
            animal = contains(interTable.ID,uID{jj});
            group = interTable.effect == k;
            theseInds = logicFind(1,int & animal & group,'==');
            % Open each file; check for enough samples during intervention
            for m = theseInds
                load(['F:/irdmRound2/processedContinuous3/',...
                    interTable.ID{m},'-',inter{ii},'_DDT',...
                    interTable.Date{m},'_all.mat'])
                % Get first stim time
                first = min([hist.eventTs.t{9}(1),hist.eventTs.t{10}(1)]);
                last = max([hist.eventTs.t{9}(end),...
                    hist.eventTs.t{10}(end)]);
                % Grab any data that comes from between stim (during actual
                % stim should already be removed)
                notNaN = squeeze(~isnan(psdTrls{1}.bandPow(1,1,1,:)))';
                % Use -1.5 to account for window, and +60 to account for
                % washout, and notNaN
                interInd = psdTrls{1}.t > first-1.5 & psdTrls{1}.t < last+60 & notNaN;
                thisStim = [reshape(psdTrls{1}.bandPow(:,:,1,interInd),...
                    36,sum(interInd))',...
                    reshape(coh{1}.mBandCoh(:,:,1,interInd),90,...
                    sum(interInd))'];
                % If there is any, grab pre and post, stim as well
                preData = [];
                postdata = [];
                % Then use n baselines to also load and use as
                % distribution to calculate z-scores from
                for b = 1:n
                    
                end
            end
        end
    end
end
%%
for n = 53:100
c = 3;
count = 1;
for interI = 4%[1,2,4]
    allPosX = [];
    allPosY = [];
    allNegX = [];
    allNegY = [];
    % Positive model
    for animI = logicFind(1,posInds(interI,:),'==')
%         posX = [null{interI,animI}(nullSamp{c,animI}(n,:),:);...
%             pos{interI,animI}(posSamp{c,animI}(n,:),:)];
%         posY = [zeros(size(posX,1)/2,1);ones(size(posX,1)/2,1)];
        posX = [null{interI,animI};pos{interI,animI}];
        posY = [zeros(size(null{interI,animI},1),1);...
            ones(size(pos{interI,animI},1),1)];
        % Add to group data - only if no NaNs
        if ~any(any(isnan(posX(:,:))))
            allPosX = [allPosX;posX];
            allPosY = [allPosY;posY];
        end
        % Check for and remove NaN
        posX(:,any(isnan(posX))) = [];
        testY = [];
        % Generate test and train sets, s.t. testY has both 0s and 1s
        while ~(any(testY==0) && any(testY==1))
            [trainX,trainY,testX,testY] = trainTest(posX,posY,0.2);
        end
        % Lasso
        cfg = lassoNetCfg({testX,testY},[],'n','n','n',100,'1se',[]);
        [~,lambda,beta,fits,acc,hist] = lassoNet(trainX,...
            trainY,'binomial','class',1,5,1,cfg);
        posAcc{c,animI} = acc;
        posHist{c,animI} = hist;
        % GLM
        mdl = fitglm(trainX,trainY,'distribution','binomial');
        p = predict(mdl,testX);
        [~,~,~,posGLMAcc(c,animI)] = perfcurve(testY,p,1);
        % Lasso permuted
        cfg = lassoNetCfg({testX,testY},[],'y','n','n',100,'1se',[]);
        [~,lambdaR,betaR,fitsR,accR,histR] = lassoNet(trainX,...
            trainY,'binomial','class',1,5,1,cfg);
        posAccR{c,animI} = accR;
        % Single feature
        for jj = 1:size(trainX,2)
            mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
            p = predict(mdl,testX(:,jj));
            [xPos{c,animI,jj},yPos{c,animI,jj},~,aPos(c,animI,jj)] =...
                perfcurve(testY,p,1);
            coeffPos(c,animI,jj) = table2array(mdl.Coefficients(2,1));
        end
        % Also build group model with this animal left out; unless NaNs
        if ~any(any(isnan(...
                null{interI,animI}(nullSamp{c,animI}(n,:),:))))...
                && ~any(any(isnan(...
                pos{interI,animI}(posSamp{c,animI}(n,:),:))))
            possible = logicFind(1,posInds(interI,:),'==');
            others = possible(~ismember(possible,animI));
            theseNull = [];
            thesePos = [];
            for oi = others
                thisNull = null{interI,oi}(nullSamp{c,oi}(n,:),:);
                thisPos = pos{interI,oi}(posSamp{c,oi}(n,:),:);
                if ~any(any(isnan(thisNull))) && ~any(any(isnan(thisPos)))
                    theseNull = [theseNull;...
                        null{interI,oi}(nullSamp{c,oi}(n,:),:)];
                    thesePos = [thesePos;...
                        pos{interI,oi}(posSamp{c,oi}(n,:),:)];
                end
            end
            trainX = [theseNull;thesePos];
            trainY = [zeros(size(theseNull,1),1);ones(size(thesePos,1),1)];
            testX = [null{interI,animI}(nullSamp{c,animI}(n,:),:);...
                pos{interI,animI}(posSamp{c,animI}(n,:),:)];
            testY = [zeros(size(...
                null{interI,animI}(nullSamp{c,animI}(n,:),:),1),1);...
                ones(size(pos{interI,animI}(posSamp{c,animI}(n,:),:),1),1)];
            if isempty(trainX)
                posLOOAcc{c,animI} = [];
                posLOOHist{c,animI} = [];
                posLOOA(c,animI) = NaN;
            else
                % Lasso
                cfg = lassoNetCfg({testX,testY},[],'n','n','n',100,'1se',[]);
                [~,lambda,beta,fits,acc,hist] = lassoNet(trainX,...
                    trainY,'binomial','class',1,5,1,cfg);
                posLOOAcc{c,animI} = acc;
                posLOOHist{c,animI} = hist;
                % GLM
                mdl = fitglm(trainX,trainY,'distribution','binomial');
                p = predict(mdl,testX);
                [~,~,~,posLOOGLMAcc(c,animI)] = perfcurve(testY,p,1);
            end
        end
        save(['F:/irdmRound2/interEffect/interEffect',num2str(n),'_',num2str(count),'.mat']...
            ,'posAcc','posAccR','posLOOAcc','posLOOHist',...
            'posGLMAcc','xPos','yPos','aPos',...
            'posHist','posLOOGLMAcc','coeffPos')
        count = count+1;
    end
    %%
    % Build group model with 80/20 test
    [trainX,trainY,testX,testY] = trainTest(allPosX,allPosY,0.2);
    % Lasso
    cfg = lassoNetCfg({testX,testY},[],'n','n','n',100,'1se',[]);
    [~,lambda,beta,fits,acc,hist] = lassoNet(trainX,...
            trainY,'binomial','deviance',1,10,1,cfg);
    allPosAcc{c} = acc;
    allPosHist{c} = hist;
    % GLM
    mdl = fitglm(trainX,trainY,'distribution','binomial');
    p = predict(mdl,testX);
    [~,~,~,allPosGLMAcc(c)] = perfcurve(testY,p,1);
  save(['F:/irdmRound2/interEffect/interEffectGroup',num2str(n),'.mat'],'allPosAcc',...
        'allPosHist','allPosGLMAcc')
    %% Only build negative model for sIL
    if interI == 1
        % Negative model
        for animI = logicFind(1,negInds(interI,:),'==')
%             negX = [null{interI,animI}(nullSamp{c,animI}(n,:),:);...
%                 neg{interI,animI}(negSamp{c,animI}(n,:),:)];
%             negY = [zeros(size(negX,1)/2,1);ones(size(negX,1)/2,1)];
            negX = [null{interI,animI};neg{interI,animI}];
            negY = [zeros(size(null{interI,animI},1),1);...
                ones(size(neg{interI,animI},1),1)];
                    % Add to group data - only if no NaNs
            if ~any(any(isnan(negX(:,:))))
                allNegX = [allNegX;negX];
                allNegY = [allNegY;negY];
            end
            negX(:,any(isnan(negX))) = [];
            testY = [];
            % Generate test and train sets, s.t. testY has both 0s and 1s
            count = 1;
            while ~(any(testY==0) && any(testY==1))
                [trainX,trainY,testX,testY] = trainTest(negX,negY,0.2);
                count = count+1;
                if count == 1000
                   warning('danger! you might be stuck here forever') 
                end
            end
            % Lasso
            cfg = lassoNetCfg({testX,testY},[],'n','n','n',100,'1se',[]);
            [~,lambda,beta,fits,acc,hist] = lassoNet(trainX,...
                trainY,'binomial','class',1,5,1,cfg);
            negAcc{c,animI} = acc;
            negHist{c,animI} = hist;
            % GLM
            mdl = fitglm(trainX,trainY,'distribution','binomial');
            p = predict(mdl,testX);
            [~,~,~,negGLMAcc(c,animI)] = perfcurve(testY,p,1);
%             Lasso permuted
            cfg = lassoNetCfg({testX,testY},[],'y','n','n',100,'1se',[]);
            [~,lambdaR,betaR,fitsR,accR,histR] = lassoNet(trainX,...
                trainY,'binomial','class',1,5,1,cfg);
            negAccR{c,animI} = accR;
            
            % Single feature
            for jj = 1:size(trainX,2)
                mdl = fitglm(trainX(:,jj),trainY,'distribution','binomial');
                p = predict(mdl,testX(:,jj));
                [xNeg{c,animI,jj},yNeg{c,animI,jj},~,aNeg(c,animI,jj)] =...
                    perfcurve(testY,p,1);
                coeffNeg(c,animI,jj) = table2array(mdl.Coefficients(2,1));
            end
            % Also build group model with this animal left out; unless
            % lonely
            if numel(logicFind(1,negInds(interI,:),'==')) ~= 1
                possible = logicFind(1,negInds(interI,:),'==');
                others = possible(~ismember(possible,animI));
                theseNull = [];
                theseNeg = [];
                for oi = others
                    theseNull = [theseNull;...
                        null{interI,oi}(nullSamp{c,oi}(n,:),:)];
                    theseNeg = [theseNeg;...
                        neg{interI,oi}(negSamp{c,oi}(n,:),:)];
                end
                trainX = [theseNull;theseNeg];
                trainY = [zeros(size(theseNull,1),1);...
                    ones(size(theseNeg,1),1)];
                testX = [null{interI,animI}(nullSamp{c,animI}(n,:),:);...
                    neg{interI,animI}(negSamp{c,animI}(n,:),:)];
                testY = [zeros(size(testX,1)/2,1);ones(size(testX,1)/2,1)];
                % Lasso          
                cfg = lassoNetCfg({testX,testY},[],'n','n','n',100,'1se',[]);
                [~,lambda,beta,fits,acc,hist] = lassoNet(trainX,...
                    trainY,'binomial','class',1,5,1,cfg);
                negLOOAcc{c,animI} = acc;
                negLOOHist{c,animI} = hist;
                % GLM
                mdl = fitglm(trainX,trainY,'distribution','binomial');
                p = predict(mdl,testX);
                [~,~,~,negLOOGLMAcc(c,animI)] = perfcurve(testY,p,1);
            else
                negLOOACC{c,animI} = [];
                negLOOGLMAcc(c,animI) = NaN;
            end
        end
        % Build group model with 80/20 test
        if ~isempty(allNegX)
            [trainX,trainY,testX,testY] = trainTest(allNegX,allNegY,0.2);
            % Lasso
            cfg = lassoNetCfg({testX,testY},[],'n','n','n',100,'1se',[]);
            [~,lambda,beta,fits,acc,hist] = lassoNet(trainX,...
                trainY,'binomial','class',1,5,1,cfg);
            allNegAcc{c} = acc;
            allNegHist{c} = hist;
            % GLM
            mdl = fitglm(trainX,trainY,'distribution','binomial');
            p = predict(mdl,testX);
            [~,~,~,allNegGLMAcc(c,animI)] = perfcurve(testY,p,1);
        end
    end
    c = c+1;
end
end