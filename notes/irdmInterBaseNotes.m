load('F:\irdmRound2\ddtFiles.mat')
% Remove files with: 'skip','neural','seizure'
skipInd = contains(processedData(:,6),'skip') | ...
    contains(processedData(:,6),'neural') | ...
    contains(processedData(:,6),'seizure') | ...
    contains(processedData(:,6),'OFC');
processedData(skipInd,:) = [];

% Get list of all processed files
files = fileSearch('F:/irdmRound2/processedContinuous3/','.mat');
% Set n, numer of baselines to use
n = 3;
% Set up indices of pow of coh
powInds = reshape(1:48,6,8)';
cohInds = reshape(49:216,6,28)';
% Get unique IDs
uID = unique(processedData(:,1));
% Go through each animal and find the recordings with interventions and
% their prior baselines; also calculate number of clean trials in data
inter = {'sIL','sNAcC','guan0.3','mph3','mph0.3','mph1','mph0.1',...
    'sNAcC-LSD','sIL-LSD'};
% Clear table
interTable = array2table(zeros(0,10));
interBaseInd = cell(numel(uID),numel(inter));
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
    % Figure out which is the first intervention by using indices (sorted
    % by date)
    firstInt = min(cell2mat(cellfun(@min,intInd(ii,:),...
        'UniformOutput',false)));
    % For each intervention, find the n closest baselines
    interBaseInd = cell(numel(inter),max(cell2mat(cellfun(@numel,...
        intInd(ii,:),'UniformOutput',0))));
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
                    datetime(processedData(theseBaseInd(end-(n-1):end),...
                    2)));
                % Store info, including samples and nDays
                interBase{ii,jj,k} = [processedData(theseBaseInd(end-...
                    (n-1):end),:),...
                    num2cell(nDays)];
                interBaseInd{jj,k} = num2str(theseBaseInd(end-(n-1):end));
            end
        end
    end
    lastBaseInd = logicFind(1,cell2mat(cellfun(@(x) any(x==firstInt),...
        intInd(ii,:),'UniformOutput',false)),'==');
    if ~isempty(interBaseInd)
        these = str2num(interBaseInd{lastBaseInd,1});
        for k = 1:n
            % Find file using animal and date
            thisFile = files(contains(files,...
                processedData{these(k),1}) & ...
                contains(files,processedData{these(k),2}));
            % Load and remove NaNed data
            load(['F:/irdmRound2/processedContinuous3/',thisFile{1}]);
            thisNotNan = ~isnan(squeeze(psdTrls{1}.relPow(1,1,1,:))) &...
                ~isnan(squeeze(coh{1}.normBandCoh(1,1,1,:)));
            notNanPow = squeeze(psdTrls{1}.relPow(:,:,1,thisNotNan));
            notNanCoh = squeeze(coh{1}.normBandCoh(:,:,1,thisNotNan));
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
            thisBase{ii,k} = [rePow,reCoh];
        end
    end
end
save('lastBase.mat','thisBase','uID')
%% Load models
for ii = 1:100
    load(['F:/irdmRound2/baseEffect/baseEffect',num2str(ii),'.mat'])
    ilA(ii) = multiClassAUC(ILacc{1}.pred{1},ILhist.cfg.naive.testY);
    coreA(ii) = coreAcc{1}.acc;
end