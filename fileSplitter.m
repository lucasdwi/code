function fileSplitter(sDir,searchStr,nChan)
%% Splits LFPTs and eventTs
%__________________________________________________________________________
% INPUTS
% sDir = source directory; where you want the program to look for files;
%    format = string, preferably ending with \
% searchStr = string you want to use to find files; format = string
% chan = number of channels (LFP data) you expect each animal/box to have
%__________________________________________________________________________
% EXAMPLE
% fileSplitter('F:\DDT\pl2\','.mat',8)
% looks for all files that contain '.mat' (easy way to grab all matlab
% files) that exist in 'F:\DDT\pl2\ and splits them with 8 channels in each
% resulting file
%__________________________________________________________________________
% N.B.: 
% - digital event channels (box or stim pulses) are hardcoded at beginning
% - assumes date is last part of name
% - uses numbers in file name to count number of animals in recording 
%   (if a number is used within the name, they are fine if included with 
%   hyphens, e.g., IRDM23-MPH3)
%%
% check sDir ends in \ or /; else add \
if ~contains(sDir(end),'\')
    sDir = [sDir,'\'];
end
files = fileSearch(sDir,searchStr);
% Set up box indices
boxInds = [25:32;33:40;41:48;49:56];
stimInds = reshape(9:24,4,4)';
% dInds = reshape(1:32,8,4)';
chans = cell(1,32);
for ii = 1:32
    if ii < 10
        chans(ii) = {['FP0',num2str(ii)]};
    else
        chans(ii) = {['FP',num2str(ii)]};
    end
end
chans = reshape(chans,8,4)';
for fI = 1:size(files,2)
    disp(['Splitting ',files{fI},': ',num2str(fI),' of ',...
        num2str(size(files,2))])
    load(files{fI})
    % Grab this filename
    thisFile = files{fI};
    % Remove file extension by finding period; '.'
    thisFile = thisFile(1:strfind(thisFile,'.')-1);
    % Split str at all underscores; '_'
    parts = strsplit(thisFile,'_');
    % Find all parts that are numbers
    nums = regexp(parts,'\d+','match');
    % Find all number parts excluding the date
    numInd = logicFind(1,~cellfun(@isempty,nums(1:end-1)),'==');
    % Check that the number of sub-files and expected channels per sub-file
    % matches the total amount of available data
    if ~isequal(size(numInd,2)*nChan,size(LFPTs.data,1))
        warning(['There are a different number of channels than expected'... 
            ' by the number of animals.'])
    end
    % Store LFPTs and eventTs
    oldLFPTs = LFPTs;
    oldEventTs = eventTs;
    % Set rat counter and default number of boxes (4)
    nRat = 1;
    nBox = 4;
    % Split, going through each potential box
    for k = 1:nBox
        % First check if any data from this box exists and get inds
        dInds = contains(oldLFPTs.label,chans(k,:));
        if any(dInds)
            % Split events
            theseInds = [boxInds(k,:),stimInds(k,:)];
            eventTs.t = oldEventTs.t(theseInds);
            eventTs.label = oldEventTs.label(theseInds);
            % Split LFP data
            LFPTs.data = oldLFPTs.data(dInds,:);
            LFPTs.label = oldLFPTs.label((dInds));
            % Create new file name
            newName = strjoin([parts{nRat},parts(end-1:end)],'_');
%             % Run ddtTrials
%             if ddt
%                 % Check for empty feeder (indicates)
%                 trials = ddtTrials(eventTs,LFPTs,0);
%                 % Add auc to table
%                 auc = [auc;{newName},trials(1,1).auc];
%             else
%                 trials = [];
%             end
%             save(['F:\irdmRound2\toProcess\',newName],'eventTs','LFPTs',...
%                 'adfreq','trials')
            if ~exist([sDir,'split'], 'dir')
                mkdir([sDir,'split'])
            end
            save([sDir,'split\',newName],'LFPTs','eventTs','adfreq')
            nRat = nRat +1;
        else
            warning(['LFP missing from box ',num2str(k)])
        end
    end
%     movefile(files{ii},'F:\irdmRound2\mat\')
end