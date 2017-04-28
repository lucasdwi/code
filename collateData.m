function [data,samp] = collateData(sdir,searchStr,vars,dat)
%% Collates data from many files into matrices.
%__________________________________________________________________________
% INPUTS:
% sdir = file source directory; format: string
% searchStr = cell array of strings to search for in file titles; format:
%   cell array with first column indicating string and second column
%   indicating criteria ('in' or 'ex')
% vars = variables to collate; format: string; options: 'pow','coh','corr'
% dat = kind of data to be collated indicating whether all trials are used
%   or just the averages; format: string, either 'avg' or 'trl'
%__________________________________________________________________________
% OUTPUTS:
% data = data cell array of matrices for each searchStr used
%__________________________________________________________________________
% USE: 
% [data] = collateData('C:\Users\matFiles\',{'foo';'bar'},{'pow','coh'},'trl') 
% Will go through ~\matFiles\ and search for any files with 'foo' in the title
% and 'bar' in the title. Then will open these files concatenating all of
% the trials for each behavior into cell arrays. So if there are two
% behaviors per file, data will be a (1,2) cell (one cell for 'foo' files
% and for 'bar' files). Within each of these cells will be a (n,2) cell
% where n is the number of files with 'foo' or 'bar' in title and the
% column dimension holds the two behaviors.
%__________________________________________________________________________
% DEPENDENCIES:
% fileSearch.m
% logicFind.m
%__________________________________________________________________________
% LLD 2017
%% Initialization
% Get number of strings going to be used in searchStr
nStr = size(searchStr,1);
% Check that dat is set to either 'trl' or 'avg'
if ~strcmpi(dat,'avg') && ~strcmpi(dat,'trl')
    error('Warning: dat needs to be set to either "avg" or "trl".')
end
% Checks that if dat is set to trl, then 'corr' is not in vars
if strcmpi(dat,'trl') && sum(strcmpi(vars,'corr')) == 1
    error('Warning: Power correlations can not be used for trialized data.')
end
% Checks if searchStr has a criterion per string (stored in second column);
% if not, assumes inclusion and appends column of 'in'
if size(searchStr,2) == 1
    searchStr = [searchStr,repmat({'in'},size(searchStr,1),1)];
end
% Preallocate files
if nStr > 1
    files = cell(nStr,1);
else
    files = cell(1,1);
end
%% Cycle through each searchStr and get file names
for sI = 1:size(searchStr,1)
    [files{sI}] = fileSearch(sdir,searchStr{sI,1},searchStr{sI,2});
end
nAllFile = sum(cellfun(@numel,files));
%% Go through each cell of files and collate data from files within into
% one structure per cell of files
% Set up file counter
count = 1;
% Set up data
data = cell(1,nStr);
% Set up samp
samp = cell(1,nStr);
for sI = 1:nStr
    % Get number of files within given cell
    nFile = size(files{sI},2);
    for fI = 1:nFile
        disp(['Adding file ',num2str(count),' of ',num2str(nAllFile)])
        load([sdir,files{sI}{fI}]);
        % Add one to counter
        count = count + 1;
        for iE = 1:size(trls,2)
            thisData{iE} = [];
            % Normalized power - reshapes into row vector with columns of
            % the following pattern: c1b1,c1b2,c1b3,c1b4,c2b1,...
            if sum(strcmpi(vars,'pow')) == 1
                if strcmpi(dat,'trl')
                    [b,c,t] = size(psdTrls{iE}.relPow);
                    % Reshape and transpose
                    thisPow = reshape(psdTrls{iE}.relPow,b*c,t)';
                else
                    [b,c] = size(psdTrls{iE}.avgRelPow);
                    % Reshape
                    thisPow = reshape(psdTrls{iE}.avgRelPow,1,b*c);
                end
                % Store power
                thisData{iE} = [thisData{iE},thisPow]; %#ok<*AGROW>
            end
            % Normalized coherence - reshapes into row vector with columns
            % of the following pattern:
            % c1c2b1,c1c2b2,c1c2b3,c1c2b4,c1c3b1...
            if sum(strcmpi(vars,'coh')) == 1
                if strcmpi(dat,'trl')
                    [cmb,b,t] = size(coh{iE}.rel);
                    % Permute coh.rel into a similar pattern as power
                    thisCoh = reshape(permute(coh{iE}.rel,[2,1,3]),cmb*b,t)';
                else
                    [cmb,b,~] = size(coh{iE}.rel);
                    % Mean and permute coh.rel
                    thisCoh = reshape(permute(mean(coh{iE}.rel,3),[2,1]),1,cmb*b);
                end
                % Store coherence
                thisData{iE} = [thisData{iE},thisCoh];
            end
            % Power correlation [only available for avg data sets] -
            % reshapes into following pattern:
            % c1c2b1,c1c2b2,c1c2b3,c1c2b4,c1c3b1...
            if sum(strcmpi(vars,'corr')) == 1
                if strcmpi(dat,'avg')
                    [b,cor] = size(rVect{iE});
                    % Reshape
                    thisCorr = reshape(rVect{iE},1,b*cor);
                end
                % Store power correlations
                thisData{iE} = [thisData{iE},thisCorr];
            end
            % Grab sampleInfo for timing
            thisSamp{iE} = trls{1,iE}.sampleinfo;
        end
        data{sI} = [data{sI};thisData];
        samp{sI} = [samp{sI};thisSamp];
    end
end
