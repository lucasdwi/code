function [varargout] = fileCycle(fun,fType,files,sdir)
%% Generic script to batch process files
% Inactivate/activate function cells with commenting

% Inputs (pick one from fType and files, set the other to []):
% fun = function group to run; options = 'scb' for spectcompbase.m; 'tab'
%   for tabulateData.m 
% fType = wildcard(s) to search for in order to populate list of files to
%   be processed; format = string array
%   N.B.: not case sensitive, but use best practice to avoid problems
% files = array of filenames to use
% sdir = source directory of files; format = string
%   N.B.: if more than one source directory is to be used, rerun
%   fileCycle.m for each

% Outputs: Depend on function group called (see those functoins for more
%   detail)

% Example:
% fileCycle('scb',{'Base','FoodDep24'},[],'C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed\');
% Runs all files in C:\...\channel_renamed that have either 'Base' or
% 'FoodDep24' in their name through spectcompbase.m

% 1st batch of files: {'H10BaseSep27','H10FoodDepSep25','H13BaseSep27','H13FoodDepSep25','H14BaseSep29','H14FoodDepSep28','H15BaseSep29','H15FoodDepSep28','I1BaseNov9','I1FoodDepNov13','I2BaseNov9','I2FoodDep24Dec16','I3BaseNov9','I3FoodDep24Nov2','I4BaseSep24','I4FoodDepSep30','I6BaseSep18','I6FoodDepSep30','I8BaseSep23','I8FoodDepSep30','I11BaseOct30','I11FoodDep24Dec16','I12BaseNov12','I12FoodDepNov13'};
% 2nd batch: {'I2BaseDec15','I3BaseNov11','I4BaseDec4','I6BaseNov24','I8BaseOct29','I11BaseOct29','I12BaseOct30'}
%   'H10BaseOct15','H13BaseNov12','H14BaseOct15','H15BaseOct11','I1BaseOct26',
%% Setup files to be proceesed, then run through spectcompbase.m
if strcmp(fun,'scb')
    % First sefup files to be processed if 'fType' rather than 'files' is
    % used
    cd(sdir);
    if ~isempty(fType)
        % Creates data structure with information on files with wildcard fType
        % in name; if >1 fType, then concatentates together
        files = [];
        for f = 1:length(fType)
            thisF = dir(strcat('*',fType{f},'*'));
            % Just pull out 'name' field and concatenate to 'files'
            files = vertcat(files,extractfield(thisF,'name')');
        end
    end
    % Run spectcompbase.m
    for i = 1:length(files)
        [LFPTs,nNaN,indSkp,trls,clnTrls,clnEvents,relPower,psdTrls,TFRs,fds,avgCoh,relCoh,~,~] = spectcompbase(sdir,files{i},'y',5,2,5,17000,3,1.5,[1 2 150],0.5,{1,[0 0.005 3];2,[0 0.005 3];3,[0 0.005 3]},[3])
        close all; clearvars -except files i;
    end
end
%% Initialize files to be tabulated, then run through tabulateData.m
if strcmp(fun,'tab')
    % Initialize files     
    % Go to directory of processed data
    cd(sdir);
    if ~isempty(fType)
        % Delete 'files' so it can be remade as structure
        clear files;
        % Creates data structure with information on files with wildcard fType
        % in name; if >1 fType, then sets up different child directories
        % in 'files'
        for f = 1:length(fType)
            files.(fType{f}) = dir(strcat('*',fType{f},'*'));
        end
    end   
    %% Run tabulateData.m
    tic
    [T,allPredict,allResponse,allGroups] = tabulateData(files,[3],2);
    toc
    % Set up varargout options
    varargout{1} = T; varargout{2} = allPredict; varargout{3} = allResponse; varargout{4} = allGroups;
end
%%
% [models] = lineReg(T,allVars);
%% Get RMSs for all channels of all files
% cd('C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\channel_renamed\');
% for ii = 1:length(files)
%     tic
%     ii
%     % Load file
%     load(strcat(files{ii},'.mat'));
%     % Threshold
%     [LFPTsNaN,nNaN,indSkp] = threshFilt(LFPTs,2,5,17000,3,1.5,adfreq,1,size(LFPTs.data,1));
%     for j = 1:size(LFPTsNaN.data,1)
%         r(ii,j) = rms(LFPTsNaN.data(j,~isnan(LFPTsNaN.data(j,:))));
%     end
%     toc
% end
