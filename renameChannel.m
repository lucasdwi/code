function renameChannel(sdir,searchStr,names)
%% Goes through files within source directory and renames channels. 
% WARNING: FUNCTION WILL OVERWRITE FILES!
%__________________________________________________________________________
% INPUTS
% sdir = source directory, path of directory in which to look for files;
%   format = string
% searchStr = string to look for in fileneames; format = string
% names = new channel names; format = cell row vector; N.B. MUST MATCH
%   DIMENSIONS OF OLD LFPTs.label, ELSE WILL GIVE WARNING AND SKIP FILE
%__________________________________________________________________________
% LLD 2018
%%
files = [];
cd(sdir)
thisF = dir(strcat('*',searchStr,'*'));
% Just pull out 'name' field and concatenate to 'files'
files = vertcat(files,extractfield(thisF,'name')');
for fi = 1:size(files,1)
    disp(['Loading file ',num2str(fi), ' of ',num2str(size(files,1))])
    load(files{fi})
    % Check for matching dimensions
    if size(LFPTs.label) == size(names)
    % Rename channel labels
    LFPTs.label = names;
    disp('Saving...')
    % Save all variables except for those used internally
    save(files{fi},'-regexp',...
        '^(?!(thisF|files|sdir|searchStr|names|fi)$).')
    else
        warning(['Dimensions of new and old channel name cell array do '...
            'not match. Skipping file!'])
    end
end