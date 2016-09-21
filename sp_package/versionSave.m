function [fileVersion] = versionSave(saveParent,name)
%% Checks if saveDir already exists and if so, renames it
% INPUTS:
% saveDir = directory to save file; format = string


% OUTPUT:
% fileVersion = new name to use for saving
%% Go to data location
cd(saveParent)
%% Check if directory exists, if not make it
if exist(name,'dir') == 0
    mkdir(name);
% Otherwise count number of folders with same name and make new directory
else
    files = dir([name,'*']);
    dirFlags = [files.isdir];
    fDir = files(dirFlags);
    fNum = size(fDir,1);
    mkdir(strcat(name,num2str(fNum)));
end
%% Go into saveDir and search for already existing file matches
cd(saveDir);
fDir = dir([fileName,'*']);
fNum = size(fDir,1) + 1;
fileVersion = strcat(fileName,num2str(fNum))