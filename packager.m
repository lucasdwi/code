function [fList] = packager(func,dest)
% INPUTS 
% func = function to package; format = string corresponding to function;
%   e.g. 'packager.m'
% dest = destination; format = path as string
% OUTPUTs
% fList = list of function dependencies
%%
% Get list of function dependencies
[fList,~] = matlab.codetools.requiredFilesAndProducts(func);
% Check if destination path exists, if not make it
if ~7 == exist(dest,'dir')
    mkdir(dest)
end
% Copy all dependencies into destination path; skip original function
% 'func'
if size(fList,2) > 1
    for fi = 2:size(fList,2)
        copyfile(fList{fi},dest);
    end
end   