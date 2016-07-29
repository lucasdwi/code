function [fd] = ConvertPl2All_Files(dir)
%   Runs Pl2tomvdm.m on all files within the directory. Written by
%   JJS, edited by LLD

% dir = directory within which to search for files
cd(dir);
fd = FindFiles('*pl2');

for iSess = 1:length(fd);
    pushdir(fileparts(fd{iSess}));
    [~, name, ~] = fileparts(fd{iSess});
    disp(iSess)
    disp(name)
    if exist(strcat(fd{iSess},'.mat')) ==2;
        disp('matfile already exists...skipping')
    else
        %[~, ~] = Pl2tomvdm(fd{iSess});
        [~, ~] = Pl2tomvdmGenFile(fd{iSess});
    end
    popdir;
end

