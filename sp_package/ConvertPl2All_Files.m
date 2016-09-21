function [fd] = ConvertPl2All_Files(dir)
%   Runs Pl2tomvdm.m on all files within the directory. Written by
%   JJS, edited by LLD

% Check if MClust and other toolboxes (van der Meer) from path
pathCell = regexp(path,pathsep,'split');
if sum(cell2mat(regexp(pathCell,'MClust'))) ~= 0
    rmpath(genpath('C:\Users\Lucas\Documents\GitHub\neuraldata-w16\toolboxes\MClust-4.3\'));
    rmpath(genpath('C:\Users\Lucas\Documents\GitHub\neuraldata-w16\toolboxes\MClust-3.5\'));
end
%%
% dir = directory within which to search for files
cd(dir);
fd = FindFiles('*pl2');

for iSess = 1:length(fd);
    pushdir(fileparts(fd{iSess}));
    [~, name, ~] = fileparts(fd{iSess});
    disp(iSess)
    disp(name)
    if exist(strcat(fd{iSess},'.mat')) == 2;
        disp('matfile already exists...skipping')
    else
        %[~, ~] = Pl2tomvdm(fd{iSess});
        [~, sd] = Pl2tomvdmGenFile(fd{iSess});
        ad = sd.ad; adfreq = sd.adfreq; eventTs = sd.eventTs; fn = sd.fn; lfpchan = sd.lfpchan; LFPTs = sd.LFPTs; n = sd.n; pl2 = sd.pl2; TimeSampEr = sd.TimeSampEr; ts = sd.ts; WBchan = sd.WBchan;
        disp('Saving file...')
        save(strcat(name, '.mat'),'ad','adfreq','eventTs','fn','lfpchan','LFPTs','n','pl2','TimeSampEr','ts','WBchan');
    end
    popdir;
end

