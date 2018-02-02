function [fNames] = ConvertPl2All_Files(dir)
%% Runs Pl2tomvdmGenFile.m on all files within the directory.
%__________________________________________________________________________
% INPUTS:
% dir = source directory of files to be converted; format = string
%__________________________________________________________________________
% OUTPUTS:
% fNames = file names that have been converted
%__________________________________________________________________________
% Written by JJS, edited by LLD (got rid of 'sd' and replaced all MClust
% functions) 2017
%%
fNames = fileSearch(dir,'pl2');
for iFile = 1:size(fNames,2)
    [~, name, ~] = fileparts(fNames{iFile});
    disp(iFile)
    disp(name)
    if exist(strcat(fNames{iFile},'.mat'),'file') == 2
        disp('.mat file already exists...skipping')
    else
        [~,sd] = Pl2tomvdmGenFile(fNames{iFile});
        % Extract variables from sd to save
        ad = sd.ad; adfreq = sd.adfreq; eventTs = sd.eventTs; fn = sd.fn; lfpchan = sd.lfpchan; LFPTs = sd.LFPTs; n = sd.n; pl2 = sd.pl2; TimeSampEr = sd.TimeSampEr; ts = sd.ts; WBchan = sd.WBchan;  %#ok<NASGU>
        save(strcat(name, '.mat'), 'ad','adfreq','eventTs','fn','lfpchan','LFPTs','n','pl2','TimeSampEr','ts','WBchan','-v7.3');
    end
end
