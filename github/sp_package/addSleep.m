%% Add 'sleep' to data file
% Get file lists
cd('C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\edited\')
dataF = dir('*.mat*');
cd('C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\markedSleep\')
sleepF = dir('*.mat*');
% Go through load each and save together
if size(dataF) == size(sleepF)
    for fi = 1:size(dataF,1)
        % Load data
        cd('C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\edited\')
        load(dataF(fi,1).name);
        % Load sleep using data file name
        cd('C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\markedSleep\')
        load([dataF(fi,1).name(1:end-4),'_InstAmp_sleep.mat']);
        % Save together
        disp(['Saving... ',dataF(fi,1).name])
        save(['C:\Users\Lucas\Desktop\GreenLab\data\twoSiteStim\editedWithSleep\',dataF(fi,1).name,'.mat'],'ad','adfreq','eventTs','hist','LFPTs','pl2','sleep','TimeSampEr');
        clearvars -except fi dataF sleepF 
    end
else
    error('Inequal number of files in directories - data will not match.')
end

