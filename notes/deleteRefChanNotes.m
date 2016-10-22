stimFiles = dir('*.mat*');
for f = 1:size(stimFiles,1)
    load(stimFiles(f).name)
    % If four reference channels exist, get rid of
    if size(LFPTs.data,1) == 8
       LFPTs.data = LFPTs.data(2:2:8,:);
    end
    % If four reference channel labels exist, get rid of
    if size(LFPTs.label,2) == 8
        LFPTs.label = LFPTs.label(:,2:2:8);
    end
    save(strcat('C:/Users/Lucas/Desktop/GreenLab/data/twoSiteStim/edited/',stimFiles(f).name),'ad','adfreq','eventTs','LFPTs','pl2','TimeSampEr','ts','WBchan');
end