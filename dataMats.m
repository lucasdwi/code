%% Create data matrices for each condition
% Initialize matrices
conds = {'Reg'};
cd('C:\Users\Lucas\Desktop\GreenLab\data\rest2');
%%
for c = 1:length(conds)
    % Find all files of condition c
    files = []; files = dir(strcat('*',conds{c},'*'));
    for f = 1:size(files,1)
        % Load file
        load(files(f).name);
        mast.(conds{c})(f,1) = bingeSize;
        mats.(conds{c})(f,2:51)=horzcat(reshape(relPower,1,20),reshape(relCoh,1,30));
        clear relPower psdTrls fds avgCoh relCoh
    end
end

        

