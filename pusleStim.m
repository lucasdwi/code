chan = [1:16];
f = repmat(12,1,16);
params = cat(1,repmat([100:50:200;-100:-50:-200;90,90,90;90,90,90;25,25,...
    25],1,1,3),permute(repmat(2:4,1,1,3),[1,3,2]));
nRep = 60;

PS_InitAllStim
stimNum = PS_GetNStim;

NChan = PS_GetNChannels(stimNum);
PS_SetMonitorChannel(1,chan(1))

combs = [];
c=1;
for ii = 1:3
    for jj = 1:3
        combs(c,:) = [ii,jj];
        c = c+1;
    end
end
% combs = repmat(combs,1,1,nRep);
for ii = 1:nRep
    allCombs(:,:,ii) = combs(randperm(9,9),:);
end
%% dual site
% Randomize parameters
for k = 1:nRep
    allCombs(:,:,k) = combs(randperm(9,9),:);
end
pause(600)
for r = 1:nRep
    disp(r)
    for ii = 1:9
        for jj = 1:size(chan,2)
            PS_SetRate(1,chan(jj),f(jj));
            PS_SetRepetitions(1,chan(jj),params(6,allCombs(ii,1,r),...
                allCombs(ii,2,r)));
            PS_SetPatternType(1,chan(jj),0);
            PS_SetRectParam(1,chan(jj),params(1:5,allCombs(ii,1,r),...
                allCombs(ii,2,r))');
            PS_LoadChannel(1,chan(jj));
        end
        % First stim site channels
        PS_StartStimChannel(1,chan(1));
        PS_StartStimChannel(1,chan(2));
        PS_StartStimChannel(1,chan(5));
        PS_StartStimChannel(1,chan(6));
        PS_StartStimChannel(1,chan(9));
        PS_StartStimChannel(1,chan(10));
        PS_StartStimChannel(1,chan(13));
        PS_StartStimChannel(1,chan(14));
%         pause(0.0417)
        % Pause for 0.0022 seconds shorter than target (0.0417) to account
        % for computer speed
        pause(0.0395)
        % Second stim site channels
        PS_StartStimChannel(1,chan(3));
        PS_StartStimChannel(1,chan(4));
        PS_StartStimChannel(1,chan(7));
        PS_StartStimChannel(1,chan(8));
        PS_StartStimChannel(1,chan(11));
        PS_StartStimChannel(1,chan(12));
        PS_StartStimChannel(1,chan(15));
        PS_StartStimChannel(1,chan(16));
        pause(6)
    end
end
pause(600)
%%
for site = 1:2
    % Randomize parameters
    for k = 1:nRep
        allCombs(:,:,k,site) = combs(randperm(9,9),:);
    end
    for r = 1:nRep
        disp(r)
        for ii = 1:9
            for jj = 1:size(chan,2)
                PS_SetRate(1,chan(jj),f(jj));
                PS_SetRepetitions(1,chan(jj),params(6,allCombs(ii,1,r,site),allCombs(ii,2,r,site)));
                PS_SetPatternType(1,chan(jj),0);
                PS_SetRectParam(1,chan(jj),params(1:5,allCombs(ii,1,r,site),allCombs(ii,2,r,site))');
                PS_LoadChannel(1,chan(jj));
            end
            % First stim site channels
            if site == 1
                PS_StartStimChannel(1,chan(1));
                PS_StartStimChannel(1,chan(2));
                PS_StartStimChannel(1,chan(5));
                PS_StartStimChannel(1,chan(6));
                PS_StartStimChannel(1,chan(9));
                PS_StartStimChannel(1,chan(10));
                PS_StartStimChannel(1,chan(13));
                PS_StartStimChannel(1,chan(14));
            % Second stim site channels
            elseif site == 2
                PS_StartStimChannel(1,chan(3));
                PS_StartStimChannel(1,chan(4));
                PS_StartStimChannel(1,chan(7));
                PS_StartStimChannel(1,chan(8));
                PS_StartStimChannel(1,chan(11));
                PS_StartStimChannel(1,chan(12));
                PS_StartStimChannel(1,chan(15));
                PS_StartStimChannel(1,chan(16));
            end
            pause(6)
        end
    end
end
