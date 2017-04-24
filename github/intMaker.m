function [overInt] = intMaker(times)
%% Makes intervals from onsets/offsets and time stamps
% INPUTS:
% times = cell row array by channel; format of each cell = row 1 is
% start times, row 2 is stop times 
%%
chans = size(times,2);
% Merge intervals
overInt = cell(1,chans); % Creates interval structure
for c = 1:chans
    if ~isempty(times{1,c})
        s = times{1,c}(1,1);
        e = times{1,c}(2,1);
        k = 0;
        for j=1:length(times{1,c})
            if times{1,c}(1,j)> e
                k = k+1;
                overInt{1,c}(1,k) = s;
                overInt{1,c}(2,k) = e;
                s = times{1,c}(1,j);
                e = times{1,c}(2,j);
            else
                e = max(e,times{1,c}(2,j));
            end
        end
        % Stack last interval
        k = k+1;
        overInt{1,c}(1,k) = s;
        overInt{1,c}(2,k) = e;
    end
end
