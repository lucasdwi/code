%% NaN contiguous data intervals less than minInt
tic
dataInt = 1;
clnTrls = cell(numMark,largest); %Preallocate clnTrls
for i = 1:length(markers)
    for j = 1:numel(eventTs.t{markers(i)})
        for k = 1:chans
            A = intTime{i,j}(k,:); A(~isnan(A)) = 1; A(isnan(A)) = 0;
            dataStart = find(diff(A)==1)+1;
            dataStop = find(diff(A)==-1);
            if ~isnan(intTime{i,j}(k,1))
                dataStart = [1, dataStart];
            end
            if ~isnan(intTime{i,j}(k,end))
                dataStop = horzcat(dataStop,length(intTime{i,j}));
            end
            
            if ~isnan(sum(intTime{i,j}(k,:))) && length(intTime{i,j}(k,:)) >= minInt %No NaNed data and longer than minInt
                numTrls = floor(length(intTime{i,j}(k,:))/minInt);
                thisTrls = intTime{i,j}(k,1:(minInt*numTrls));
                clnTrls{i,j}(k,:) = thisTrls;
            else for intInd = 1:length(dataStart) %Run through data intervals
                    intLen = dataStop(intInd) - dataStart(intInd);
                    thisTrls = [];
                    if intLen >= minInt %Keep if big enough
                        numTrls = floor(intLen/minInt);
                        thisTrls = horzcat(thisTrls,intTime{i,j}(k,dataStart(intInd):(dataStart(intInd)+minInt*numTrls)-1));
                    end
                end
                if ~isempty(thisTrls)
                    clnTrls{i,j}(k,:) = thisTrls;
                end
            end
        end
    end
end
toc                
                