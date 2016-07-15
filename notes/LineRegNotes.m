%% All combinations
for f = 1:numel(cond)
    mI =0;
    for i = 1:length(allGroups)
        for j = 1:length(allResponse)
            for k = 1:length(allPredict)
                thismd = fitlm(T.(cond{f})(T.(cond{f}).(allGroups{i}) == 1,:),strcat(allResponse{j},'~',allPredict{k}));
                mI = mI +1;
                thisind = mI;
                models.(cond{f}){1,thisind} = thismd;
                models.(cond{f}){2,thisind} = allGroups{i};
            end
        end
    end
end