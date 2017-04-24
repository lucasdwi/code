%% Old - several models   
% models = cell(1,length(allPredict)*length(allResponse));
% 
% 
% if modelNum == 1
%     for i = 1:length(allResponse)
%         for j = 1:length(allPredict)
%             thismd = fitlm(T,strcat(allResponse{i},'~',allPredict{j}));
%             thisind = j+(i-1)*length(allPredict);
%             models{thisind} = thismd;
%         end
%     end
% else if modelNum == 2
%         for i = 1:length(allResponse)
%            for j = 1:length(allPredict)
%                thismd = fitlm(T,strcat(allResponse{i},'~',allPredict{j},'+group'));
%                thisind = j+((i-1)*length(allPredict));
%                models{thisind} = thismd;
%            end
%         end
%     else if modelNum == 3
%             grouping = {'neither','core','shell'};
%             for i = 1:length(allResponse)
%                 for j = 1:length(grouping)
%                     for k = 1:length(allPredict)
%                         thismd = fitlm(T(T.group == grouping{j},:),strcat(allResponse{i},'~',allPredict{k}));
%                         thisind = (k+(j-1)*length(allPredict))+(length(grouping)*length(allPredict)*(i-1));
%                         models{1,thisind} = thismd;
%                         models{2,thisind} = grouping{j};
%                     end
%                 end
%             end
%         end
%     end
% end
% % Create stats table
% % Setup empty table
% Tstats = cell2table({});
% %% Fill table with stats from models
% for i = 1:length(models)
%     if modelNum == 3
%         thisName = strcat(models{2,i},'_',models{1,i}.ResponseName,'~',cellstr(models{1,i}.Formula.LinearPredictor));
%     else 
%         thisName = strcat(models{1,i}.ResponseName,'~',cellstr(models{1,i}.Formula.LinearPredictor));
%     end
%     thisTstat = table(models{1,i}.Rsquared.Ordinary',models{1,i}.Rsquared.Adjusted',coefTest(models{1,i})',i,'RowNames',thisName);
%     Tstats(thisName,:) = thisTstat; %vertcat(Tstats,thisTstat);
% end
% % Set RowNames to model formula and set VariableNames
% Tstats.Properties.VariableNames = {'R_Ord','R_Adj','pValue','ModelNumber'};
% % Sort table by p-values
% Tstats = sortrows(Tstats,'pValue','ascend');
% %% Find model numbers with p-values below alpha corresponding to either core or shell response
% either = Tstats(cellfun(@isempty,strfind(cellstr(Tstats.Properties.RowNames),'neither')),:);