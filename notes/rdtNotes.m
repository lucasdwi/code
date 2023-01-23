files = fileSearch('H:\Shared drives\Green Lab 3\RDM\','.xlsx');
thisData = cell(1,4);
c = 1;
for ii = 1:numel(files)
    [~,sheetName] = xlsfinfo(files{ii});
    for jj = 1:numel(sheetName)
        if ~contains(sheetName{jj},'Summary') && ...
                ~contains(sheetName{jj},'Bias')
            thisSheet = readtable(files{ii},'sheet',sheetName{jj});
            % Get indices of rows that are not: 'no pl2', 'pre-implant',
            % don't have date, 'skip'
            inds = logicFind(1,...
                cellfun(@isempty,strfind(thisSheet.notes,'no pl2')) & ...
                cellfun(@isempty,strfind(thisSheet.notes,'pre-implant'))...
                & ~isnat(thisSheet.Date) & ...
                cellfun(@isempty,strfind(thisSheet.notes,'skip')) & ...
                cellfun(@isempty,strfind(thisSheet.notes,'T3')),'==');
            for k = inds
                if isdatetime(thisSheet{k,1}) && ~isnan(thisSheet{k,2})
                    thisData(c,:) = {sheetName{jj},...
                        char(datetime(thisSheet{k,1},'format',...
                        'yyyy-MM-dd')),thisSheet.AUC(k),...
                        thisSheet.notes(k)};
                    % Grab percent delay choice
                    percentDelay(c,:) = thisSheet{k,[3,5,7,9,11]}./100;
                    trapzAUCs(c) = trapz([0,.25,.5,.75,1],percentDelay(c,:));
                    thisFile(c,:) = thisData(c,1:2);
                    notes(c,1) = thisSheet.notes(k);
                    c = c+1;
                end
            end
        end
    end
end
% Fill in NaNs in empty spots of thisFile so that it can be searched
thisFile(cellfun(@isempty,thisFile(:,1)),1) = deal({'nemo'});
thisFile(cellfun(@isempty,thisFile(:,2)),2) = deal({'nowhen'});
%% Plot all AUCs of shock 0.3
these = cellfun(@(x) contains(x,'0.3'),notes);
figure
histogram(trapzAUCs(these),50)
title('shock = 0.3')
xlabel('AUC')
ylabel('sessions')
box off
%% Plot RDM1 AUCs through time across last intervention
these = [90:102];
figure
hold on
plot(1:5,trapzAUCs(these(1:5)),'ok');
plot(6:13,trapzAUCs(these(6:13)),'or');
xlabel('session')
ylabel('AUC')