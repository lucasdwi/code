% Input sex info
load('F:/irdmRound2/irdmSex.mat')
files = fileSearch('F:/irdmRound2/processedContinuous3','.mat');
for ii = 1:numel(files)
    load(files{ii},'hist')
    parts = strsplit(files{ii},'_');
    name = parts{1};
    if ~isempty(strfind(files{ii},'-'))
        nameParts = strsplit(name,'-');
        name = nameParts{1};
        condition = join(nameParts{2:end},'-');
    end
%     task = parts{2};
%     recordingMode = parts{3};
%     date = parts{4};
%     
%     hist.name = name;
%     hist.condition = condition;
%     hist.recordingMode = recordingMode;
    hist.sex = sex{logicFind(name,sex,'=='),2};
    save(files{ii},'hist','-append')
end