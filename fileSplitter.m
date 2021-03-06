function fileSplitter(sdir,searchStr,chan,nAnimal,nEvent)
files = fileSearch(sdir,searchStr);
for fI = 1:size(files,2)
    load(files{fI})
    % Grab this filename
    thisFile = files{fI};
    % Remove file extension by finding period; '.'
    thisFile = thisFile(1:strfind(thisFile,'.')-1);
    % Split str at all underscores; '_'
    parts = strsplit(thisFile,'_');
    % Grab date
    d = parts{end};
    % Check that the number of sub-files and expected channels per sub-file
    % matches the total amount of available data
    if ~isequal(size(inds,2)*chan,size(LFPTs.data,1))
        error(['There are a different number of channels than expected'... 
            ' by the number of animals.'])
    end
    % Store LFPTs
    oldLFPTs = LFPTs;
    % Go through each sub-file, pull out LFPTs data and save
    for ii = 1:size(nAnimal,2)
        clear LFPTs;
        LFPTs.type = oldLFPTs.type;
        LFPTs.tvec = oldLFPTs.tvec;
        LFPTs.data = oldLFPTs.data(chan*ii-chan+1:chan*ii,:);
        LFPTs.label = oldLFPTs.label(chan*ii-chan+1:chan*ii);
        LFPTs.cfg = oldLFPTs.cfg;
%         save(strjoin([parts{1},nums{inds(ii)},parts{size(inds,2)+ii+1},...
%             parts{end}],'_'),'LFPTs','adfreq','eventTs','pl2')
        save(strjoin([parts{1},nums{inds(ii)},parts{size(inds,2)+ii+1},...
            parts{end}],'_'),'LFPTs','adfreq','eventTs','pl2')
    end
end