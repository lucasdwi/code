% Used for files with 4 animals and 3 stim sessions; splits apart the
% animals and the stim blocks from baseline through washout with stim
% intervals included
files = fileSearch('G:\GreenLab\data\lsdStim\mat\','NAcs');
for k = 4:numel(files)
    load(files{k})
    parts = strsplit(files{k},'_');
    these = sort([eventTs.t{9};eventTs.t{11}],'ascend');
    % Stop stim block
    blockInds = logicFind(600,round(diff(these)),'>=');
    blockStops = [these(blockInds);these(end)]+0.005;
    % Start stim block  
    blockStarts = [these(1);these(blockInds+1)];
    % Stim chunks
    stimStarts = these(logicFind(10,round(diff(these)),'=='));
    % Baselines
    baseStart = [LFPTs.tvec(1),LFPTs.tvec(nearest_idx3(blockStops(1:2),...
        LFPTs.tvec)+1)];
    baseStop = LFPTs.tvec(nearest_idx3(blockStarts,LFPTs.tvec)-1);
    % Washouts
    washStart = LFPTs.tvec(nearest_idx3(blockStops,LFPTs.tvec)+1);
    washStop = [LFPTs.tvec(nearest_idx3(blockStarts(2:3),LFPTs.tvec)-1),...
        LFPTs.tvec(end)];
    % EventTs
    oldE = eventTs;
    oldL = LFPTs;
    for ii = 1:3
        eventTs = [];
        eventTs.label{1} = 'Base (Start)'; eventTs.label{2} = 'Base (End)';
        eventTs.t{1} = baseStart(ii); eventTs.t{2} = baseStop(ii);
        eventTs.label{3} = 'Stim (Start)'; eventTs.label{4} = 'Stim (End)';
        eventTs.t{3} = stimStarts(29*ii-28:29*ii);
        eventTs.t{4} = stimStarts(29*ii-28:29*ii)+10;
        eventTs.label{5} = 'Wash (Start)'; eventTs.label{6} = 'Wash (End)';
        eventTs.t{5} = washStart(ii); eventTs.t{6} = washStop(ii);
        for jj = 1:4
            LFPTs = [];
            LFPTs.tvec = oldL.tvec(nearest_idx3(baseStart(ii),oldL.tvec)...
                :nearest_idx3(washStop(ii),oldL.tvec));
            LFPTs.data = oldL.data(8*jj-7:8*jj,nearest_idx3(...
                baseStart(ii),oldL.tvec):nearest_idx3(washStop(ii),...
                oldL.tvec));
            LFPTs.label = oldL.label(8*jj-7:8*jj);
            LFPTs.cfg = oldL.cfg;
            save(['G:\GreenLab\data\lsdStim\mat\split\',strjoin(parts(...
                [jj,5,5+ii,9]),'_')],'LFPTs','eventTs','adfreq')
        end
    end
end