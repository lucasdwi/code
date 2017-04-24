
% Creates data structure with information on files with wildcard fType
% in name; if >1 fType, then concatentates together
files = [];
cd('C:\Users\Lucas\Desktop\GreenLab\data\paper2\to_analyze\')
thisF = dir(strcat('*mat*'));
% Just pull out 'name' field and concatenate to 'files'
files = vertcat(files,extractfield(thisF,'name')');

%% Run spectcompbase.m
for i = 1:length(files)
    tic
    disp(['Running spectcompbase on file ',num2str(i),' out of ',num2str(length(files)),': ',files{i}])
    [sdir,file,filter,dsf,thresh,onset,offset,minInt,foi,bands,cycles,ftimwin,overlap,cohMethod,eoi,saveParent] = scbParamsMulti(files{i});
    cd(sdir)
    load(file)
    chans = size(LFPTs.data,1);
    [LFPTs,chk_nan,zeroedChannel,clnTrls,clnEvents,trls,adfreq] = preProcess(LFPTs,adfreq,dsf,thresh,onset,offset,minInt,eoi,eventTs);
    overlap = 1 - overlap;
    cmb = nchoosek(1:chans,2);
    for c = 1:size(cmb,1)
        channelCmb(c,:) = LFPTs.label(cmb(c,:));
    end
    for e = 1:size(eoi,1)
        for b = 1:size(bands,1)
            tic
            disp(['Computing CSD with ft_freqanalysis.mat for ',bands{b,1},' band...'])
            cfg              = [];
            cfg.output       = 'powandcsd';
            cfg.method       = 'mtmconvol';
            cfg.taper        = 'hanning';
            cfg.foi          = foi(1):foi(2):foi(3); % frequencies to use
            % Use frequency dependent windows (n cycles per window, computed at
            % start, 'cycFtimwin') with (x%) overlap
            if ~isempty(cycles)
                cfg.t_ftimwin    = ones(size(cfg.foi)).*(cycles/bands{b,2}(1));
                minTWin          = min(cfg.t_ftimwin)*overlap;
                cfg.toi          = eoi{e,2}(1):minTWin:eoi{e,2}(2);
                % Or use a constant size for windows with (x%) overlap to compute
                % cycles and apply forward
            else
                cfg.t_ftimwin    = ones(size(cfg.foi)).*(cycFtimwin/bands{b,2}(1));
                cfg.toi          = eoi{e,2}(1):ftimwin*overlap:eoi{e,2}(2);
            end
            cfg.keeptrials   = 'yes';
            cfg.channel      = LFPTs.label;
            cfg.channelcmb   = channelCmb;
            
            powCorrTFR{e,b} = ft_freqanalysis(cfg,trls{1,e});
            toc
        end
        % Run powerCorr
        disp('Running powerCorr.m...')
        tic
        cfg.trialwindows = 'yes';
        [STDCorr{e},MeanCorr{e},TWCorr{e},powerCorrSort{e},~,~,~,~] = powerCorr(powCorrTFR(e,:),bands,cfg);
        toc
    end
    [~,name,~] = fileparts(strcat(sdir,file));
    save(['C:\Users\Lucas\Desktop\GreenLab\data\paper2\to_analyze\',name,'_powCorr.mat'],'powerCorrSort')
    clearvars -except i files
end
%% Merge files
[~,fNames] = fileSearch({'C:\Users\Lucas\Desktop\GreenLab\data\paper2\withPowCorr\new\'},{'powCorr','binge'});
for fi = 1:size(fNames{1,1},2)
    disp(num2str(fi))
    load(fNames{1,2}{1,fi})
    load(fNames{1,1}{1,fi})
    save([fNames{1,1}{1,fi},'extra.mat'],'powerCorrSort','LFPTs','trls','relPower','psdTrls','coh','stdPower','hist','powerEventComp')
end