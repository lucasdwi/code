function DBSPowerPlot
%% Spectrogram
cfg = [];
cfg.layout = 'ordered';
cfg.channel = TFRs{1,1}.label;

layout = ft_prepare_layout(cfg,TFRs{1,1});
%%
cfg = [];
cfg.channel = TFRs{1,1}.label;
cfg.layout = layout;

ft_multiplotTFR(cfg,TFRs{1,1})