function [spectroPlots] = spectroComp(trl1,trl2,TFR_event1,TFR_event2,eventLabel)
%% Event 1
% Define plotting layout
cfg = [];
cfg.layout = 'ordered'; cfg.channel = trl1.label; cfg.showlabels = 'yes';
layout = ft_prepare_layout(cfg,trl1);
% Plot
spectroPlots{1} = figure;
cfg = []; cfg.channel = trl1.label; cfg.layout = layout;
ft_multiplotTFR(cfg,TFR_event1);
title(eventLabel{1});
%% Event 2
% Define plotting layout
cfg = [];
cfg.layout = 'ordered'; cfg.channel = trl2.label; cfg.showlabels = 'yes';
layout = ft_prepare_layout(cfg,trl2);
% Plot
spectroPlots{2} = figure;
cfg = []; cfg.channel = trl2.label; cfg.layout = layout;
ft_multiplotTFR(cfg,TFR_event2);
title(eventLabel{2});
