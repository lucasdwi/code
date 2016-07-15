%% Plot top x sigmodels
top = 10;
f = figure;
for i = 1:top
    subplot(2,5,i)
    h{i} = plot(models.Base{1,sigModels.Base.ModelNumber(i)});
    title(sigModels.Base.Group{i});
    l = legend('show'); set(l,'visible','off');
    text(0.5,0.1,strcat('p \approx ',num2str(round(sigModels.Base.pValue(i),4))),'Units','normalized','HorizontalAlignment','center');
end
tightfig(f)
%% Plot sigmodels 1 and 6
figure;
i = 1;
h{i} = plot(models.Base{1,sigModels.Base.ModelNumber(i)});
title(sigModels.Base.Group{i});
l = legend('show'); set(l,'visible','off');
text(0.5,0.1,strcat('p \approx ',num2str(round(sigModels.Base.pValue(i),4))),'Units','normalized','HorizontalAlignment','center');
figure;
i = 6;
h{i} = plot(models.Base{1,sigModels.Base.ModelNumber(i)});
title(sigModels.Base.Group{i});
l = legend('show'); set(l,'visible','off');
text(0.5,0.1,strcat('p \approx ',num2str(round(sigModels.Base.pValue(i),4))),'Units','normalized','HorizontalAlignment','center');