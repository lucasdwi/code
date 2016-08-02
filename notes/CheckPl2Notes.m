%% Open and plot data
cd('C:\Users\Lucas\Desktop\GreenLab\data\WilderBinge\second_batch')
name = 'I4BaseDec15';
load(strcat(name,'.mat'));
figure;
for ii = 1:size(LFPTs.data,1)
    subplot(2,2,ii)
    plot(LFPTs.data(ii,:))
    ylim([-1 1]);
    if ii == 1
        title(name);
    end
end