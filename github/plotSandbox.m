%%
for ii = 1:4
    minusPSD{ii} =  (wAvgPSD{ii}-restPSD{ii})./restPSD{ii};
    minusCoh{ii} = (wAvgCoh{ii}-restCoh{ii})./restCoh{ii};
    minusRel{ii} = (wAvgRelPSD{ii}-restRelPSD{ii})./restRelPSD{ii};
end
%%
figure
for ii = 1:4
    subplot(2,2,ii)
    hold on
    for j = 1:4
        plot(minusPSD{j}(ii,:))
        title(['Channel ',num2str(ii)])
    end
end
legend('Base','Dep24','Dep48','Chow')
%%
figure
for ii = 1:4
    subplot(2,2,ii)
    hold on
    bar([minusRel{1}(:,ii),minusRel{2}(:,ii),minusRel{3}(:,ii),minusRel{4}(:,ii)]);
    title(['Channel ',num2str(ii)])
end
legend('Base','Dep24','Dep48','Chow')
%%
cmb = nchoosek(1:4,2);
figure
for ii = 1:6
    subplot(2,3,ii)
    hold on
    for j = 1:4
        plot(minusCoh{j}(ii,:))
    end
    title(['Channel ',num2str(cmb(ii,1)),'-Channel ',num2str(cmb(ii,2))]);
end
legend('Base','Dep24','Dep48','Chow')