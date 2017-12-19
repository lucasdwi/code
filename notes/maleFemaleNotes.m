[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'estrus';'di'},{'pow','coh'},'avg','');
estrus = cat(1,data{1,1}{:});
diestrus = cat(1,data{1,2}{:});
data = [estrus;diestrus];
y = [ones(size(estrus,1),1);zeros(size(diestrus,1),1)];
%%
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'di';'male'},{'pow','coh'},'avg','');
female = cat(1,data{1,1}{:});
male = cat(1,data{1,2}{:});
data = [female;male];
y = [ones(size(female,1),1);zeros(size(male,1),1)];
%%
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'estrus';'pro'},{'pow','coh'},'avg','');
estrus = cat(1,data{1,1}{:});
pro = cat(1,data{1,2}{:});
data = [estrus;pro];
y = [ones(size(estrus,1),1);zeros(size(pro,1),1)];
%%
[data,samp,files] = collateData('C:\Users\Pythia\Documents\GreenLab\data\maleFemale\processed\',{'di';'pro'},{'pow','coh'},'avg','');
diestrus = cat(1,data{1,1}{:});
pro = cat(1,data{1,2}{:});
data = [diestrus;pro];
y = [ones(size(diestrus,1),1);zeros(size(pro,1),1)];
%%
load('C:\Users\Lucas\Desktop\maleFemale\phase1.mat')
realAcc = (1-real{1}.err).*100;
rM = mean(realAcc);
rS = std(realAcc);
permAcc = (1-perm{1}.err).*100;
pM = mean(permAcc);
pS = std(permAcc);
[d] = distES(realAcc,permAcc);
%%
figure
histogram(realAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','k','EdgeColor','w')
hold on
histogram(permAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','w')
yticks = get(gca,'Ytick').*100;
set(gca,'YTickLabel',yticks);
ylabel('Percent of Models (%)')
xlabel('Accuracy')
title('Diestrus vs. Estrus')
legend({['Real: ',num2str(round(rM)),'\pm',num2str(round(rS)),'%'],['Permuted: ',num2str(round(pM)),'\pm',num2str(round(pS)),'%']},'Location','northwest')
text(22,.086,['d = ',num2str(round(d,2))])
box off

%%
load('C:\Users\Lucas\Desktop\maleFemale\maleFemale1.mat')
realAcc = (1-real{1}.err).*100;
rM = mean(realAcc);
rS = std(realAcc);
permAcc = (1-perm{1}.err).*100;
pM = mean(permAcc);
pS = std(permAcc);
[d] = distES(realAcc,permAcc);
%%
figure
histogram(realAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','k','EdgeColor','w')
hold on
histogram(permAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','w')
title('Male vs. Female (Estrus)')
ylabel('Percent of Models (%)')
yticks = get(gca,'Ytick').*100;
set(gca,'YTickLabel',yticks);
xlabel('Accuracy (%)')
legend({['Real: ',num2str(round(rM)),'\pm',num2str(round(rS)),'%'],['Permuted: ',num2str(round(pM)),'\pm',num2str(round(pS)),'%']},'Location','northwest')
text(12,.076,['d = ',num2str(round(d,2))])
box off
%%
load('C:\Users\Lucas\Desktop\maleFemale\maleDi1.mat')
realAcc = (1-real{1}.err).*100;
rM = mean(realAcc);
rS = std(realAcc);
permAcc = (1-perm{1}.err).*100;
pM = mean(permAcc);
pS = std(permAcc);
[d] = distES(realAcc,permAcc);
%%
figure
histogram(realAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','k','EdgeColor','w')
hold on
histogram(permAcc,'Normalization','probability','BinWidth',1,'FaceAlpha',1,'FaceColor','w')
title('Male vs. Female (Diestrus)')
ylabel('Percent of Models (%)')
yticks = get(gca,'Ytick').*100;
set(gca,'YTickLabel',yticks);
xlabel('Accuracy (%)')
legend({['Real: ',num2str(round(rM)),'\pm',num2str(round(rS)),'%'],['Permuted: ',num2str(round(pM)),'\pm',num2str(round(pS)),'%']},'Location','northwest')
text(22,.086,['d = ',num2str(round(d,2))])
box off
%%
load('maleDi.mat')
cmbs = nchoosek(1:26,5);
inds = 1:26;
for ii = 1:20
    chk = 1;
    while chk == 1
        samp = cmbs(randi(65780,1,1),:);
        chk = size(unique(y(samp)),1);
    end
    test(:,ii) = samp;
end

for ii = 1:size(data,2)
    for jj = 1:20
        train = inds(~ismember(inds,test(:,jj)));
        mdl = fitglm(data(train,ii),y(train),'distribution','binomial');
        prob = predict(mdl,data(test(:,jj),ii));
        [~,~,~,a(ii,jj)] = perfcurve(y(test(:,jj)),prob,1);
    end
end
%%
mA = mean(a,2);
sA = std(a,[],2);
[smA,inds] = sort(mA,'descend');
ssA = sA(inds);
scatterErr(1:60,smA,ssA,1)
hold on
% Plot 50 line
plot([0 60],[0.5 0.5],'--k','LineWidth',2)

xlabel('Feature')
ylabel('AUC')
title('Average AUC from Univariate Logistic')
%%
load('phase.mat')
cmbs = nchoosek(1:26,5);
inds = 1:26;
for ii = 1:20
    chk = 1;
    while chk == 1
        samp = cmbs(randi(65780,1,1),:);
        chk = size(unique(y(samp)),1);
    end
    test(:,ii) = samp;
end

for ii = 1:size(data,2)
    for jj = 1:20
        train = inds(~ismember(inds,test(:,jj)));
        mdl = fitglm(data(train,ii),y(train),'distribution','binomial');
        prob = predict(mdl,data(test(:,jj),ii));
        [~,~,~,a(ii,jj)] = perfcurve(y(test(:,jj)),prob,1);
    end
end
%%
mA = mean(a,2);
sA = std(a,[],2);
[smA,inds] = sort(mA,'descend');
ssA = sA(inds);
scatterErr(1:60,smA,ssA,1)
hold on
% Plot 50 line
plot([0 60],[0.5 0.5],'--k','LineWidth',2)

xlabel('Feature')
ylabel('AUC')
title('Average Univariate AUC: Estrus vs. Diestrus')