load('C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\analyzed\miaPerm')
masterPerm = [];
for ii = 1:10
    masterPerm = [masterPerm,perm.allLambda{ii}.allErr];
end
%%
figure
hold on
histogram(100-real.allLambda{1}.allErr.*100,'normalization','probability','binwidth',1)
histogram(100-reshape(masterPerm,1,1000).*100,'normalization','probability','binwidth',1)
legend({'Real','Permuted'},'location','nw')
xlabel('Accuracy (%)')
ylabel('Proportion of Models')
%%
s = data([1:15,35:36],:);
poly = data(16:34,:);
[~,p,bulkP] = bulkT(s,poly,0,'fdr');
[tp,ind] = sort(bulkP,'ascend');

label = {'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'};
nameVect = names(label,{'d','t','a','b','lg','hg'});

feats = nameVect(ind);
% Direction 1 = poly>s; 0 = poly<s
for ii = 1:7
   dir(ii) = mean(s(:,ind(ii)))<mean(poly(:,ind(ii)));
end
%%
for ii = 1:7
   figure
   boxplot([poly(:,ind(ii));s(:,ind(ii))],[ones(19,1);zeros(17,1)])
   set(gca,'XTickLabel',{'Saline','MIA'})
   title(nameVect(ind(ii)))
   box off
end