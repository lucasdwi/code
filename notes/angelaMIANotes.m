load('C:\Users\Pythia\Documents\GreenLab\data\angelaMIA\analyzed\avgData.mat')
nameVect = names({'ILL','CA1L','PLL','NAcL','PLR','CA1R','ILR','NAcR'},{'d','t','a','b','lg','hg'});
%%
inds = 49:6:216;
s = data(1:15,:);
poly = data(16:36,:);
%%
for ii = 1:size(inds,2)
    [~,p(ii)] = ttest2(s(:,inds(ii)),poly(:,inds(ii)));
    if p(ii)*28 < 0.05
       figure
       plot(ones(1,size(poly,1)),poly(:,inds(ii)),'ob')
       hold on
       plot([0.75 1.25 ],repmat(mean(poly(:,inds(ii))),1,2),'-k')
       
       plot(zeros(1,size(s,1)),s(:,inds(ii)),'or')
       plot([-.25 .25],repmat(mean(s(:,inds(ii))),1,2),'-k')
       
       xlim([-0.5 1.5])
       title(nameVect{inds(ii)})
       set(gca,'XTick',[0 1],'XTickLabel',{'S','Poly'})
       box off
    end
end

