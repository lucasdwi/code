function [tpr,fpr] = rocCurve(target,output,thresh)
%%
% targetT = logicFind(1,target,'==');
% targetF = logicFind(0,target,'==');
for ti = 1:length(thresh)
    thisRound(ti,:) = threshRound(output,thresh(ti),[0 1]);
    ind1 = logicFind(1,thisRound(ti,:),'==');
    ind0 = logicFind(0,thisRound(ti,:),'==');
    fp(ti) = sum(target(ind1) == 0);
    tp(ti) = sum(target(ind1) == 1);
    fn(ti) = sum(target(ind0) == 1);
    tn(ti) = sum(target(ind0) == 0);
    fpr(ti) = fp(ti)/(fp(ti)+tn(ti));
    tpr(ti) = tp(ti)/(tp(ti)+fn(ti));
end
figure
plot(fpr,tpr)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
