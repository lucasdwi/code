%% Bootstrapping
for a = 2:5
    for c = 2:size(T3{1,a},2)
        [bootStat{a-1}(:,c-1),bootdist{a-1}(:,:,c-1)] = bootstrp(10000,@mean,T3{1,a}(:,c)-mean(T3{1,a}(:,c)));
        %[thisStat,bootdist{a-1}(:,:,c-1)] = bootstrp(10000,@(x)[mean(x) std(x)],T3{1,a}(:,c));
        %bootMean{a-1}(:,c-1) = thisStat(:,1);
        %bootStd{a-1}(:,c-1) = thisStat(:,2);
    end
    disp(a-1)
end
%% T-Statistic
for a = 1:size(bootStat,2)
   for b = 1:size(bootStat{1,a},2)
       t{1,a}(b) = (mean(bootStat{1,a}(:,b))-mean(T3{1,a+1}(:,b+1)))/std(bootStat{1,a}(:,b));
       p{1,a}(b) = 2*tcdf(-abs(t{1,a}(b)),size(bootStat{1,a},1)-1);
   end
   pAdj(a,:) = mafdr(p{1,a}(1,:)','BHFDR','true')';
end
%%
figure
for a = 1:4
    for ii = 1:144
        if basePs(a,ii) <= 0.05
            hold on; plot(ii,a*0.5,'ok')
            if sum(basePs(:,ii)<=0.05) == 1
                hold on; plot(ii,a*0.5,'or')
            end
            if sum(basePs(:,ii)<=0.05) == 2
                hold on; plot(ii,a*0.5,'og')
            end
            if sum(basePs(:,ii)<=0.05) == 3
                hold on; plot(ii,a*0.5,'ob')
            end
        end
    end
end
ylim([0 2.5])
set(gca,'YTickLabel',{'','Base a','Base b','Base c','Base Boot',''})

