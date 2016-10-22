coreGroup = {'coreRespond','coreStrict'};
nonCI = logicFind(0,T2.Base.coreRespond,'==');
respCI = logicFind(1,T2.Base.coreRespond,'==');
strictCI = logicFind(1,T2.Base.coreStrict,'==');
shellGroup = {'shellRespond','shellStrict'};
nonSI = logicFind(0,T2.Base.shellRespond,'==');
respSI = logicFind(1,T2.Base.shellRespond,'==');
strictSI = logicFind(1,T2.Base.shellStrict,'==');
%%
% T-tests
for ii = 6:55;
    [h{ii,1},p{ii,1}] = ttest2(table2array(T2.Base(nonCI,ii)),table2array(T2.Base(respCI,ii)));
    [h{ii,2},p{ii,2}] = ttest2(table2array(T2.Base(nonCI,ii)),table2array(T2.Base(strictCI,ii)));
    [h{ii,3},p{ii,3}] = ttest2(table2array(T2.Base(respCI,ii)),table2array(T2.Base(strictCI,ii)));
end
ch = h(6:55,:); cp = p(6:55,:);
for ii = 6:55;
    [h{ii,1},p{ii,1}] = ttest2(table2array(T2.Base(nonSI,ii)),table2array(T2.Base(respSI,ii)));
    [h{ii,2},p{ii,2}] = ttest2(table2array(T2.Base(nonSI,ii)),table2array(T2.Base(strictSI,ii)));
    [h{ii,3},p{ii,3}] = ttest2(table2array(T2.Base(respSI,ii)),table2array(T2.Base(strictSI,ii)));
end
sh = h(6:55,:); sp = p(6:55,:);

ps = {cp,sp};
%% Plots
for thisp = 1:2    
    for ii = 1:size(ps{thisp},1)
        if ((ps{thisp}{ii,1}<=0.05) + (ps{thisp}{ii,2}<=0.05) + (ps{thisp}{ii,3}<=0.05)) >= 1
            if thisp == 1
                data(:,1) = [table2array(T2.Base(nonCI,ii+5));table2array(T2.Base(respCI,ii+5));table2array(T2.Base(strictCI,ii+5))];
                len1 = length(nonCI); len2 = length(respCI); len3 = length(strictCI);
            else if thisp == 2
                    data(:,1) = [table2array(T2.Base(nonSI,ii+5));table2array(T2.Base(respSI,ii+5));table2array(T2.Base(strictSI,ii+5))];
                    len1 = length(nonSI); len2 = length(respSI); len3 = length(strictSI);
                end
            end
            data(1:len1,2) = 1; data(len1+1:len1+len2,2) = 2; data(len1+len2+1:len1+len2+len3,2) = 3;
            figure;
            b{ii} = boxplot(data(:,1),data(:,2),'Labels',{'Non-Responders','Responders','Strict'});
            yt = get(gca,'Ytick'); xt = get(gca,'Xtick');
            yTick = abs(yt(2)-yt(1));
            yMax = [get(b{ii}(3,1),'YData');get(b{ii}(3,2),'YData');get(b{ii}(3,3),'YData')];
            hold on;
            for j = 1:3
                bar = {[1 2],[1 3],[2 3]};
                if ps{thisp}{ii,j} <= 0.05
                    if j == 1 || j == 3
                        plot(xt(bar{j}),[1 1]*max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+0.5*yTick,'-k',mean(xt(bar{j})),max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+0.8*yTick,'*k');
                        text(mean(xt(bar{j}))+0.1,max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+0.8*yTick,strcat('p \approx ',num2str(round(ps{thisp}{ii,j},4))));
                    end
                    if j == 2
                        plot(xt(bar{j}),[1 1]*max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+yTick,'-k',mean(xt(bar{j})),max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+1.3*yTick,'*k');
                        text(mean(xt(bar{j}))+0.1,max(yMax(bar{j}(1),1),yMax(bar{j}(2),1))+1.3*yTick,strcat('p \approx ',num2str(round(ps{thisp}{ii,j},4))));
                    end
                end
                plot(j,mean(data(data(:,2)==j,1)),'gsq');
            end
            lim = ylim;
            if lim(2) < max(yMax(:,1))+yTick*1.5
                ylim([lim(1) lim(2)+yTick*1.5]);
            end
            if thisp == 1
                title([strrep(T2.Base.Properties.VariableNames{ii+4},'_',' '),': Core']);
            else if thisp == 2
                    title([strrep(T2.Base.Properties.VariableNames{ii+4},'_',' '),': Shell']);
                end
            end
        end
    end
end