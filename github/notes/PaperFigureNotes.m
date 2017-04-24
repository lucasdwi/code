%% Raw signal figure
figure; plot(LFPTs.tvec,LFPTs.data(1,:),'-k'); xlim([1 2]); ylim([-.31 .31]);
xlabel('Time (s)'); ylabel('mV'); title('Raw Signal'); 
%%
figure;
F = psdTrls.F;
notchInd = [nearest_idx(59,F);nearest_idx(61.5,F)];
c=axes('Position',[.1 .1 .8 1e-12],'Color','none');
set(c,'Units','normalized');

b=axes('Position',[.1 .1 .8 1e-12]);
set(b,'Units','normalized');
set(b,'Color','none');

a=axes('Position',[.1 .2 .8 .7]);
set(a,'Units','normalized');

events = {'event1','event2'};
clrs = {'r','k'};
for i = 1:2
    for j = 1
        hold on; ax{i} = plot(F,psdTrls.(events{i}).Overall(j,:),'color',clrs{i});
        % Fill NaN data in plot
        thisstart = [F(notchInd(1)-1),F(notchInd(2)+1)];
        thisend = [psdTrls.(events{i}).Overall(j,notchInd(1)-1),psdTrls.(events{i}).Overall(j,notchInd(2)+1)];
        plot(thisstart,thisend,'color',clrs{i});
    end
end
% xlim([0 100]);
% ylabel('Power (dB)'); xlabel('Frequency (Hz)');

set(a,'xlim',[0 100]);
set(b,'xlim',[0 100]);
set(c,'xlim',[0 100]);


cTick = [5.5,10.5,22.5,55,80];
set(c,'XTick',cTick,'XTickLabel',{'\theta','\alpha','\beta','low \gamma','high \gamma'},'TickLength',[0 0]);
xlabel(c,'Frequency (Hz)');

bTick = [4,7,8,13,15,30,45,65,70,90];
set(b,'XTick',bTick,'XTickLabel',{});

%vline([4,7,8,13,15,30,45,65,70,90],'-.-k');

ylabel('Power (dB)');
title('Rest and Binge PSD');
legend([ax{1},ax{2}],'Binge','Rest');
%%
for i = 1:5
    relPower(i,5) = mean(relPower(i,1:4)).*100;
    relPower(i,6) = std(relPower(i,1:4)).*100;
end
%% For Core NASL power
vars = {'pwr_NASL_thet','pwr_NASL_alph','pwr_NASL_bet','pwr_NASL_lgam','pwr_NASL_hgam'};
for i = 1:length(vars)
    resp(:,i) = T.Base.(vars{i})(T.Base.coreRespond == 1);
    non(:,i) = T.Base.(vars{i})(T.Base.coreRespond ~= 1);
    respm(1,i) = mean(resp(:,i)); respm(2,i) = std(resp(:,i));
    nonm(1,i) = mean(non(:,i)); nonm(2,i) = std(non(:,i));
end

empty = zeros(1,5);
figure;
barwitherr([empty',respm(2,:)',nonm(2,:)'].*100,[resp(2,:)',respm(1,:)',nonm(1,:)'].*100);
set(gca,'XTickLabel',{'\theta','\alpha','\beta','low \gamma','high \gamma'});
xlabel('Frequency Band'); ylabel('Percent Change in Power from Rest')
title('Percent Change in Binge Power from Rest: NASL');
legend('Sample','Group','Other','Location','southeast');
%%
clrs = {'r','k'};
figure;
fds = {fd1,fd2};
fd1notch = fd1;
% N.B. fds freq axes need to be the same
notchInd = [nearest_idx(59,fds{1}.freq);nearest_idx(61.5,fds{1}.freq)];
c=axes('Position',[.1 .1 .8 1e-12],'Color','none');
set(c,'Units','normalized');

b=axes('Position',[.1 .1 .8 1e-12]);
set(b,'Units','normalized');
set(b,'Color','none');

a=axes('Position',[.1 .2 .8 .7]);
set(a,'Units','normalized');

for i = 1:length(fds)
    temp = nanmean(sq(fds{i}.cohspctrm(6,:,:)),2);
    temp(notchInd(1):notchInd(2)) = NaN;
    h{i} = plot(fds{i}.freq,temp,'color',clrs{i});
    % Fill NaN data in plot
    hold on;
    thisstart = [fds{i}.freq(notchInd(1)-1),fds{i}.freq(notchInd(2)+1)];
    thisend = [temp(notchInd(1)-1),temp(notchInd(2)+1)];
    plot(thisstart,thisend,'color',clrs{i});
end
set(a,'xlim',[0 100]);
set(b,'xlim',[0 100]);
set(c,'xlim',[0 100]);


cTick = [5.5,10.5,22.5,55,80];
set(c,'XTick',cTick,'XTickLabel',{'\theta','\alpha','\beta','low \gamma','high \gamma'},'TickLength',[0 0]);
xlabel(c,'Frequency (Hz)');

bTick = [4,7,8,13,15,30,45,65,70,90];
set(b,'XTick',bTick,'XTickLabel',{});
ylabel('Coherence');
legend([h{1},h{2}],'Binge','Rest');
title('Binge and Rest Coherence: NACL-NACR');
%% Get averages for NACL_NACR coherence (channel pair = 6)
%nonInd = [1,3,6,9,11]; respInd = [2,4,5,7,8,10,12];
vars = {'coh_NACL_NACR_thet','coh_NACL_NACR_alph','coh_NACL_NACR_bet','coh_NACL_NACR_lgam','coh_NACL_NACR_hgam'};
for i = 1:length(vars)
    resp(:,i) = T.Base.(vars{i})(T.Base.shellRespond == 1);
    non(:,i) = T.Base.(vars{i})(T.Base.shellRespond ~= 1);
    respm(1,i) = mean(resp(:,i)); respm(2,i) = std(resp(:,i));
    nonm(1,i) = mean(non(:,i)); nonm(2,i) = std(non(:,i));
end

empty = zeros(1,5);

figure;
barwitherr([empty',respm(2,:)',nonm(2,:)'].*100,[resp(2,:)',respm(1,:)',nonm(1,:)'].*100);
set(gca,'XTickLabel',{'\theta','\alpha','\beta','low \gamma','high \gamma'});
xlabel('Frequency Band'); ylabel('Percent Change in Coherence from Rest')
title('Percent Change in Binge Coherence from Rest: NACL-NACR');
legend('Individual','Group','Other','Location','southeast');