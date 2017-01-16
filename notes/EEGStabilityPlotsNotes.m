
%% Get differnce between variables of different response groups
for r = 1:4
    for ii = 1:50
        
    end
end
%% Get NACL theta and NASR-NACR high gamma
for s = 1:length(first)
        theta_sub(s) = table2array(T.Base(second(s),15)) - table2array(T.Base(first(s),15));
        hgamma_sub(s) = table2array(T.Base(second(s),49)) - table2array(T.Base(first(s),49));
end
%% Plot NACL Theta
figure;
plot(day,theta_sub.*100,'ks');
ylim([-2 2]);
lsline;
thetMd = fitlm(day,theta_sub.*100);
text(10,1.5,strcat('p = ',num2str(thetMd.Coefficients.pValue(2),2),' R^2 = ',num2str(thetMd.Rsquared.Ordinary,2)))
title('Stability of NACL Theta Power over Time'); xlabel('Number of Days'); ylabel('Change in Percent Power');
%% Plot NASR-NACR high gamma
figure;
plot(day,hgamma_sub,'ks');
ylim([-.3 .3]);
lsline;
lgammaMd = fitlm(day,hgamma_sub);
text(10,.2,strcat('p = ',num2str(lgammaMd.Coefficients.pValue(2),2),' R^2 = ',num2str(lgammaMd.Rsquared.Ordinary,2)))
title('Stability of NASR-NACR High Gamma Coherence over Time'); xlabel('Number of Days'); ylabel('Change in Coherence Index');
%% Setup power and coherence indices
% Gets indices of first entry of each channel or channel-pair
powInd = [1:5:20];
cohInd = [21:5:50];
%% Model and plot all power per channel
marker = {'ok','sk','dk','+k','xk'};
lineSty = {'-',':','-.','--','-'};
ts = {'NASL','NASR','NACL','NACR'};
figure;
for c = 1:4
    subplot(2,2,c); hold on
    for b = 1:5
        plot(day,sub(:,powInd(c)+(b-1)).*100,marker{b})
        md1{b,c} = fitlm(day,sub(:,powInd(c)+(b-1)).*100);
    end
    ls = lsline;
    for ii = 1:5;
        set(ls(ii),'LineStyle',lineSty{ii});
        if ii == 5
            set(ls(ii),'LineWidth',1.25);
        end
    end
    ylim([-3 3]);
    title(ts{c}); xlabel('Days'); ylabel('% Change in Power');
end
legend('Theta','Alpha','Beta','L Gamma','H Gamma','Theta','Alpha','Beta','L Gamma','H Gamma');
% Create table for p-values and R^2
stat = [];
for m = 1:numel(md1)
    pStat(m,1) = md1{m}.Coefficients.pValue(2);
    pStat(m,2) = md1{m}.Rsquared.Ordinary;
end
%% Model and plot coherence values per channel-pair
marker = {'ok','sk','dk','+k','xk'};
lineSty = {'-',':','-.','--','-'};
ts = {'NASL-NASR','NASL-NACL','NASL-NACR','NASR-NACL','NASR-NACR','NACL-NACR'};
figure;
for c = 1:6
    subplot(2,3,c); hold on
    for b = 1:5
        plot(day,sub(:,cohInd(c)+(b-1)),marker{b})
        md2{b,c} = fitlm(day,sub(:,cohInd(c)+(b-1)).*100);
    end
    ls = lsline;
    for ii = 1:5;
        set(ls(ii),'LineStyle',lineSty{ii});
        if ii == 5
            set(ls(ii),'LineWidth',1.25);
        end
    end
    ylim([-1 1]);
    title(ts{c}); xlabel('Days'); ylabel('Change in Coherence Index');
end
legend('Theta','Alpha','Beta','L Gamma','H Gamma','Theta','Alpha','Beta','L Gamma','H Gamma');

%% Create table for p-values and R^2
for m = 1:numel(md2)
    cStat(m,1) = md2{m}.Coefficients.pValue(2);
    cStat(m,2) = md2{m}.Rsquared.Ordinary;
end
