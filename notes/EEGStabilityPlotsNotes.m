%% Setup day information and perform subtractions
% Days between recordings of H10, H13, H14, H15, I11, I12, I1, I2, I3, I4, I6, I8
day = [18,46,16,12,1,13,14,34,2,71,67,36];
% Get every other row index for subtraction
first = [1:2:size(T.Base,1)];
second = [2:2:size(T.Base,1)];
for s = 1:length(first)
    sub(s,:) = table2array(T.Base(second(s),8:end)) - table2array(T.Base(first(s),8:end));
end
absSub = abs(sub);
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
