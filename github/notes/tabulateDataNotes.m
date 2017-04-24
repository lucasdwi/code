%% Notes


%% Calculate percent change from baseline to food dep in power and coherence
for i = 1:length(files.Base)
    PSDs.DepvsBase{i} = (PSDs.Dep.Bands{i,3}-PSDs.Base.Bands{i,3})./abs(PSDs.Base.Bands{i,3});
    Cohs.DepvsBase{i} = (Cohs.rel.Dep{i}-Cohs.rel.Base{i})./abs(Cohs.rel.Base{i});
end
%% Prep channel-band data for table
extract = @(C,k) cellfun(@(c)c(k),C);
for i = 1:numel(PSDs.DepvsBase{1})
    cb{i} = extract(PSDs.DepvsBase,(i));
end
cbNames = {'pwr_NASLthet','pwr_NASLalph','pwr_NASLbet','pwr_NASLlgam','pwr_NASLhgam','pwr_NASRthet','pwr_NASRalph','pwr_NASRbet','pwr_NASRlgam','pwr_NASRhgam','pwr_NACLthet','pwr_NACLalph','pwr_NACLbet','pwr_NACLlgam','pwr_NACLhgam','pwr_NACRthet','pwr_NACRalph','pwr_NACRbet','pwr_NACRlgam','pwr_NACRhgam'};
%cbNames = {'pwr_1thet','pwr_1alph','pwr_1bet','pwr_1lgam','pwr_1hgam','pwr_2thet','pwr_2alph','pwr_2bet','pwr_2lgam','pwr_2hgam','pwr_3thet','pwr_3alph','pwr_3bet','pwr_3lgam','pwr_3hgam','pwr_4thet','pwr_4alph','pwr_4bet','pwr_4lgam','pwr_4hgam'};
%% Prep channelpair-band data for table
for i = 1:numel(Cohs.DepvsBase{1})
    pb{i} = extract(Cohs.DepvsBase,(i));
end
%% Create names
chlComb = {'NASL_NASR'; 'NASL_NACL'; 'NASL_NACR'; 'NASR_NACL';'NASR_NACR';'NACL_NACR'};
bands = {'thet','alph','bet','lgam','hgam'};
pbNames = {};
for j = 1:length(chlComb)
    for k = 1:length(bands)
        thisName = strcat(chlComb{j},'_',bands{k});
        pbNames = horzcat(pbNames,thisName);
    end
end
%% Create data table
% Initialize with animal and group name
% 1 = neither, 2 = core, 3 = shell [N.B. numbers only used to put neither at top for modeling - makes neither the default]
resp = {'2';'3';'2';'3';'1';'1';'3';'3';'2';'1';'1';'2'};
T = cell2table(resp,'VariableNames',{'group'},'RowNames',{'H10','H13','H14','H15','I1','I2','I3','I4','I6','I8','I11','I12'}); 
% Fill in power data (channel-frequency band)
for i = 1:numel(PSDs.DepvsBase{1})
    this = array2table(cb{i}','VariableNames',{cbNames{i}});
    T = horzcat(T,this);
end
% Fill in coherence data (channel combo-frequncy band)
for i = 1:numel(Cohs.DepvsBase{1})
    this = array2table(pb{i}','VariableNames',{pbNames{i}});
    T = horzcat(T,this);
end
% Add binge response columns
shellReduct = [-0.069767442,-0.4140625,0.06993007,-0.564343164,-0.499866986,-0.20657277,-0.4353683,-0.452188799,-0.060344828,-0.593134139,0.016771488,-0.350282486];
coreReduct = [-0.444861215,-0.15625,-0.27972028,0.105898123,-0.396212933,-0.241622575,-0.073514602,0.093333333,-0.340761374,-0.568500539,0.228710462,-0.249753208];
T = horzcat(T,table(coreReduct','VariableNames',{'coreReduct'}),table(shellReduct','VariableNames',{'shellReduct'}));
% Sort table by group
T = sortrows(T,'group','ascend');
% Give categories names
T.group = categorical(T.group,{'1','2','3'},{'neither','core','shell'});
%% Define predictor and response variables
allPredict = horzcat(cbNames,pbNames);
allResponse = {'shellReduct','coreReduct'};