[data,~,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\data\'...
    'irdm\processed\base\'],{'IRDM11';'IRDM14';'IRDM15';'IRDM16';...
    'IRDM18';'IRDM21';'IRDM22';'IRDM23'},{'pow','coh'},'avg','rel');
save('C:\Users\Pythia\Documents\GreenLab\data\irdm\irdmBaseData.mat',...
    'data','files')
%%
load('C:\Users\Pythia\Documents\GreenLab\data\irdm\irdmBaseData.mat')
% Pre-surgey median split; animals 15,16,21,23 = low (0); 11,14,18,22 =
% high (1)
catData = cat(1,data{:});
x = cat(1,catData{:});
% Pre-surgery
y = [ones(9,1);zeros(14,1);ones(4,1);zeros(4,1);ones(4,1)];
cfg = lassoNetCfg([],[],'n','y','n',100,'1se',[]);
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x,y,...
    'binomial','class',1,5,1,cfg);
save('C:\Users\Pythia\Documents\GreenLab\data\irdm\basePreReal.mat',...
    'catData','x','y','accArray','allAlpha','allBeta','allLambda',...
    'cvFitsArray','hist')
% Post/Thethered; 15,16,18,21 = low (0); 11,14,22,23 = high (1); based on
% average of last 5 tethered sessions (used post-surgery value for 23)
y = [ones(9,1);zeros(18,1);ones(8,1)];
cfg = lassoNetCfg([],[],'n','y','n',100,'1se',[]);
[allAlpha,allLambda,allBeta,cvFitsArray,accArray,hist] = lassoNet(x,y,...
    'binomial','class',1,5,1,cfg);
save('C:\Users\Pythia\Documents\GreenLab\data\irdm\basePostReal.mat',...
    'catData','x','y','accArray','allAlpha','allBeta','allLambda',...
    'cvFitsArray','hist')
%% Pre-Surgery Prediction - Using even split permuted
permErr = [];
for ii = 1:18
    load(['C:\Users\Pythia\Documents\GreenLab\data\irdm\basePerm\evenPre\'...
        'irdmBasePerm',num2str(ii),'.mat']) 
    permErr = [permErr;allRndLambda{1}.allErr];
end
load('C:\Users\Pythia\Documents\GreenLab\data\irdm\baseReal.mat')
doubleHist((1-allLambda{1}.allErr).*100,(1-permErr).*100,'main',...
    'Pre-Surgery','xlab','Accuracy')
%% Post-Surgery Prediction - Using even split permuted
permErr = [];
for ii = 1:18
    load(['C:\Users\Pythia\Documents\GreenLab\data\irdm\basePerm\'...
        'evenPost\irdmPostBasePerm',num2str(ii),'.mat'])
    permErr = [permErr;allRndLambda{1}.allErr];
end
load('C:\Users\Pythia\Documents\GreenLab\data\irdm\basePostReal.mat',...
    'allLambda')
doubleHist((1-allLambda{1}.allErr).*100,(1-permErr).*100,'main',...
    'Post-Surgery','xlab','Accuracy')
%%
[data,~,files] = collateData(['C:\Users\Pythia\Documents\GreenLab\data'...
    '\irdm\processed\ddt\'],{'IRDM11';'IRDM14';'IRDM15';'IRDM16';...
    'IRDM18';'IRDM21';'IRDM22'},{'pow','coh'},'avg','rel');
allZ = cell(1,size(data,2));
for ii = 1:size(data,2)
    thisCat = cat(1,data{1,ii}{:});
    thisZ = zscore(thisCat);
    
    allZ{ii} = thisZ;
end

%% timings
eventTs.label = {'levIM','levDEL','head','feeder','ltDEL','ltIM','tone',...
    'TrialInit'};
% Use TrialInit as start time for trial
for ii = 1:size(eventTs.t{8},1)
    % Determine if one or three pellets were rewarded
    fInd = nearest_idx3(eventTs.t{1,8}(ii),eventTs.t{1,4},1);
    % First check if last index
    if fInd == size(eventTs.t{1,4},1)
        trial{ii}.reward = 'im';
        trial{ii}.leverPress = eventTs.t{1,1}(...
                nearest_idx3(eventTs.t{1,8}(ii),eventTs.t{1,1}));
    else
        if eventTs.t{1,4}(fInd+1)-eventTs.t{1,4}(fInd) < 1
            trial{ii}.reward = 'del';
            trial{ii}.leverPress = eventTs.t{1,2}(...
                nearest_idx3(eventTs.t{1,8}(ii),eventTs.t{1,2}));
        else
            trial{ii}.reward = 'im';
            trial{ii}.leverPress = eventTs.t{1,1}(...
                nearest_idx3(eventTs.t{1,8}(ii),eventTs.t{1,1}));
        end
    end
    trial{ii}.toneStart = eventTs.t{1,7}(nearest_idx3(eventTs.t{1,8}(ii)...
        ,eventTs.t{1,7},-1));
    trial{ii}.headStart = eventTs.t{1,8}(ii);
    trial{ii}.feeder = eventTs.t{1,4}(nearest_idx3(eventTs.t{1,8}(ii),...
        eventTs.t{1,4}));
    trial{ii}.headReward = eventTs.t{1,3}(nearest_idx3(...
        eventTs.t{1,8}(ii),eventTs.t{1,3}));
end
%%
    % Determine which lever was pressed first (Immediate [1] vs. Delay [2]
    im = eventTs.t{1,1}(nearest_idx3(eventTs.t{1,8}(ii),...
        eventTs.t{1,1},1))-evenTs.t{1,8}(ii);
    de = eventTs.t{1,2}(nearest_idx3(eventTs.t{1,8}(ii),...
        eventTs.t{1,2},1))-eventTs.t{1,8}(ii);
    if im < de || (im>0 && de<0)
        reward{ii} = 'im';
%        this.reward = 'im';
    else
        reward{ii} = 'de';
%         this.reward = 'de';
    end
    % Calculate
end