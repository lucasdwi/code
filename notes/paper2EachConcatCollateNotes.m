%% Load concatLog files
clear
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\',...
    'concatLog2\'])
% Preallocate
[ccLogA] = deal(zeros(1,20));
ceLogA = zeros(20,12);
[ccLogX,ccLogY] = deal(cell(1,20));
[ceLogX,ceLogY] = deal(cell(20,12));
for ii = 1:20
   load([num2str(ii),'.mat'])
   % Extract concat->concat data
   ccLogA(ii) = concatData.auc;
   ccLogX{ii} = concatData.rocX;
   ccLogY{ii} = concatData.rocY;
   % Extract concat->each data
   ceLogA(ii,:) = eachData.auc;
   ceLogX(ii,:) = eachData.rocX;
   ceLogY(ii,:) = eachData.rocY;
end
% Preallocate
[ceLogXM,ceLogYM,ceLogXS,ceLogYS] = deal(zeros(12,size(ceLogX{1},1)));
% Get average ROC for 'each'
for ii = 1:12
    ceLogXM(ii,:) = mean(cat(2,ceLogX{:,ii}),2);
    ceLogYM(ii,:) = mean(cat(2,ceLogY{:,ii}),2);
    ceLogXS(ii,:) = std(cat(2,ceLogX{:,ii}),[],2);
    ceLogYS(ii,:) = std(cat(2,ceLogY{:,ii}),[],2);
end
% Get average ROC for concat-concat models
ccLogXM = mean(cat(2,ccLogX{:}),2);
ccLogYM = mean(cat(2,ccLogY{:}),2);
% Get average 'fill' for concat-concat models
[ccLogXfill,ccLogYfill] = avgFill(cat(2,ccLogX{:}),cat(2,ccLogY{:}),2,1);
% Get average auc and 95% CI using 19 dof
ccLogAM = mean(ccLogA);
ceLogAM = mean(ceLogA,1);
ccLogACI = conf(ccLogA,0.95,'tail',2);
ceLogACI = conf(ceLogA',0.95,'tail',2);
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLog.mat'],'ccLogA','ccLogX','ccLogY','ccLogAM','ccLogACI',...
    'ccLogXM','ccLogYM','ceLogA','ceLogX','ceLogY','ceLogAM','ceLogACI',...
    'ceLogXM','ceLogYM','ccLogXfill','ccLogYfill')
%% Load concat files
clear
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\',...
    'concat\'])
% Preallocate
[ccA] = deal(zeros(1,20));
ceA = zeros(20,12);
[ccX,ccY] = deal(cell(1,20));
[ceX,ceY] = deal(cell(20,12));
betas = zeros(20,58);
for ii = 1:20
   load([num2str(ii),'.mat'])
   % Extract concat->concat data
   ccA(ii) = concatData.auc;
   ccX{ii} = concatData.rocX;
   ccY{ii} = concatData.rocY;
   % Extract concat->each data
   ceA(ii,:) = eachData.auc;
   ceX(ii,:) = eachData.rocX;
   ceY(ii,:) = eachData.rocY;
   % Extract betas
   betas(ii,:) = concatData.allBeta{1,1}.survBeta;
end
% Preallocate
[ceXM,ceYM,ceXS,ceYS] = deal(zeros(12,size(ceX{1},1)));
% Get average ROC for 'each'
for ii = 1:12
    ceXM(ii,:) = mean(cat(2,ceX{:,ii}),2);
    ceYM(ii,:) = mean(cat(2,ceY{:,ii}),2);
    ceXS(ii,:) = std(cat(2,ceX{:,ii}),[],2);
    ceYS(ii,:) = std(cat(2,ceY{:,ii}),[],2);
end
% Get average ROC for concat-concat models
ccXM = mean(cat(2,ccX{:}),2);
ccYM = mean(cat(2,ccY{:}),2);
% Get average 'fill' for concat-concat models
[ccXfill,ccYfill] = avgFill(cat(2,ccX{:}),cat(2,ccY{:}),2,1);
% Get average auc and 95% CI 
ccAM = mean(ccA);
ceAM = mean(ceA,1);
ccACI = conf(ccA,0.95,'tail',2);
ceACI = conf(ceA',0.95,'tail',2);
% Get average number of features
featM = mean(sum(betas~=0,2));
% Get average feature survival
betaM = mean(betas,1)';
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLasso.mat'],'ccA','ccX','ccY','ccAM','ccACI','ccXM','ccYM',...
    'ceA','ceX','ceY','ceAM','ceACI','ceXM','ceYM','ccXfill','ccYfill',...
    'featM','betaM','betas')
%% Load concatLogRand files
clear
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\',...
    'concatLogRand2\'])
% Preallocate
ccLogRandA = zeros(1,20);
ceLogRandA = zeros(20,12);
[ccLogRandX,ccLogRandY] = deal(cell(1,20));
[ceLogRandX,ceLogRandY] = deal(cell(20,12));
for ii = 1:20
   load([num2str(ii),'.mat'])
   % Check if ROC is 50 line, if so interpolate it; store curve
   if isequal(concatData.rocX,[0;1])
      ccLogRandX{ii} = (0:1/1200:1)'; 
      ccLogRandY{ii} = (0:1/1200:1)';
   else
       ccLogRandX{ii} = concatData.rocX;
       ccLogRandY{ii} = concatData.rocY;
   end
   % Check if ROC is 50 line, if so interpolate it; store curve
   for jj = 1:12
       if isequal(eachData.rocX{jj},[0;1])
           ceLogRandX{ii,jj} = (0:1/125:1)';
           ceLogRandY{ii,jj} = (0:1/125:1)';
       else
           ceLogRandX{ii,jj} = eachData.rocX;
           ceLogRandY{ii,jj} = eachData.rocY;
       end
   end
   % Store auc value
    ccLogRandA(ii) = concatData.auc;
    ceLogRandA(ii,:) = eachData.auc;
end
% Get average ROC for concat-concat random models
ccLogRandXM = mean(cat(2,ccLogRandX{:}),2);
ccLogRandYM = mean(cat(2,ccLogRandY{:}),2);
% Get average 'fill' for concat-concat random models
[ccLogRandXfill,ccLogRandYfill] = avgFill(cat(2,ccLogRandX{:}),...
    cat(2,ccLogRandY{:}),2,1);
% Get average auc and 95% CI using 19 dof
ccLogRandAM = mean(ccLogRandA);
ceLogRandAM = mean(ceLogRandA,1);
ccLogRandACI = conf(ccLogRandA,0.95,'tail',2);
ceLogRandACI = conf(ceLogRandA',0.95,'tail',2);
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLogRand.mat'],'ccLogRandA','ccLogRandX','ccLogRandY',...
    'ccLogRandAM','ccLogRandACI','ccLogRandXM','ccLogRandYM',...
    'ceLogRandA','ceLogRandX','ceLogRandY','ceLogRandAM','ceLogRandACI',...
    'ccLogRandXfill','ccLogRandYfill')
%% Load concatRand files
clear
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\',...
    'concatRand\'])
% Preallocate
ccRandA = zeros(1,20);
ceRandA = zeros(20,12);
[ccRandX,ccRandY] = deal(cell(1,20));
[ceRandX,ceRandY] = deal(cell(20,12));
for ii = 1:20
   load([num2str(ii),'.mat'])
   % Check if ROC is 50 line, if so interpolate it; store curve
   if isequal(concatData.rocX,[0;1])
      ccRandX{ii} = (0:1/1200:1)'; 
      ccRandY{ii} = (0:1/1200:1)';
   else
       ccRandX{ii} = concatData.rocX;
       ccRandY{ii} = concatData.rocY;
   end
   % Check if ROC is 50 line, if so interpolate it; store curve
   for jj = 1:12
       if isequal(eachData.rocX{jj},[0;1])
           ceRandX{ii,jj} = (0:1/125:1)';
           ceRandY{ii,jj} = (0:1/125:1)';
       else
           ceRandX{ii,jj} = eachData.rocX;
           ceRandY{ii,jj} = eachData.rocY;
       end
   end
   % Store auc value
    ccRandA(ii) = concatData.auc;
    ceRandA(ii,:) = eachData.auc;
end
% Get average ROC for concat-concat random models
ccRandXM = mean(cat(2,ccRandX{:}),2);
ccRandYM = mean(cat(2,ccRandY{:}),2);
% Get average 'fill' for concat-concat random models
[ccRandXfill,ccRandYfill] = avgFill(cat(2,ccRandX{:}),...
    cat(2,ccRandY{:}),2,1);
% Get average auc and 95% CI using 19 dof
ccRandAM = mean(ccRandA);
ceRandAM = mean(ceRandA,1);
ccRandACI = conf(ccRandA,0.95,'tail',2);
ceRandACI = conf(ceRandA',0.95,'tail',2);
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\concatLassoRand.mat'],'ccRandA','ccRandX','ccRandY','ccRandAM',...
    'ccRandACI','ccRandXM','ccRandYM','ceRandA','ceRandX','ceRandY',...
    'ceRandAM','ceRandACI','ccRandXfill','ccRandYfill')
%% Load eachLog files
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\',...
    'eachLog\'])
% Preallocate
eeLogA = zeros(12,12,20);
[eeLogX,eeLogY] = deal(cell(1,12));
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    % Get each and self auc and rocs
    eachData.auc(animal) = selfData.auc;
    eeLogA(animal,:,iter) = eachData.auc;
    eachData.rocX{animal} = selfData.rocX;
    eeLogX{animal}(:,:,iter) = cat(2,eachData.rocX{:});
    eachData.rocY{animal} = selfData.rocY;
    eeLogY{animal}(:,:,iter) = cat(2,eachData.rocY{:});
    % Get concat auc and rocs
%     ecLogA(animal,iter) = concatData.auc;
%     ecLogX(animal,:,iter) = concatData.rocX;
%     ecLogY(animal,:,iter) = concatData.rocY;
end
eeSelfLogA = diag(mean(eeLogA,3));
% Get average self auc and 95% CI using 11 dof
eeSelfLogAM = mean(eeSelfLogA);
eeSelfLogACI = tinv(0.975,11)*(std(diag(mean(eeLogA,3)))/sqrt(11));
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\eachLog.mat'],'eeLogA','eeLogX','eeLogY','eeSelfLogAM',...
    'eeSelfLogACI','eeSelfLogA')
%% Load each files
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\',...
    'each\'])
% Preallocate
eeA = zeros(12,12,20);
[eeX,eeY] = deal(cell(1,12));
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    % Get each and self auc and rocs
    eachData.auc(animal) = selfData.auc;
    eeA(animal,:,iter) = eachData.auc;
    eachData.rocX{animal} = selfData.rocX;
    eeX{animal}(:,:,iter) = cat(2,eachData.rocX{:});
    eachData.rocY{animal} = selfData.rocY;
    eeY{animal}(:,:,iter) = cat(2,eachData.rocY{:});
    % Get concat auc and rocs
%     ecA(animal,iter) = concatData.auc;
%     ecX(animal,:,iter) = concatData.rocX;
%     ecY(animal,:,iter) = concatData.rocY;
end
eeSelfA = diag(mean(eeA,3));
% Get average self auc and 95% CI using 11 dof
eeSelfAM = mean(eeSelfA);
eeSelfACI = tinv(0.975,11)*(std(diag(mean(eeA,3)))/sqrt(11));
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\eachLasso.mat'],'eeA','eeX','eeY','eeSelfAM','eeSelfACI','eeSelfA')
%% Load eachLogRand files 
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\',...
    'eachLogRand\'])
% Preallocate
[ecLogRandX,ecLogRandY,eeLogRandX,eeLogRandY] = deal(cell(12,20));
eeLogRandA = zeros(12,12,20);
ecLogRandA = zeros(12,20);
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    % Add self data to each
    eachData.rocX{animal} = selfData.rocX;
    eachData.rocY{animal} = selfData.rocY;
    % Check if roc is 50 line, if so interpolate it; store curve
    if isequal(concatData.rocX,[0;1])
        ecLogRandX{animal,iter} = (0:1/1500:1)';
        ecLogRandY{animal,iter} = (0:1/1500:1)';
    else
        ecLogRandX{animal,iter} = concatData.rocX;
        ecLogRandY{animal,iter} = concatData.rocY;
    end
    % Check if roc is 50 line, if so interpolate it; store curve
    if isequal(eachData.rocX{animal},[0,1])
        eeLogRandX{animal,iter} = (0:1/125:1)';
        eeLogRandY{animal,iter} = (0:1/125:1)';
    else
        eeLogRandX{animal,iter} = eachData.rocX;
        eeLogRandY{animal,iter} = eachData.rocY;
    end
    % Get each and self AUC
    eachData.auc(animal) = selfData.auc;
    eeLogRandA(animal,:,iter) = eachData.auc;
    % Get concat AUC
    ecLogRandA(animal,iter) = concatData.auc;
end
% Get average self AUC and CI
eeSelfLogRandA = diag(mean(eeLogRandA,3));
eeSelfLogRandAM = mean(eeSelfLogRandA);
eeSelfLogRandACI = conf(eeSelfLogRandA',0.95,'tail',2);
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\eachLogRand.mat'],'ecLogRandX','ecLogRandY','eeLogRandX',...
    'eeLogRandY','eeLogRandA','ecLogRandA','eeSelfLogRandA',...
    'eeSelfLogRandAM','eeSelfLogRandACI')
%% Load eachRand files 
cd(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew\',...
    'eachRand\'])
% Preallocate
[ecRandX,ecRandY,eeRandX,eeRandY] = deal(cell(12,20));
eeRandA = zeros(12,12,20);
ecRandA = zeros(12,20);
for ii = 1:240
    load([num2str(ii),'.mat'])
    animal = ceil(ii/20);
    iter = rem(ii,20);
    if iter == 0
        iter = 20;
    end
    % Add self data to each
    eachData.rocX{animal} = selfData.rocX;
    eachData.rocY{animal} = selfData.rocY;
    % Check if roc is 50 line, if so interpolate it; store curve
    if isequal(concatData.rocX,[0;1])
        ecRandX{animal,iter} = (0:1/1500:1)';
        ecRandY{animal,iter} = (0:1/1500:1)';
    else
        ecRandX{animal,iter} = concatData.rocX;
        ecRandY{animal,iter} = concatData.rocY;
    end
    % Check if roc is 50 line, if so interpolate it; store curve
    if isequal(eachData.rocX{animal},[0,1])
        eeRandX{animal,iter} = (0:1/125:1)';
        eeRandY{animal,iter} = (0:1/125:1)';
    else
        eeRandX{animal,iter} = eachData.rocX;
        eeRandY{animal,iter} = eachData.rocY;
    end
    % Get each and self AUC
    eachData.auc(animal) = selfData.auc;
    eeRandA(animal,:,iter) = eachData.auc;
    % Get concat AUC
    ecRandA(animal,iter) = concatData.auc;
end
eeSelfRandA = diag(mean(eeRandA,3));
save(['C:\Users\Pythia\Documents\GreenLab\data\paper2\analyzed\finalNew'...
    '\eachLassoRand.mat'],'ecRandX','ecRandY','eeRandX','eeRandY',...
    'eeRandA','ecRandA','eeSelfRandA')