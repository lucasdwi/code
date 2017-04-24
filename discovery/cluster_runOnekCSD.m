function cluster_runOnekCSD(fn_in)
% function cluster_runOnekCSD(fn_in)

fprintf('cluster_runOnekCSD called, fn_in is %s...\n',fn_in);

%% set paths
restoredefaultpath;

arch = getenv('ARCH');
switch arch
    case '' % running on Windows, for testing (?, why doesn't this return anything?)
        addpath(''); % kCSD code
        data_path = '';
    case 'glnxa64' % running on cluster
        addpath(genpath('/ihome/mvdm/code/kCSDv1')); % kCSD code
        addpath(genpath('/ihome/mvdm/code/github/vandermeerlab/code-matlab/shared')); % shared code; should become git-ified
        data_path = '/ihome/mvdm/data';
end

cd(data_path);

try
    load(fn_in); % loads variable (struct) called gammaEvt
catch
   fprintf('Failed to load file %s.\n',fn_in); 
end

% actually compute CSD here
gammaEvt.CSD = [];

save(fn_in);
exit;
