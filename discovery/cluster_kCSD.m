%% master script for running event-based kCSD on Discovery cluster
%
% make sure to set data_path correctly
%
% run from command line: matlab -nojvm -nosplash -nodisplay -r cluster_kCSD.m

%% set paths
restoredefaultpath;

arch = getenv('ARCH');
switch arch
    case '' % running on Windows, for testing (?, why doesn't this return anything?)
        addpath(''); % kCSD code
        data_path = '';
    case 'glnxa64' % running on cluster
        addpath(genpath('/ihome/mvdm/code/kCSDv1')); % kCSD code
        addpath(genpath('/ihome/mvdm/code/github/vandermeerlab/code-matlab/shared')); % shared code
        data_path = '/ihome/mvdm/data';
end




%% find files to process

% go into data folder
cd(data_path);

% get list of matching filenames
do_fd = FindFiles('*gammaEvt.mat');
fd_status = nan(size(do_fd));

%% main loop over files 

for iF = 1:length(do_fd)
   
   fn_in = do_fd{iF};
   load(fn_in); % loads a struct called gammaEvt
   
   % check if output variable already exists for this file (DONE)
   if isfield(gammaEvt,'csd')
       fd_status(iF) = 1;
       fprintf('File %s (%d/%d): COMPLETED\n',fn_in,iF,length(do_fd));
   end
   
   % check if job has been submitted; if not, submit job and update file
   switch gammaEvt.submitted
       case 1
           fd_status(iF) = 2;
           fprintf('File %s (%d/%d): previously submitted\n',fn_in,iF,length(do_fd));
           
       case 0 % submit
           fd_status(iF) = 3;
           fprintf('File %s (%d/%d): generating job file...\n',fn_in,iF,length(do_fd));
           
           % call job generator here
           job_fn = cluster_writekCSDjobFile(fn_in);
           
           fprintf('File %s (%d/%d): submitting...\n',fn_in,iF,length(do_fd));
           COMMAND = cat(2,'qsub ',job_fn);
           submit_out = system(COMMAND);
           
           % actual job command here: call matlab function which loads this
           % file... but how to set up environment (modules) correctly on node
           % running it?
           %
           % command is: cluster_runOnekCSD(fn_in);
           % so becomes something like matlab -nojvm -nosplash -nodisplay -r "cluster_runOnekCSD('/ihome/mvdm/data/test_gammaEvt.mat');"
           
           gammaEvt.submitted = 1;
           save(fn_in,'gammaEvt');
   end
   

end
exit;
