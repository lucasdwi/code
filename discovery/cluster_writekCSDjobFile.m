function fn_out = cluster_writekCSDjobFile(fn_in)
% write PBS job file for clusterified kCSD

cfg = [];
cfg.q = 'testing'; % default, testing, largeq
cfg.nnodes = 1;
cfg.walltime = '02:00:00'; % h-m-s
cfg.codepath = '/ihome/mvdm/code/kCSDv1';

% figure out correct filenames
[fp,fn,fe] = fileparts(fn_in); % fn_in is the FULL path to .mat filename
arch = getenv('ARCH');
switch arch
    case '' % running on Windows, for testing (?, why doesn't this return anything?)
        filesep = '\';
    case 'glnxa64' % running on cluster
        filesep = '/';
end

if ~isempty(fp)
    fn_out = cat(2,fp,filesep,fn,'.pbs');
else
    fn_out = cat(2,fn,'.pbs');
end
cfg.jobname = fn;

fh = fopen(fn_out,'w');

% define what to write in job file (clunky, should just get text template)
lines{1} = '#!/bin/bash -l';
lines{2} = '# declare a name for this job';
lines{3} = sprintf('# PBS -N %s',cfg.jobname);
lines{4} = '# request the queue (enter the possible names, if omitted, default is the default)';
lines{5} = '# if more then 600 jobs use the largeq';
lines{6} = sprintf('# PBS -q %s',cfg.q);
lines{7} = '# request 1 core on 1 node';
lines{8} = '# ensure you reserve enough cores for the projected memory usage';
lines{9} = '# figuring 4G/core';
lines{10} = sprintf('# PBS -l nodes=%d:ppn=1',cfg.nnodes);
lines{11} = '# request 4 hours and 30 minutes of wall time';
lines{12} = sprintf('# PBS -l walltime=%s',cfg.walltime);
lines{13} = '# mail is sent to you when the job begins and when it exits or aborts';
lines{14} = '# you can use all or some or none. If you don''t want email leave this';
lines{15} = '# and the following (#PBS -M) out of the script.';
lines{16} = '#'; % '#PBS -m bea';
lines{17} = '# specify your email address';
lines{18} = '#'; % '#PBS -M John.Smith@dartmouth.edu';
lines{19} = '# By default, PBS scripts execute in your home directory, not the';
lines{20} = '# directory from which they were submitted. The following line';
lines{21} = '# places you in the directory from which the job was submitted.';
lines{22} = 'cd $PBS_O_WORKDIR';
lines{23} = '# run the program';
lines{24} = 'module add matlab/r2014b';
lines{25} = sprintf('cd %s',cfg.codepath);
lines{26} = sprintf('matlab -nojvm -nosplash -nodisplay -r "cluster_runOnekCSD(''%s'');"',fn_in);

for iL = 1:length(lines)
   
    fprintf(fh,cat(2,lines{iL},'\n'));
    
end

fclose(fh);



 