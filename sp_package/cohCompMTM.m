function [coh] = cohCompMTM(trls,adfreq,eoi)
%%
%%
% Get number of channels
chans = size(trls{1,1}.trial,1);
% Get set of all combinations
cmb = nchoosek(1:chans,2);
% Get number of combinations
nCmb = size(cmb,1);
% Get number of events
nE = size(eoi,1);
% If more than one event, preallocate output arrays
coh = cell(1,nE);
% Cycle through each event of interest
for iE = 1:size(eoi,1)
    % Preallocate for event
    % Check notes at bottom for derivation of this size
    thisC = zeros(size(trls{1,iE}.trial,3),size(trls{1,iE}.trial,2)/2,nCmb);
    thisPh = zeros(size(thisC));
    thisPhi = zeros(size(thisC));
    thisCi = zeros(1,nCmb);
    % Cycle through each combination
    for iC = 1:nCmb
        % Cycle through all trials
        for iT = 1:size(trls{1,iE}.trial,3)
            % Grab this x
            x = trls{1,iE}.trial(cmb(iC,1),:,iT);
            % Grab this y
            y = trls{1,iE}.trial(cmb(iC,2),:,iT);
            % Calculate coherence
            [thisF,thisC(iT,:,iC),thisPh(iT,:,iC),dummyCi,thisPhi(iT,:,iC)] = cmtm(x,y,1/adfreq,8,0,0,0);
            thisCi(iC) = dummyCi(1);
        end
    end
    % Plop everything into cell arrays
    coh{1,iE}.f = thisF;
    coh{1,iE}.coh = thisC;
    coh{1,iE}.ph = thisPh;
    coh{1,iE}.ci = thisCi;
    coh{1,iE}.phi = thisPhi;
end

%% Derivation Notes:
% The size used for preallocating the outputs of cmtm.m is from noticing
% that cmtm.m takes into account the Nyquist frequency when constructing
% the frequency output using the given sampling rate (adfreq). As such
% cmtm.m uses the following to create the frequency vector:
% 0:1/(N*1/dt):dt-1/(N*1/dt)
% where dt is the sampling interval (1/adfreq) and N is the length of the
% data. As such this can be simplified to: 0:adfreq/N:1/adfreq-adfreq/N
% which is then cut in half. This means that the length is
% (adfreq/2*N)/adfreq = N/2.
%%
% test = c>=ci(1,1);
% b = bwboundaries(test');
% figure
% imagesc([1:360].*5,s,c')
% hold on
% %  figure
% for bi = 1:length(b)
%     bound = b{bi};
%     hold on
%     plot(bound(:,2).*5,bound(:,1),'k')
% end