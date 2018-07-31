function [mGC,p,q] = permGC(data,fs,iterates,q)
%%
% INPUTS
% data = raw data; format = channel X time X trial
% fs = sampling frequency
% iterates = number of iterates to use for permutation testing
% q = critical value to use; equates to percentile; i.e. 0.05 will find
%   the value equivalent to a p < 0.05 significance
%__________________________________________________________________________
% OUTPUTS
% mGC = maximum GC values; used to determine cutoff
%%
% Check that q is usable
if rem(q*iterates,1)~=0
    q = round(q*iterates)/iterates;
    warning([num2str(q),'*',num2str(iterates),...
        ' does not equal an integer; rounding up to ',num2str(q),'.'])
end
if q*iterates-1==0
    % If q will lead to 0 indexing determine minimial q.
    q = 2/iterates;
    warning(['Chosen q will lead to 0 indexing; increasing to ',...
        num2str(q),'.']) 
end
% Determine all pairs; flip and duplicate to account for reversed direction
cmbs = [nchoosek(1:size(data,1),2)];%fliplr(nchoosek(1:size(data,1),2))];
for ii = 1:size(cmbs,1)
    disp([num2str(ii), ' of ',num2str(size(cmbs,1))])
    % Keep first channel ordered; permute second
    thisData = [data(cmbs(ii,1),:,:);...
        data(cmbs(ii,2),:,randperm(size(data,3)))];
    % Flip data around for condGCnPar
    thisData = permute(thisData,[2,3,1]);
    for jj = 1:iterates
        [thisGC,~] = condGCnPar(thisData,fs);
        thisMGC = squeeze(max(thisGC,[],1));
        mGC(cmbs(ii,1),cmbs(ii,2),jj) = thisMGC(1,2);
        mGC(cmbs(ii,2),cmbs(ii,1),jj) = thisMGC(2,1);
    end
    % Sort mGC values in descending order
    smGC12 = sort(mGC(cmbs(ii,1),cmbs(ii,2),:),'descend');
    smGC21 = sort(mGC(cmbs(ii,2),cmbs(ii,1),:),'descend');
    % Grab mGC value that is one below q 
    p(cmbs(ii,1),cmbs(ii,2)) = smGC12(iterates*q-1);
    p(cmbs(ii,2),cmbs(ii,1)) = smGC21(iterates*q-1);
end
