function [out] = threshRound(in,thresh,target)
%% Rounds numbers with custom threshold from in to targets with threshold
% Uses following logic: if in(ii) >= threshold it is rounded to max(target)
out = zeros(1,numel(in));
for ii = 1:numel(in)
   if in(ii) >= thresh
       out(ii) = max(target);
   else
       out(ii) = min(target);
   end
end