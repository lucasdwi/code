function ci = conf(x,perc,varargin)
%% Calculates confidence interval
%__________________________________________________________________________
% INPUTS
% x = data; if matrix, will calculate by rows
% perc = confidence percent; format = decimal (e.g. .95 for 95%)
% tail = number of tails to use for calculating t; either 1 or 2
%__________________________________________________________________________
% OUTPUTS
% ci = confidence interval value; N.B. does not give mean, rather the
%   margin
%__________________________________________________________________________
% LLD 2017
%% Check for proper inputs
p = inputParser;
p.CaseSensitive = false;
addRequired(p,'x',@isnumeric);
addRequired(p,'perc');
addParameter(p,'tail',2);
parse(p,x,perc,varargin{:});

assert(perc<1,[num2str(perc),' is not in the correct format. Make sure it is'...
       'a decimal < 1'])
tail = p.Results.tail;
if tail ~= 1 && tail ~= 2
   warning([num2str(tail),' is not a proper input for tails; defaulting'...
       '  to 2-tailed.']) 
   tail = 2;
end
%% Calculate confidence margins
for ii = 1:size(x,1)
    % Get sample size
    n = numel(x(ii,:));
    % Get standard deviation
    s = std(x(ii,:));
    % Get t-value to use
    if tail == 2
        % For two tailed, shift 'perc' up by 1/2 of the difference from 1
        t = tinv(perc+(1-perc)/2,n-1);
    elseif tail == 1
        % Otherwise use 'perc' given
        t = tinv(perc,n-1);
    end
    % Calculate confidence interval
    ci(ii) = t*(s/sqrt(n));
end
