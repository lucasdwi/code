function [r,c] = subDim(X,maxCol,dim)
%% Provide number of subplots, max number of subplots per column, and preferred larger dimension 
% dim = x or y; if x then will opt for wider array over taller if
% necessary
% Calculate possible numbers of rows given number of overall plots and columns
% Create vector of all possible numbers of columns above 1
nCs = 2:maxCol;
possR = [];
for i = 1:length(nCs)
    possR(i) = X/nCs(i);
end

possInt = possR(rem(possR,1)==0);
% Check if there are no possible full subplot arrays
if isempty(possInt)
    if dim == 'y'
        % Calculate and minimize remainders
        remainX = rem(possR,1);
        maxRem = max(remainX);
        r = floor(possR(remainX == maxRem)); % Might have problem if more than one option has the same remainder
        c = ceil(X/r);
    else if dim == 'x'
        remainX = rem(possR,1);
        minRem = min(remainX);
        r = floor(possR(remainX == minRem)); % Might have problem if more than one option has the same remainder
        c = ceil(X/r);
        end
    end
end
numX = (nCs+1).*floor(possR)-13;
numY = nCs.*ceil(possR)-13;



