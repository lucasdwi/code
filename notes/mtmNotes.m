function cohCompMTM(trls,nWin,dt,)
nCmbs = size(cmbs,1);
for ei = 1:size(eoi,1)
    nTrls = size(trls{ei}.trial,3);
    % Normalize data s.t. mean(x) = 0
    normData = trls{ei}.trial - mean(trls{ei}.trial,2);
    % Get number of samples in each trial
    nSamp = size(trls{ei}.trial,2);
%     % Calculate parameter w 
%     w = nWin/nSamp;
    % Calculate parameter k s.t. k < the Shannon number (2NW); 
    % optimizes the set of k eigenspectra with good leakage properties    
    k = 1-(2*nWin);
    % Construct set of discrete prolate slepian sequences
    [E,lambda] = dpss(nSamp,nWin,k);
    for ti = 1:nTrls
        for ci = 1:nCmbs
        
        end
    end
end

function mtmCoh(x,y,k,E,lambda,N)
% Replicate x and y s.t. there are k repeats along 2nd dimension
x = x(:,ones(1,k));
y = y(:,ones(1,k));
% Calculate DFTs
fkx = fft(E(:,1:k).*x,N);
fky = fft(E(:,1:k).*y,N);
% Calculate PSD
Pkx = abs(fkx).^2;
Pky = abs(fky).^2;
for ii = 1:2
    % Calculations for x
    if ii == 1
        vari = x'*x/N;
        Pk = Pkx;
    end
    % Calculations for y
    if ii == 2
        vari = y'*y/N;
        Pk = Pky;
    end
    % Set up initial spectral estimate and placeholders
    P = (Pk(:,1)+Pk(:,2))/2;
    Ptemp = zeros(N,1);
    P1 = zeros(N,1);
    % Calculate tolerance
    tol = 0.0005*vari/N;
    % Iterate through until tolerance is reached
    while sum(abs(P-P1)/N)>tol
        % Calculate weights (368a)
        b = P./(lambda*P+(1-lambda)*xVar);
        % Use adaptive weights (370a)
        P1 = (sum((b.^2).*(ones(N,1)*lambda')'.*Pk')./...
            sum((b.^2).*(ones(N,1)*lambda')))';
        % Swap P, P1, and Ptemp
        Ptemp = P1;
        P1 = P;
        P = Ptemp;
    end
    % Adapt DFT for x
    if ii == 1
       fkx = sqrt(k)*sqrt(b.^2).*(ones(N,1)*lambda').*fkx./...
           repmat(sum(sqrt((b.^2).*(ones(N,1)*lambda')'))',1,k); 
    end
    % Adapt DFT for y
    if ii == 2
        fky = sqrt(k)*sqrt(b.^2).*(ones(N,1)*lambda').*fky./...
           repmat(sum(sqrt((b.^2).*(ones(N,1)*lambda')'))',1,k); 
    end
end
% Calculate cross spectral density
Cxy = sum([fkx.*conj(fky)]');
% Calculate coherence
c = abs(Cxy)./sqrt(sum(abs(fkx').^2).*sum(abs(fky').^2));
end