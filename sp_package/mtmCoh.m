function c = mtmCoh(x,y,k,E,lambda,N)
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
    % Set up variance term with lambda
    a = vari*(1-lambda);
    % Iterate through until tolerance is reached
    while sum(abs(P-P1)/N)>tol
        % Calculate weights (368a)
        b = (P*ones(1,k))./(P*V'+ones(N,1)*a');
        % New spectral estimate
        wk = (b.^2).*(ones(N,1)*lambda');
        % Use adaptive weights (370a)
        P1 = (sum(wk'.*Pk')./ sum(wk'))';
        % Swap P, P1, and Ptemp
        Ptemp = P1;
        P1 = P;
        P = Ptemp;
    end
    % Adapt DFT for x
    if ii == 1
       fkx = sqrt(k)*sqrt(wk).*fkx./repmat(sum(sqrt(wk'))',1,k);
       Fx = P;
    end
    % Adapt DFT for y
    if ii == 2
        fky = sqrt(k)*sqrt(wk).*fky./repmat(sum(sqrt(wk'))',1,k);
        Fy=P;
    end
end
% Calculate cross spectral density
Cxy = sum([fkx.*conj(fky)]');
% Calculate coherence
c = abs(Cxy)./sqrt(sum(abs(fkx').^2).*sum(abs(fky').^2));