%FUNCTION: Computes conditional Granger causality using non-parametric
%analysis. For each pair-wise connection, it factors out all possible spurious
%connections arising from other nodes in the data set. This
%method removes 1) Zero-lag interactions between nodes arising from
%observed and unobserved sources 2) Lagged interactions arising from common
%drive originating at recorded nodes. Note that this program does NOT
%control for spurious interactions due to lagged influence from unobserved
%sources. Spectra are computed using the multitaper method with NW=3. The
%last taper is discarded and this version does not use a weighted average.
%The single-trial spectra are the direct averages of the first 5
%eigenspectra.
%
%INPUT: X is a times by trials by channels array
%
%       fs is the sampling rate in Hz
%
%OUTPUT: GC is a frequencies by channels by channels array of Granger spectra
%
%        pow contains power spectra
%
%        coh contains coherence spectra
%
%        f is the frequency vector determined by the sampling rate and
%        length of the trials
%
%EXAMPLE:
%
%If Y is a continous recording with time in rows and channels in columns
%and each trial has the same length use:
%
%   for ind = 1:size(Y,2)
%       X(:,:,ind) = reshape(Y(:,ind),lenTrials,numTrials)
%   end
%   [GC,pow,coh,f] = condGCnPar(X,fs)
%
%   plot(f,GC(:,1,2)) plots the Granger spectrum from channel 1 to channel 2
%
%   plot(f,pow(1,:)) plots the PSD of the first channel
%
%   plot(f,coh(:,1,2)) plots the coherence from channel 1 to channel 2
%
%A. Nakhnikian, Jan 2018. Uses subroutines provided by M. Dhamala (see
%second reference)
%
%Refs:
%
%Wen, X. Rangarajan, G. Ding, M. (2013) "Multivariate Granger causality: an
%estimation framework based on factorization of the spectral density
%matrix" Phil. Trans. Royal Academy, 371 (1997)
%
%Dhamala, M. Rangarajan, G. Ding, M. (2008) "Estimating Granger causality
%from Fourier and wavelet transforms of time series data" Phys. Rev. Lett.
%100
%
%Other notes: Dhamala's code supports increased frequency resolution using
%zero padding. To change the frequency resolution, add a 3rd input to
%"sig2mtSpect"

%----------MAIN FUNCTION------------------
function [GC,pow,coh,f] = condGCnPar(X,fs)

Nc = size(X,3);

[S,f]= sig2mTspect_nv(X,fs); %freq res can be 3rd input
nFreq = length(f);
pow = zeros(Nc,nFreq);
%Extract power
for chInd = 1:Nc
    pow(chInd,:) = squeeze(S(chInd,chInd,:));
end
%Compute Coherence
spectra = permute(S,[3 1 2]);
coh = S2coh(spectra);

GC = zeros(Nc,Nc,nFreq);
%Compute GC
chans = 1:Nc;
for ii = 1:Nc %drive is in j->i direction
    for jj = 1:Nc
        if ii==jj, continue, end
        W = setdiff(chans,[ii,jj]); %channel indices excluding current pair
        
        %Rearrange S so that the first 2x2 diag block contains i and j
        Stemp = zeros(size(S));
        for fInd = 1:nFreq
            Stemp(:,:,fInd) = [[S(ii,ii,fInd),S(ii,jj,fInd),S(ii,W,fInd)];...
                [S(jj,ii,fInd),S(jj,jj,fInd),S(jj,W,fInd)];...
                [S(W,ii,fInd),S(W,jj,fInd),S(W,W,fInd)]];
        end
        
        %Compute and normalize Full Transfer Matrix
        W2 = 3:Nc; %channels below i and j
        [H,Z,~] = sfactorization_wilson(Stemp,fs,f); %Full transfer mat and noise cov
        Ct = Z(W2,1:2); sigm = Z(1:2,1:2); %Geweke's original partition of Z
        P = [[1,zeros(1,Nc-1)];[-Z(2,1)/Z(1,1),1,zeros(1,Nc-2)];...
            [-Ct/sigm,eye(Nc-2)]]; %Geweke transformation for j->i|W
        %-This is a condensed form of Eq. 2.31 in Wen et al. See pg 305
        %of Geweke (1982) J. Amer. Stat. Assoc. 77(378)
        Hn = zeros(size(H));
        for fInd = 1:nFreq
            Hn(:,:,fInd) = H(:,:,fInd)/P; %Full transfer matrix with uncorrelated noise
        end
        
        %Compute and Normalize Reduced Transfer Matrix
        W3 = 2:Nc-1; %channels below i
        reduS = Stemp;
        reduS(2,:,:) = []; reduS(:,2,:) = []; %Remove j from S
        [reduH,reduZ,~] = sfactorization_wilson(reduS,fs,f);
        P2 = [[1,zeros(1,Nc-2)];[-reduZ(W3,1)/reduZ(1,1),eye(Nc-2)]];
        reduHn = zeros(size(reduH));
        for fInd = 1:nFreq
            reduHn(:,:,fInd) = reduH(:,:,fInd)/P2; %Reduced transfer matrix with uncorrelated noise
        end
        
        %Get conditional GC
        for fInd = 1:nFreq
            invG = inv([[reduHn(1,1,fInd),0,reduHn(1,W3,fInd)];[0,1,zeros(1,Nc-2)];...
                [reduHn(W3,1,fInd),zeros(Nc-2,1),reduHn(W3,W3,fInd)]]);
            Qii = invG(1,:)*squeeze(Hn(:,1,fInd)); %Don't need full matrix
            GC(jj,ii,fInd) = real(log(reduZ(1,1)/(Qii*Z(1,1)*conj(Qii))));
        end
    end
end
GC = permute(GC,[3,1,2]);


%---------SUBROUTINES--------------------
%Programs to convert spectral matrices to coherence and factor them into
%transfer matrices and noise covariance. Written by M. Dhamala.
%Minor changes by A. Nakhnikian (commented)

%-------Get auto and cross-spectra--------
function [S,f]= sig2mTspect_nv(X,fs,fRes)
%Usage: [S, f] = sig2mTspect_nv(X,fs,fRes);
%This function computes auto- & cross- spectra by using multitapers
%Inputs: X is multichannel data (a 3D-matrix in the form of time x trial x channel)
%               fs = sampling rate in Hz
%               nv stands for 'not vectorized program'
%               fRes = desired (lower) frequency resolution (e.g. 1), achieved with zero-padding
%              default frequency resolution is fs/datalength
%Outputs: S = 3D matrix: m by m spectral matrix at each frequency point of f
%Note:     One can change nw (half the number of tapers) below and see the effect
%Written by M. Dhamala, USA, August 2006.

[N,Ntr,m] = size(X); % N = timepoints, Ntr = trials, m = channels

if nargin<3||fRes>fs/N %AN added ||
    npad = 0;  fRes = fs/N;
end
if (nargin==3) && (fRes<=fs/N) %AN added &&
    npad = round((fs/fRes-N)/2);  %These many zeros will be padded on each side of the data
end
f = fs*(0:fix((N+2*npad)/2))/(N+2*npad);% upto Nyquist-f

nw = 3; % number of tapers = 2*nw-1 .......good nw are 1.5, 2, 3, 4, 5, 6, or 7...
[tapers,~] = dpss(N+2*npad,nw,2*nw-1);

S = zeros(m,m,N+2*npad); %AN
s = S; %AN
Xft = zeros(N+2*npad,2*nw-1,2*nw-1); %AN
for itrial = 1: Ntr
    for ii = 1: m, Xft(:,:,ii) = mtfft(squeeze(X(:,itrial,ii)),tapers,fs,npad); end
    for ii = 1:m
        for jj = 1:m
            s(ii,jj,:) = squeeze(mean(Xft(:,:,ii).*conj(Xft(:,:,jj)),2));
            %averaging over tapers
        end
    end
    S = S + s;
end
S = S/Ntr; %averaging over trials
S = 2*S(:,:,1:fix(end/2)+1)/fs;%factor 2 for make one-sided spectra
S(:,:,1) = S(:,:,1)/2; %dc-power doesn't double for one-sided case

function xf  = mtfft(data,tapers,~,npad)
%Usage: xf = mtfft(data,tapers,fs,npad);
%Written by M. Dhamala (August 2006)

x0 = zeros(npad,size(data,2));
data = cat(1,x0,data); data= cat(1,data,x0);
data = data(:,ones(1,size(tapers,2)));
data = data.*tapers;xf = fft(data,[],1);

function coh = S2coh(S)
%Input: S auto-& cross pectra in the form: frequency. channel. channel
%Output: coh (Coherence) in the form: frequency. channel. channel
%M. Dhamala, August 2006.

coh = zeros(size(S));
Nc = size(S,2);
for ii = 1: Nc
    for jj = 1: Nc
        coh(:,ii,jj) = real(abs(S(:,ii,jj)).^2./(S(:,ii,ii).*S(:,jj,jj)));
    end
end

%----Factorization-----
function [H, Z, psi] = sfactorization_wilson(S,fs,freq)
%Usage: [H, Z, psi] = sfactorization_wilson(S,fs,freq);
%Inputs: S (1-sided, 3D-spectral matrix in the form of Channel x Channel x frequency)
%            : fs (sampling frequency in Hz)
%            : freq (a vector of frequencies) at which S is given
%Outputs: H (transfer function)
%       : Z (noise covariance)
%       : psi (left spectral factor)
%This function is an implemention of Wilson's algorithm (Eq. 3.1) for spectral matrix factorization
%Ref: G.T. Wilson,"The Factorization of Matricial Spectral Densities," SIAM J. Appl. Math.23,420-426(1972).
%Written by M. Dhamala & G. Rangrajan, USA, Aug 3-4, 2006.
%Email address: mdhamala@gsu.edu

m = size(S,1);N=length(freq)-1; tol = 1E-12; %tol is error-tolerence

%Step 1: Forming 2-sided spectral densities for ifft routine in matlab

f_ind=0;
Sarr = zeros(m,m,2*N); %AN
for f=freq
    f_ind=f_ind+1;
    Sarr(:,:,f_ind)=S(:,:,f_ind);
    if(f_ind>1)
        Sarr(:,:,2*N+2-f_ind)=S(:,:,f_ind).';
    end
end

% Step 2: Computing covariance matricies

gam = zeros(size(Sarr)); %AN
for k1=1:m
    for k2=1:m
        gam(k1,k2,:)=real(ifft(squeeze(Sarr(k1,k2,:)))*fs);
    end
end

%Step 3: Initializing for iterations
gam0 = gam(:,:,1);h = chol(gam0);

%%%%h = rand(m,m); h = triu(h); %arbitrary initial condition

psi = zeros(size(Sarr));
for ind = 1: size(Sarr,3)
    psi(:,:,ind) = h;
end

I = eye(m); % Defining m x m identity matrix

Niterations = 100; % Maximum number of iterations

% Step 4: Iterating to get spectral factors

for iter = 1: Niterations
    g = zeros(size(psi));
    for ind = 1: size(Sarr,3)
        %AN: Replaced "inv" with "\" and "/". Probably doesn't make much
        %difference, but I like to keep Matlab happy
        g(:,:,ind)=psi(:,:,ind)\Sarr(:,:,ind)/psi(:,:,ind)'+I;%Eq 3.1
    end
    gp = PlusOperator(g,m,freq); %gp constitutes of positive and half of zero lags
    psi_old=psi;
    psierr = zeros(1,size(Sarr,3)); %AN
    for k = 1: size(Sarr,3)
        psi(:,:,k) = psi(:,:,k)*gp(:,:,k);
        psierr(k)=norm(psi(:,:,k)-psi_old(:,:,k),1);
    end
    psierrf=mean(psierr);   if(psierrf<tol),break;end % checking convergence
end

%for k = 1: length(freq),
%      Snew(:,:,k) = psi(:,:,k)*psi(:,:,k)'; % Snew: new spectral density
%end

%Step 5: Getting covariance matrix from spectral factors

gamtmp = zeros(size(psi)); %AN
for k1=1:m
    for k2=1:m
        gamtmp(k1,k2,:)=real(ifft(squeeze(psi(k1,k2,:))));
    end
end

% Step 6: Getting noise covariance & transfer function (see Example pp. 424)

A0=gamtmp(:,:,1);

Z = A0*A0.'*fs; %Noise covariance matrix

H = zeros(size(psi)); %AN
for k = 1: length(freq)
    H(:,:,k) = psi(:,:,k)/A0; %Transfer function
end

function gp = PlusOperator(g,m,freq)
%This function is for [ ]+operation:
%   to take the positive lags & half of the zero lag and reconstitute
% M. Dhamala, August 2006

gam = zeros(size(g)); %AN
for k1=1:m
    for k2=1:m
        gam(k1,k2,:)= ifft(squeeze(g(k1,k2,:)));
    end
end

% taking only the positive lags and half of the zero lag

gamp = gam;beta0 = 0.5*gam(:,:,1);
gamp(:,:,1) = triu(beta0);  %this is Stau
gamp(:,:,length(freq)+1:end) = 0;

% reconstituting
gp = zeros(size(gamp));
for k1=1:m
    for k2=1:m
        gp(k1,k2,:)= fft(squeeze(gamp(k1,k2,:)));
    end
end

