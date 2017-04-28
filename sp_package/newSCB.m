x = LFPTs.data(:,1:5000);
[nChan,lData] = size(x);
k = 10; gam = 3; bet = 6; d = 10;
high = 100*(2*pi/adfreq);
low = 1*(2*pi/adfreq);
%% Continuous
if strcmpi(aType,'cont')
   % Either use wavelets with frequency dependent windows or STFT with one
   % window size
   if strcmpi(fType,'wave')
       % Determine energy concentration of wavelet with parameters beta,
       % gamma, and D
       % r from Brittain et al. 2007
       r = (2*bet+1)/gam;
       % g and c, reorganized from Brittain et al. 2007 [10]
       g = gamma(r+1-1/gam)*gamma(r+1/gam)*gam^-1*gamma(r)^-2;
       c = (d+g)/g;
       % Find top n eigenvalues corresponding to their energy ratio
       % relative to original signal; only use these n ks for further
       % analysis and smoothing
       % Start with 0th eigenspectra
       k = 0:k-1;
       % Extract eigenvalues using incomplete beta function
       I = betainc((c-1)/(c+1),k+1,r-1);
       % Transform into energy ratio; equivalent to energy of
       % eigenspectrum/energy of original signal
       eR = I.^2;
       % Find top n eigenvalues which exceed energy cutoff
       eThresh = 0.95;
       eInds = logicFind(eThresh,round(eR,2),'>=');
       if isempty(eInds)
          error(['Warning: None of the first ',num2str(k(end)+1),' '...
              'eigenspectra using parameters given (beta = ',num2str(bet),...
              ', gamma = ',num2str(gam),', D = ',num2str(d),') have energy'...
              ' ratios above the cutoff (',num2str(eThresh),'). Either'...
              ' change parameters or change energy threshold'])
       end
       % Prepare smoothing operator - vector of weights
       eW = eR(eInds)./sum(eR(eInds));
       % Replace k, with number of eigenvalues above threshold
       k = numel(eInds);
       % Construct log spaced frequencies for Morse wavelets
       fLog = morsespace(gam,bet,high,low,d);
       % Construct wavelet bank
       psi = morsewave(lData,k,gam,bet,fLog,'energy');
       % Conduct wavelet transform - time X frequency X channel X
       % eigenspectra
%        w = wavetrans(x(1,:)',psi);
       w = wavetrans(x',psi);
       % Calculate auto- and cross-spectra
       for ii = 1:size(w,3)
           auto(:,:,ii,:) = abs(w(:,:,ii,:)).^2;
           for k = 1:size(w,3)
               cross{ii}(:,:,k,:) = w(:,:,ii,:).*conj(w(:,:,k,:)); 
           end
       end
       %% Normalize
       % Reshape weight vector - put each weight on its own page
       eWr = reshape(eW,1,1,1,numel(eW));
       % Transform weight vector into appropriately sized matrix
       weight = repmat(eWr,size(w,1),size(w,2),size(w,3));
       sAuto = sum(weight.*auto,4);
       for ii = 1:size(cross,2)
           sCross{ii} = sum(weight.*cross{ii},4);
       end
       %% COI
       % Get standard deviation in time domain 
       sigmaT = sqrt(1/(gammaHat((2*bet+1)/gam))*((bet^2)*gammaHat((2*bet-1)/gam)+gam^2*gammaHat((2*bet+2*gam-1)/gam)-2*bet*gam*gammaHat((2*bet+gam-1)/gam)));
       % Get peak frequency
       [~,eF,~,~] = morsefreq(gam,bet);
       % Calculate Fourier Factor - converts from scale to period (1/f)
       fF = (2*pi)/eF;
       % Calculate COI scalar, then propogate across time vector
       coiScalar = fF/sigmaT;
       dt = 1/adfreq;
       coi = coiScalar*dt*[1E-5,1:((lData+1)/2-1),fliplr((1:(lData/2-1))),1E-5];
       %% Plot Power
       wPSD = sAuto;
       % Plot power
       % Convert fLog back to Hz.
       f = (fLog)*(adfreq/(2*pi));
       % Set up time vector
       t = LFPTs.tvec(1:5000);
       % Remove last index
       figure
       for ii = 1:nChan
           subplot(2,2,ii)
           imagesc(t,f,wPSD(:,:,ii)')
           colormap('viridis')
           set(gca,'YDir','normal')
           hold on
           % Convert COI from period to frequency and shade
           area(t,1./coi,'FaceColor','w','EdgeColor','w','FaceAlpha',0.5)
           title(LFPTs.label{ii})
           xlabel('Time (sec)')
           ylabel('Frequency (Hz)')
       end
       %% Calculate Coherence
       for ii = 1:size(sCross,2)
           
       end
       test = abs(sCross{1,1}(:,:,2)).^2./(sAuto(:,:,1).*sAuto(:,:,2));
   elseif strcmpi(fType,'fft')
       % Power
       [s,f,t,ps] = spectrogram(data,);
       % Coherence
       
   end
end
%% Trialized
if strcmpi(aType,'trial')
    if strcmpi(fType,'wave')
       
   elseif strcmpi(fType,'fft')
       
   end
end