function [MIs,fSpace,Hvals,probVects,statsData] = gmwMI(data,fs,lowerFreq,upperFreq,nBins,numTrials,lenTrial,varargin)

% FUNCTION: Computes the modulation index of Tort et al. (PNAS, 2008)
% using the generalized Morse wavelets in lieu of an ad hoc bandpass.  This
% approach enables us to rapidly and efficiently compute the MI across a
% large range of frequencies with principled selection the of time-frequency
% resolution trade-off and. Requires the open source jlab package
% available at http://www.jmlilly.net/jmlsoft.html.

%________________________________________________________
%OUTPUT:

%MIs: A numFrequencies X numFrequencies X numChannels X 
%numChannels array with phase on the y axis and amplitude on the x axis.  
%The i,j,k,t entry is the MI from the phase of frequency i to the amplitude of 
%frequency j when we are gauging the influence of channel k's phase over the
%amplitude of channel t.

%fSpace: The frequency vector in *radians/unit*.  To convert to radians/sec
%take fSpace*fs.  For Hz, use (fSpace*fs)/(2*pi).

%Hvals: An array of entropy values.  Each entry in Hvals is the entropy of
%the phase to amplitude distribution for the corresponding frequencies

%probVects: This array contains the probability densities associated with
%each pair of frequencies, i.e. the probability that amplitude is high or
%low in a particular phase bin

%histBin: A numBins X numFrequencies X numFrequencies X numTrials array
%that contains the phase and amplitude series for individual trials.  Save
%this to be past to the bootstrapping program.

%statsData: A cell array containing amplitude and phase information for
%each individual trial.  This variable must be saved to be passed to the
%program multiTrialBootstrap to set confidence limits.

%________________________________________________________
%INPUT:

%data: This variable can be either a vector or a matrix.  This program
%computes either intra-signal CFC when data is a vector, or inter-signal
%CFC when data is a matrix.

%fs: the sampling rate in Hz

%lowerFreq and upperFreq at the bounds of the frequency range to be
%analyzed.  Enter these values in *Hz*, the program converts them to
%rads/unit before passing them to jlab

%nBins is the number of phase bins.  18 is established in the literature
%(Tort et al., 2008).  360/nBins must be an integer.

%numTrials is the number of trials to analyze (set this to one for
%single-trial)

%lenTrial is the length of each trial in discrete temporal units.  i.e. if
%fs is 500 and each trial is 2 seconds long, then lenTrial = 2*500.


%varargin: An optional set of input variables that allow manual control of
%the wavelet parameters.
%   varargin{1} and varargin{2} are the wavelet parameters gamma and beta.
%   See Olhede and Walden (2002) and Lilly and Olhede (2009).  Both
%   parameters must exceed one.  Gamma is typically 3.  Beta can vary.
%   Higher values of beta mean that the frequency of the windowed carrier
%   wave is higher.  The result is increased spectral resolution at the
%   cost of decreased temporal resolution.  varargin{3} is the wavelet
%   space parameter D (enter "help morsespace" with jlab installed for more
%   information).  Default values for gamma, beta, and D are 3, 5, 4


%________________________________________________________________

%License
%
%  This software is distributed under the "Creative Commons Attribution
%  Noncommercial-Share Alike License"
%
%     Version 3.0, available at
%
%         http://creativecommons.org/licenses/by-nc-sa/3.0/us/
%
%     You are free:
%
%         To Share -- To copy, distribute and transmit the work
%
%         To Remix -- To adapt the work
%
%     Under the following conditions:
%
%         Attribution -- You must attribute the work in the manner specified
%            by the author or licensor (but not in any way that suggests that
%            they endorse you or your use of the work).
%
%         Noncommercial -- You may not use this work for commercial purposes.
%
%         Share Alike -- If you alter, transform, or build upon this work,
%            you may distribute the resulting work only under the same or
%            similar license to this one.
%
%     See the above link for the full text of the license.
%     _______________________________________________________________________
%
%     Disclaimer
%
%     This software is provided 'as-is', without any express or implied
%     warranty. In no event will the author be held liable for any damages
%     arising from the use of this software.
%____________________________________________________
%EXAMPLES:

%Compute MI with default parameters.  Data sampled at 500 Hz, analyzed from
%1 to 200 Hz with 18 bins.  There are 10 trials each 
%lasting 2 secs

%[MI,~,~,~] = gmwMI(X,500,1,200,18,10,1000);

%Compute MI with user defined parameters.  Data sampled at 500 Hz, analyzed from
%1 to 200 Hz with 18 bins.  There are 10 trials each lasting 2 secs

%[MI,~,~,~] =gmwMI(X,500,1,200,18,10,1000,3,8,2);

%_______________________________________________________
%VISUALIZE RESULTS (assumes frequency vector contains 48 entries):

%The following lines return a comodulogram with frequency intervals plotted
%in increments of 4 Hz.  Note that frequencies are logarithmically spaced
%reflecting the reduced frequency resolution of wavelets at smaller scales.

%fSpaceH = (fSpace)*(fs/(2*pi));
%figure, hold on, imagesc(MI'), colorbar
%set(gca,'ytick',1:4:length(fSpace)), set(gca,'xtick',1:4:length(fSpace))
%set(gca,'xticklabel',round(fliplr(fSpaceH([1:4:length(fSpaceH)]))'))
%set(gca,'yticklabel',round((fSpaceH([1:4:length(fSpaceH)]))'))

%______________________________________________________

%References:
%
%   Lilly, J. M. (2016),  jLab: A data analysis package for Matlab, 
%   http://www.jmlilly.net/jmlsoft.html.
%
%   Nakhnikian, A., Ito, S., Dwiel, L.L., Grasse, L.M., Rebec, G.V., 
%   Lauridsen, L.N. and Beggs, J.M., (2016) A novel cross-frequency 
%   coupling detection method using the generalized Morse wavelets. Journal 
%   of neuroscience methods, 269, pp.61-73.

%A. Nakhnikian 2014.  Shinya Ito, Jefferson Davis,
%and Jonathan Lilly suggested valuable improvements to the code.



javaaddpath('C:\Users\Pythia\Documents\GreenLab\code\outside\AlexanderCFC\crossfrequencycoupling.jar');

%suitability checks
if size(data,2)<size(data,1)
    data = data';
end
if length(data)~=numTrials*lenTrial;
    error(['The given values for the trial number and length do not'...
        'match the size of the input vector'])
end

numChannels = size(data,1); %the number of recording sites
%convert to rad/sample
upperFreq = upperFreq*(2*pi/fs);
lowerFreq = lowerFreq*(2*pi/fs);
if nargin == 8
    %Lilly's logarithmic space with default parameters
    gammaa = 3; betaa = 5;
    [fSpace] = morsespace(gammaa,betaa,upperFreq,lowerFreq,4); %frequency vector in rad/samp
else
    gammaa = varargin{1}; betaa = varargin{2}; D = varargin{3};
    %With user defined parameters
    if gammaa~=3
        warning(['The Morse wavelet has time-frequency resolution near'...
            ' the Morlet wavelet when gamma = 3.  Setting gamma to'...
            'another value is not recommended.  Enter command ''warning off'''...
            'to disable this message'])
    end
    [fSpace] = morsespace(gammaa,betaa,upperFreq,lowerFreq,D); %frequency vector in rad/samp
end

rect = fSpace/morsefreq(3,6); %used to rectify wavelet bias (Liu et al, 2007)
MIs = zeros(length(fSpace),length(fSpace),numChannels,numChannels);
Hmax = log(nBins); %entropy of uniform dist.
[psi] = morsewave(lenTrial,1,gammaa,betaa,fSpace,'energy');
channelCount = 0; %used to track time in next loop
statsData = cell(2,numChannels);
for channelInd1 = 1:numChannels
    for channelInd2 = 1:numChannels
        tic
        fprintf(['Computing CFC from Channel ',num2str(channelInd1),...
            ' to Channel ',num2str(channelInd2),'\n']);
        %extract amplitude and phase information
        x = reshape(data(channelInd1,:),lenTrial,numTrials);
        y = reshape(data(channelInd2,:),lenTrial,numTrials);
        [amplitudeArray,phaseArray] = waveletBandPass(x,y,psi,rect);
        statsData{1,channelInd1} = amplitudeArray; statsData{2,channelInd2} = phaseArray;
       
        %Construct the joint time series (phase and energy)
        [histBin] = histBuilder(phaseArray,amplitudeArray,nBins);


        
            trialAverages = mean(histBin,4);
            MIsTemp = zeros(size(trialAverages,2),size(trialAverages,3));
            Hvals = MIsTemp;
            probVects = zeros(size(trialAverages));
            trialAverages(trialAverages == 0) = eps; %log of true zero is undefined
            for phaseIndex = 1:size(trialAverages,2)
                for energyIndex = 1:size(trialAverages,3)
                    histVect = squeeze(trialAverages(:,phaseIndex,energyIndex));
                    probVect = histVect./sum(histVect); %histogram to PDF
                    Hobtained = -(probVect'*log(probVect)); %entropy via inner product
                    modulationIndex = (Hmax-Hobtained)/Hmax; %Get the MI
                    MIsTemp(phaseIndex,energyIndex) = modulationIndex;
                    Hvals(phaseIndex,energyIndex) = Hobtained;
                    probVects(:,phaseIndex,energyIndex) = probVect;
                end
            end
       

        MIs(:,:,channelInd1,channelInd2) = MIsTemp;
        channelCount = channelCount+1;
        t = toc;
        fprintf(['Approximate time remaining is ',...
            num2str(((numChannels^2*t)-(t*channelCount))/60^2),'hours \n'])
    end
end

end
function [histBin] = histBuilder(x,y,nBins)
if mod(360,nBins)~=0
    error('The number of bins must be an integer factor of 360')
end
phaseInt = 360/nBins;
histBin = zeros(nBins,size(x,2),size(x,2),size(x,3));
for trialIndex = 1:size(x,3)
    phaseSeries = x(:,:,trialIndex);
    phaseSeries = (180*phaseSeries)./pi; %convert the phase series to degrees
    phaseSeries = phaseSeries+180; %set the bound to [0,360]
    amplitude = y(:,:,trialIndex);
    phaseSeriesFreqIndex = ceil(phaseSeries / phaseInt);
    phaseSeriesFreqIndex(phaseSeriesFreqIndex == 0) = nBins; % because 0 degree is the same as 360 degrees, avoids possible 
    %Java error if the value is "exactly" zero
    %Java script called here
    histBin(:,:,:,trialIndex) = crossfrequencycoupling.CFCalc.fastbuilderjava(phaseSeriesFreqIndex, amplitude, nBins);
end
end

function [amplitude,phaseSeries] = waveletBandPass(x,y,psi,rect)
X = wavetrans(x,psi);
Y = wavetrans(y,psi);
if numel(size(X)) == 4||numel(size(x))<3&&numel(size(psi))>2
    %make each eigenspectrum a "trial"
    X = reshape(X,size(X,1),size(X,2),size(X,3)*size(X,4));
    Y = reshape(Y,size(Y,1), size(Y,2), size(Y,3)*size(Y,4));
    %Y = repmat(Y(:,:,:,1),1,1,size(psi,3));
end
amplitude = abs(X);
phaseSeries = angle(Y);
end
