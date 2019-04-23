function [sigMIs,h,pvals,critP,adjP] = gmwCFCpermFDR(MI,statsData,nBins,iterates,q)

%FUNCTION: Computes p vlaues corresponding to each MI and tests for statistical
%significance using FDR control for multiple comparisons. q is the expected
%value of true discovery/false discovery.  See the comments in fdr_bh for
%further details.

%INPUT:

%MI: Modulation indices returned by gmwCFC

%statsData: Data randomized to generate null distributions of CFC. This
%variable is returned by gmwCFC

%nBins: Number of phase bins. This should be identical to the number of
%bins used to generate the raw data

%iterates: The number of iterations used to generate randomized CFC
%distributions. 200 generally returns reasonable results. If the user
%suspects there is an unacceptable number of false positives it is
%advisable to rerun the analysis with a larger number of iterates

%q: The proportion of incorrectly rejected true null hypotheses. i.e. if 1 =
%0.5 then at most 5% of true null hypotheses are rejected. The actual
%type I error rate might be lower, q is an upper bound

%OUTPUT:

%sigMIs: A modulation index array containing only MIs above chance

%h: A binary array containing ones for significant MIs and zeros elsewhere

%pvals: Empirical p-values for each MI determined from the emipirical
%cumulative distribution of each set of randomized MI values

%critP: The p value corresponding to MIs above chance at the desired false
%discovery rate

%adjP: Effective p values after FDR correction

%DEPENDENCIES: 

%Matlab statistics (or statistics and machine learning) toolbox

%fdr_bh availble on Mathworks exchange and included in the GMW CFC package
%Author: David M. Groppe

%License:

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

% References:
%
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%   rate: A practical and powerful approach to multiple testing. Journal
%   of the Royal Statistical Society, Series B (Methodological). 57(1),
%   289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%   rate in multiple testing under dependency. The Annals of Statistics.
%   29(4), 1165-1188.
%
%   Nakhnikian, A., Ito, S., Dwiel, L.L., Grasse, L.M., Rebec, G.V., 
%   Lauridsen, L.N. and Beggs, J.M., (2016) A novel cross-frequency 
%   coupling detection method using the generalized Morse wavelets. Journal 
%   of neuroscience methods, 269, pp.61-73.

%A. Nakhnikian, 2016
javaaddpath('C:\Users\Pythia\Documents\GreenLab\code\outside\AlexanderCFC\crossfrequencycoupling.jar');
randMIs = zeros([size(MI),iterates]);
if numel(size(MI))<3
    numChan =1;
else
    numChan = size(MI,3);
end
Hmax = log(nBins);
for randInd = 1:iterates
    disp(['Computing Randomization ' num2str(randInd) ' of ' num2str(iterates)])
    for chanInd1 = 1:numChan
        for chanInd2 = 1:numChan
            %Construct the joint time series (phase and energy)
            amplitudeArray = statsData{1,chanInd1};
            phaseArray = statsData{2,chanInd2};
            % Use circshift for phaseArray and keep amplitudeArray the same
            % - Edits LLD 2018
            phaseArrayRand = zeros(size(phaseArray));
            for ii = 1:size(phaseArray,3)
               phaseArrayRand(:,:,ii) = circshift(phaseArray(:,:,ii),randi(size(phaseArray,1),1)); 
            end
            amplitudeArrayRand = amplitudeArray;
%             randInds1 = randperm(size(phaseArray,3),size(phaseArray,3));
%             randInds2 = randperm(size(phaseArray,3),size(phaseArray,3));
%             phaseArrayRand = phaseArray(:,:,randInds1);
%             amplitudeArrayRand = amplitudeArray(:,:,randInds2);
            [histBin] = histBuilder(phaseArrayRand,amplitudeArrayRand,nBins);
            %[histBinTest] = histBuilder(phaseArray,amplitudeArray,nBins);
            
            
            trialAverages = mean(histBin,4);
            modIndsTemp = zeros(size(trialAverages,2),size(trialAverages,3));
            Hvals = modIndsTemp;
            probVects = zeros(size(trialAverages));
            trialAverages(trialAverages == 0) = eps; %log of true zero is undefined
            for phaseIndex = 1:size(trialAverages,2)
                for energyIndex = 1:size(trialAverages,3)
                    histVect = squeeze(trialAverages(:,phaseIndex,energyIndex));
                    probVect = histVect./sum(histVect); %histogram to PDF
                    Hobtained = -(probVect'*log(probVect)); %entropy via inner product
                    modulationIndex = (Hmax-Hobtained)/Hmax; %Get the MI
                    modIndsTemp(phaseIndex,energyIndex) = modulationIndex;
                end
            end
            if numChan == 1
                randMIs(:,:,randInd) = modIndsTemp;
            else
                randMIs(:,:,chanInd1,chanInd2,randInd) = modIndsTemp;
            end
        end
    end
    
end

disp('Computing p values and applying FDR correction')
if numChan == 1
    pvals = zeros(size(randMIs,1),size(randMIs,2));
    for ind = 1:size(randMIs,1)
        for ind2 = 1:size(randMIs,2)
            [f,y] = ecdf(squeeze(randMIs(ind,ind2,:))');
            critVal = f(find(y<=MI(ind,ind2),1,'last'));
            if isempty(critVal)
                pvals(ind,ind2) = 0;
            else
                pvals(ind,ind2) = 1-critVal;
            end
        end
    end
    [h,critP,adjP] = fdr_bh(pvals,q,'dep',0);
    sigMIs = MI; sigMIs(~h) = 0;
else
    sigMIs = zeros(size(MI));
    h = sigMIs;
    adjP = sigMIs;
    pvals = sigMIs;
    critP = zeros(numChan,numChan);
    for chanInd = 1:numChan
        for chanInd2 = 1:numChan
            randMIsTemp = squeeze(randMIs(:,:,chanInd,chanInd2,:));
            pvals1 = zeros([size(randMIs,1),size(randMIs,2)]);
            for ind = 1:size(randMIsTemp,1)
                for ind2 = 1:size(randMIsTemp,2)
                    [f,y] = ecdf(squeeze(randMIsTemp(ind,ind2,:))');
                    critVal = f(find(y<=(MI(ind,ind2,chanInd1,chanInd2)),1,'last'));
                    if isempty(critVal)
                        pvals1(ind,ind2) = 0;
                    else
                        pvals1(ind,ind2) = 1-critVal;
                    end
                end
            end
            [h1,critP1,adjP1] = fdr_bh(pvals1,q,'dep',0);
            sigMIstemp = squeeze(MI(:,:,chanInd,chanInd2));
            sigMIstemp(~h1) = 0;
            pvals(:,:,chanInd,chanInd2) = pvals1;
            h(:,:,chanInd,chanInd2) = h1;
            critP(chanInd,chanInd2) = critP1;
            adjP(:,:,chanInd,chanInd2) = adjP1;
            sigMIs(:,:,chanInd,chanInd2) = sigMIstemp;
        end
        
    end
end