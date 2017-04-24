%% Fourier Notes
% Before we dive into Fourier analysis we will briefly go over some of the
% basic terminology and features of the kind of data we are interested in.
% In doing so, we will hopefully have a better understanding of our data
% and the general motivation behind Fourier analysis.
 
% From the earliest recordings of the brain it has been clear that
% oscillations are a key feature of brain activity. Briefly, an oscillation
% can be understood as a periodic, or rhythmic, pattern; when thought of in
% the context of time, this would mean a pattern that repeats every x units
% of time. A useful example of a periodic function of time are those from
% trigonometry: sin(t) and cos(t), where t is time.
 
% As we've seen before, to plot these functions in Matlab
% requires two vectors corresponding to the x and y axes, or time and
% amplitude respectively. 
 
% fs is a conventional notation indicating sampling frequency, i.e. number
% of cycles per second (Hz)
fs = 100; 
% Set up the start and end times for our time vector, here in seconds
tStart = 0; tEnd = 1;
% Use these start and end points alongside the sampling frequency to
% construct our time vector; note that this line is equivalent to asking
% Matlab to make a long series of numbers that starts at 0, ends at 1, and
% takes steps of 0.01 (1/100)
tvec = tStart:1/fs:tEnd;
% Next we will set the frequency of sine wave that we wish to plot; also in
% Hz (cycles/second), N.B. the inverse of Hz is also called the period or
% the number of seconds required for a single oscillation, in this case
% 0.25 seconds.
f = 4; 
% The next line of code takes a little bit to unpack. First note that we
% are using the Matlab function sin() which if you look at its
% documentation expects a certain argument, namely a vector of angles in
% radians. The angles can be thought of as the angle between a point moving
% on a unit circle and the x axis. (INSERT ANIMATION HERE) All we need to
% do is multiply each of time values by our frequency, 4 Hz, and then
% convert into radians using 2*pi. This will give back a vector
% representing a periodic function going around the unit circle twice every
% second.
y = sin(2*pi*f*tvec);
% We can then visualize this function by first setting up an empty figure
% and then feeding our time vector and new sine vector in place for x and y
% in the plot() function.
figure
plot(tvec,y)
% A good plot always has labels which we can easily add in Matlab. To do
% this we take advantage of some more simple functions and a new notation
% called strings. Strings are nothing more than little packets of letters
% that Matlab recognizes through single apostrophes around the text.
xlabel('Time (s)')
ylabel('Amplitude')
title('4 Hz Oscillation')
 
% As our y label suggests, there are a few more ways to describe a
% oscillation. These are amplitude and phase, both of which are easy to
% manipulate in Matlab. Amplitude is a measure describing the distance
% along y axis at a given point and in the case of the wave we just plotted
% the math works out such that the amplitude is actually a unitless ratio
% that goes from 1 to -1; in case of natural signals the amplitude depends
% on the system in question, for many neural signals this is measured in
% millivolts. N.B. amplitude is a relative measure by which I mean that it
% is the distance along the y axis from some reference point; people use
% different reference points although a common on is the x axis itself
% (where y = 0) and in this context the amplitude at a given point is
% simply its y coordinate.
 
% Amplitude of a sine wave is manipulated through simple multiplication of
% the entire function. Remember that plotting in Matlab is done using sets
% of Cartesian coordinates (x,y). So if we wish to double the distance from
% our y axis at every instance of x, we simply multiply all of our y values
% by 2.
a = 2;
y2 = sin(2*pi*f*tvec)*a;
figure
plot(tvec,y2,'r')
% We will also plot our first signal for comparison by using the command
% 'hold on' which tells Matlab to overlay the next plot rather than getting
% rid of the previous plot.
hold on
plot(tvec,y,'b')
title('4 Hz Oscillation: Double Amplitude')
ylabel('Amplitude')
xlabel('Time (s)')
% Note that Matlab automatically rescales the y axis to fit the entire
% signal. I've also told Matlab what color to plot the lines to match the
% first plot by adding a string as the third argument; 'r' or 'b' to
% indicate a red or blue line respectively.
 
% The second feature of periodic signals, phase, is a measure of how far in
% a single oscillation the signal currently is at. Similarly to amplitude,
% the measure of phase is relative to some reference and thus arbitrary.
% Some common conventions are the distance in a complete oscillation from
% zero (it is important to remember 'in a complete oscillation' since a
% sine wave for example will cross the zero line once before a complete
% oscillation), or the distance from peak to peak or trough to trough.
% Typically phase has units of either radians or degrees. We can think
% about manipulating phase as simply shift our signal along the x axis,
% implying that all we need to do is add a number equivalent to how much we
% wish to shift our signal to each value of our input to sin().
a = 1;
phi = pi/2; 
y3 = sin(2*pi*f*tvec + pi/2)*a;
figure
plot(tvec,y3,'g')
hold on
plot(tvec,y,'b')
% Note that to use Greek letter in plots, use a backslash before the letter
title('4 Hz Oscillation: \pi/2 Phase Shift')
xlabel('Time (s)')
ylabel('Amplitude')
% You might be able to imagine that by changing these two features
% alongside the frequency of oscillation we can make any pure oscillator we
% desire.  But what about more complex waveforms?
 
% All we have to do to generate a more complex waveform is add up pure
% oscillators point by point. For example, we can add a 1 Hz wave to our 4
% Hz wave. 
f = 1;
y4 = sin(2*pi*f*tvec);
yCombine = y + y4;
figure
% To draw multiple plots in the same figure use subplot to tell Matlab
% which panel you are about to plot in using the following notation:
% subplot(m,n,p) where m = number of rows; n = number of columns; and p =
% which position to plot.
 
subplot(3,1,1);
plot(tvec,y)
% Using the command ylim() with a upper and lower bound pair in as the
% argument to manually set the limits to the y axis for visualization
% purposes
ylim([-2 2])
title('4 Hz Oscillation')
subplot(3,1,2);
plot(tvec,y4)
ylim([-2 2])
title('1 Hz Oscillation')
ylabel('Amplitude')
subplot(3,1,3);
plot(tvec,yCombine);
xlabel('Time (s)')
title('Combination')
 
% By manipulating phase, amplitude, and frequency of many oscillators and
% adding them up it looks like we can create complicated waveforms. This
% simple idea is actually the groundbreaking theory behind Fourier
% analysis. To play with this idea I recommend either writing your own code
% to make the process faster, or using the script I wrote, sigGen, which
% will be used throughout the rest of this write-up to generate signals.
% For example the following line will generate a combination of the
% following waves using a sampling rate of 1000 over 2 seconds: (1) 2 Hz,
% 1.5 amplitude, 0 phase offset; (2) 7 Hz, 1 amplitude, 90 degree phase
% offset; and (3) 9 Hz, 0.5 amplitude, 45 degree phase offset.
[x,t] = sigGen(1000,2,[2,7,9],[1.5,1,0.5],[0,90,45]);
figure
plot(t,x)
 
% However, before we can jump into the fun of analyzing our data, we need
% to understand the data collection process.
%% Sampling
% Signals in the natural world are analog rather than digital, i.e. analog
% signals are continuous with an infinite set of possible values, whereas a
% digital signal is discrete and has a finite or bounded set of possible
% values. Another way to think of this is that an analog signal has an
% infinite amount of information in both time (x axis) and amplitude (y
% axis), i.e. no matter how small of a gap there is between x1 and x2 (time
% 1 and time 2), there will be an infinite number of steps between them and
% the same is true for y1 and y2. N.B.: a digital signal can have a very
% small, short, time step and very large set of possible values.
 
% Visual example of analog vs. digital signal.
fs = 100; % Sampling frequency
t = 0:1/fs:1; % Time vector with steps of 1/fs
f = 5; % Frequency of signal (Hz)
x = sin(2*pi*f*t); % Signal
% N.B.: in the future this process will be simplified by using the function
% sigGen to generate simple signals
 
% Plot
figure
subplot(2,1,1)
title('Continuous/Analog')
hold on 
plot(t,x) % Plot continuous (analog) signal
subplot(2,1,2)
stairs(t,x)
title('Discrete/Digital')
xlabel('Time (s)')
 
% However, collecting and analyzing an infinite amount of data is
% impossible. As such we convert our natural, analog, signal into a digital
% one. The first step is to reduce the amount of information in the
% temporal dimension (x axis), i.e. sampling. Sampling is the process of
% regularly pulling individual values, samples, from the analog signal. The
% rate at which we pull these values is the sampling rate; conventionally:
% the number of samples per second.
 
% The goal of sampling is to model the original signal with relatively high
% fidelity, i.e. we want to use as few samples as possible while still
% being able to faithfully recreate the original signal. As such we will be
% limited by a few factors: the stability of the original signal (over what
% timescale does the signal make large deviations). 
 
% Since samples are being taken from a finite time interval, one can think
% of the process as holding, or 'freezing' the signal every T seconds and
% obtaining a single value (in the case of LFPs, a voltage).  This can be
% visualized as multiplying a series of sampling impulses. 
 
% Generate signal with 3 and 5 Hz components - analog (input) data
[x,t] = sigGen(100,1,[4,5],1);
% Set up sampling impulse vector with 10 Hz fs
fs = 10;
pulseT = 0:1/fs:t(end);
pulseVect = ones(length(pulseT));
% Represent multiplication of impulses with signal
sampVect = x(1:floor(length(t)/fs):end);
sampT = t(1:floor(length(t)/fs):end);
% Plot
figure
hold on
subplot(4,1,1)
plot(t,x)
title('Input Signal')
subplot(4,1,2)
stem(pulseT,pulseVect)
title(['Sampling Impulses: ',num2str(fs),' Hz'])
ylim([0 1.5])
subplot(4,1,3)
stem(sampT,sampVect);
title('Sampled Signal')
subplot(4,1,4)
stairs(sampT,sampVect);
title('Reconstructed Signal')
xlabel('Time (s)')
% Note the difference between the reconstructed signal and the input
% signal. Change the sampling rate, what happens as it approaches the
% length of the original signal?
%% Quantization
% However, there is a problem with the process described above. We've
% already mentioned that there is infinite precision along the y axis as
% well as x axis, so what value exactly are we pulling out of the signal
% when we sample? Without going into too much detail, as we take samples
% from our signal we also apply a kind of threshold or rounding to the
% exact value. How much rounding occurs depends on the bit resolution of
% the converter used; e.g. Plexon tends to have 12-bit resolution between
% 10 and -10 V. This means that it can resolve 4.9 mV: (20 volt range)/2^12
%% Aliasing and the Nyquist Frequency
% Going back to sampling; how do we choose the right sampling rate?
% Ultimately this is measured by how accurately you are able to represent
% the original signal, i.e we want to minimize the average error between
% the reconstructed signal and original signal. A good starting place is to
% determine the minimum sampling rate.
 
% To do this we need to know something about the question we are trying to
% answer. Specifically, we need to know what the highest frequency is that
% we wish to be able to reconstruct. The reason for this is something
% called aliasing - which is the misidentification of a wave, usually by a
% slower wave, because our sampling rate is too low. Let's go through an
% example:
figure
% Create a 10 Hz signal
[x,t] = sigGen(1000,1,10,1);
plot(t,x)
title('10 Hz Wave')
% Lets sample this signal at 8 Hz
fs = 8;
sampT = 0:1/fs:1;
sampVect = interp1(t,x,sampT,'nearest');
title('10 Hz Wave Sampled at 8 Hz')
hold on 
stem(sampT,sampVect)
% Note that if we sample at 8 Hz, it looks like there is more than one
% signal that the samples could come from, i.e. . In other words, when we
% go to reconstruct the signal there would be more than one option to
% choose from and the simpler might be a slower wave (~2 Hz judging by the
% number of peaks in one sec) than the original 10 Hz wave. We can double
% check this by plotting a 2 Hz wave.
[x2,t2] = sigGen(1000,1,2,1);
plot(t2,x2,'--r')
% But maybe this is just because we sampled at a lower rate than the
% frequency of the wave we are interested; let's next try 15 Hz and see if
% that works.
figure
plot(t,x)
title('10 Hz Wave Sampled at 15 Hz')
hold on
fs = 15;
sampT2 = 0:1/fs:1;
sampVect2 = interp1(t,x,sampT2,'nearest');
stem(sampT2,sampVect2)
% At first glance this looks a little better, but still slower than the 10
% Hz wave we know is there. In fact it looks like 5 Hz.
[x3,t3] = sigGen(1000,1,-5,1);
plot(t3,x3,'--r')
% If you were to slowly increase the sampling rate you would notice that
% once fs approaches 20 Hz, this problem would go away. This is what is
% called the Nyquist frequency - the highest frequency able to reproduced
% without error is half of the sampling rate. Another way to phrase this is
% the Nyquist rate - the minimum sampling rate usable to accurately
% reproduce the highest frequency in the signal.
% Let's look at an example:
figure
plot(t,x)
title('10 Hz Wave Sampled at 25 Hz')
hold on
fs = 25;
sampT3 = 0:1/fs:1;
sampVect3 = interp1(t,x,sampT3,'nearest');
stem(sampT3,sampVect3)
% However, it is not the best idea to choose the Nyquist frequency exactly.
% To visualize why, try sampling the signal at exactly 10 Hz. Instead it is
% better practice to err on the side of a higher sampling rate, say 3-4x
% the highest frequency you are interested in.
 
% A few last notes on aliasing. It is typical to apply what is called an
% anti-aliasing filter before sampling the data. This filter - a low pass
% filter - is designed to let low frequencies through and filter out high
% frequencies in a systematic way to limit the effects of aliasing from
% high frequency (> Nyquist frequency) components. The value in doing this
% is to control, somewhat, what might be happening between samples. 
%%
% Generate a signal with 3 and 10 Hz components with a sampling rate of
% 1200 over one second; in this example 
[x,t] = sigGen(1200,1,[3,10],[1 0.5]);
figure
subplot(2,1,1)
plot(t,x)
% If we sample this signal at 12 Hz we will get aliasing due to twelve
% being less than our Nyquist frequency (2*10 Hz = 20 Hz). N.B.: sampling
% our signal (1200 Hz sampling rate) at 12 Hz is equivalent to taking every
% 100th sample.
fs = 12;
sampT = 0:1/fs:1;
sampVect = interp1(t,x,sampT,'nearest');
subplot(2,1,1)
hold on
stem(sampT,sampVect,'.','MarkerSize',10)
title('Signal without Anti-Alias Filter')
% Note that we do not capture our pure 3 Hz component. Next use decimate()
% which has a built in low-pass filter.
% Here we will take every 100th sample of both our time vector and signal
sampT2 = decimate(t,100);
sampVect2 = decimate(x,100);
subplot(2,1,2)
% Plot pure 3 Hz sine wave for visualization
[x2,t2] = sigGen(1200,1,3,1);
plot(t2,x2)
title('Signal with Anti-Alias Filter')
hold on
stem(sampT2,sampVect2,'.','MarkerSize',10)
ylabel('Time (s)')
% Note that this is better than without the built in anti-aliasing filter
% of decimate, although you may notice some edge effects (first and last
% samples) which we will discuss further when we dive into filters. The
% take home message here is that there exists a limit in one's ability to
% parse apart the components of a signal - defined by the trade-off between
% the sampling rate of data acquisition and highest frequency of interest -
% which is defined as the Nyquist limit. To account for this, we typically
% apply an anti-alias filter to cut out frequencies above the Nyquist
% frequency and then oversample the analog signal.
 
% A practical example from actual data acquisition: Plexon's DigiAmp
% (analog-digital converter and preamplifier) has a sampling rate of 40 kHz
% (40,000 samples a second; 25 microsecond resolution) which in practice
% means that the highest frequency resolvable is 20 kHz (40kHz/2). Usually
% an anti-aliasing filter is designed to remove components > Nyquist
% frequency (20 kHz in this case). However, since neural signals aren't
% thought to have naturally occurring components at this frequency, Plexon
% instead uses a 8 kHz low-pass filter. By doing this we also get the
% benefits of oversampling our signal - at the lowest we will be sampling a
% 8 kHz wave at 40 kHz, 5x oversampling. For the cases of physiological
% signals, usually capped around 250 Hz, we are oversampling by 1600x. The
% benefits of oversampling are increased resolution of our components of
% interest, reduced impact of noise (higher signal to noise ratio), and
% more relaxed anti-aliasing filter construction (larger possible
% transition band).
%% Fourier Analysis
% Now that we have a better understanding of our data in the context of
% oscillations and can appreciate how we obtained the data we can begin to
% delve into the magic of Fourier analysis.
 
% As mentioned in the first section, Fourier's simple, but groundbreaking,
% theory was that using a infinite set of sine wave which are integer
% multiples of one another any periodic waveform, complex or simple, could
% be represented as an infinite sum of integer multiple sine waves. A
% corollary to this is that any periodic waveform can be decomposed into
% its component parts - a Fourier series. Although cool, in the realm of
% neuroscience it would be nice to not need a infinite amount of anything
% (especially since we've already discovered that we are dealing with
% discrete and not continuous data) nor be constrained to integer multiples
% (nature is rarely so neatly ordered, although music theory offers a
% interesting exception). Luckily for us, mathematicians have found a way
% around these obstacles: the discrete Fourier transform (DFT).
 
% The DFT not only discretizes along time, but also frequency. What this
% means for us is that we do not need a continuous signal and we are no
% longer constrained by only looking at integer harmonics. Instead, we can
% look at as many discrete frequencies that we want. The only limit is that
% depending on our desired resolution in the frequency domain, we will
% need higher resolution in the time domain (i.e. more samples). This will
% become clearer through examples.
 
% Let's start with a pure sine wave for which we know a priori the spectral
% components for.
% We are using 1024 as our sampling rate because the fast Fourier
% transform (the method utilized in the DFT to do the heavy lifting) works
% fastest when the number of samples is a factor of 2.
[x,t] = sigGen(1024,1,4,1,0);
xF = fft(x);
% The output of the FFT is a series of complex numbers which gives us both
% the magnitudes and phase components of our signal as a complex number. To
% extract the magnitudes we just take the absolute value of our transformed
% signal (X).
magX = abs(xF);
% To extract the phase information we use the function angle() which
% converts a complex number into radians.
phiX = angle(xF);
figure
stem(magX)
hold on
plot(magX)
% However, plotting our magnitude component as is does not easy to
% interpret. It is unclear what the axes should be (where is frequency? the
% frequency of our signal is 4 Hz and not 0 or 1000...) and why there are
% two peaks when we looking at a pure 4 Hz sine wave. The reason why the
% frequency axis is not immediately apparent is that fft doesn't deal with
% frequencies per se, but rather samples. Further, the fft automatically
% cuts off at the Nyquist frequency (fs/2). 
nPoints = length(x);
F = (-nPoints/2:nPoints/2-1);
% The reason for the double peaks is that the Fourier transform itself is
% complex and so when we give it a real signal it still gives back a
% complex output where the negative frequencies correspond to the imaginary
% component. Since we are only dealing with real signals, we can ignore the
% negative component. First we will use fftshift() to swap the real and
% imaginary components around so that the imaginary portion comes first and
% the real portion second to match our frequency axis.
%%
xFshift = fftshift(xF);
magXshift = abs(xFshift);
figure
stem(F,magXshift)
hold on
plot(F,magXshift)
xlabel('Frequency (Hz)')
% Usually we won't be using fft() manually to compute the Fourier transform
% so we don't have to worry about clipping out the negative side of this
% plot. We've plotted the Fourier coefficients in terms of their magnitude
% by using abs() and removing the imaginary/negative components. However,
% you may be aware that in signal processing it is conventional to plot in
% power (e.g. Power Spectral Densities). Power is merely the square of the
% amplitude.
powX = magXshift.^2;
% When we want to compare two signals, just squaring our outputs will make
% it difficult to interpret; two components with slightly different
% amplitudes can have very different powers, e.g. an amplitude of 4 is
% twice as large as an amplitude of 2, but the power is 4x larger (4^2 = 16
% vs. 2^2 = 4), the difference in power will always be the square of the
% difference in amplitude. All of this is to say that we need a better way
% to plot our power so that it is visually easy to interpret. The
% convention here is to use decibels which can be obtained using a base 10
% logarithm and scaling up it up 10 times since most of our differences
% will be relatively small.
dbX = 10*log10(powX);
figure
plot(F,dbX)
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
ylim([-50 75])
xlim([0 500])
% Notice that once we've converted our magnitude into decibels (dB) we now
% have negative values. Don't worry, this is because dBs are not an
% absolute measure, rather they are relative to a baseline (0 dB) which
% corresponds to the smallest amount of pressure (quietest sound) a human
% ear can detect. So 10 dB would be 10 times louder than 0 dB and -20 dB is
% 100 times weaker than 0 dB.
 
% If you zoom in on our plot, you can see that the estimate is not exactly
% at 4 Hz 


