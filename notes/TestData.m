%% Create signals for spectcompbase
adfreq = 2000; % Sampling rate (Hz)
t0 = 0; t1 = 25; % Start and end times, seconds
tvec = t0:1./adfreq:t1; % Time axis

mag = [.1;.15;.2;.25;.3;.35;.4;.45;.5;.55;.6;.61];
f = [80;80;50;50];
% mag = [.1 .2]; 
% f = [90,8;91,9;92,10;93,11];
%this_syn = zeros(size(tvec));
noiseAmp = [.01];
for ii = 1:numel(mag)
    for c = 1:4
        this_syn = zeros(size(tvec));
        this_sig{ii,c} = mag(ii)*cos(2*pi*f(c)*tvec);
        this_syn = this_syn + this_sig{ii,c};
        synData{ii}(c,:) = this_syn + noiseAmp*rand(1,length(this_syn));
    end
end
LFPTs.data = synData{1};
LFPTs.tvec = tvec;
LFPTs.label = {'NASL','NASR','NACL','NACR'};
LFPTsFilt = LFPTs;
eventTs.label = {'RSTART','RSTOP','Orientation','Approach (Start)','Approach (End)','Binge (Start)','Binge (End)','Rest (Start)','Rest (End)'};
eventTs.t = {[],[],[],[],[],[1;13],[5;17],[7;19],[11;23]};


for c = 1:4
    this_syn = zeros(size(tvec));
    for ii = 1:numel(mag)
        this_sig{ii} = mag(ii)*cos(2*pi*(f(c,ii))*tvec);
        %plot(tvec,this_sig,'r:'); hold on;
         this_syn = this_syn + this_sig{ii};
    end
    synData(c,:) = this_syn + noiseAmp*rand(1,length(this_syn));
end

%% Create test data for tabulate
nfile = 12;
for ii = 1:nfile
    files{ii} = strcat('syn',num2str(ii));
end
%%
psd = reshape([.1:.1:2],5,4);
coh = reshape([21:50],5,6);
y1 = reshape([1:20],5,4); y2 = reshape([21:50],5,6);
%%
for ii = 1:nfile
    PSDs.Base.Bands{ii,2} = psd;
    PSDs.Base.Bands{ii,1} = y1.*psd+psd;
    Cohs.rel.Base{ii,1} = coh;
    PSDs.Base.Bands{ii,3} = (PSDs.Base.Bands{ii,1}-PSDs.Base.Bands{ii,2})./abs(PSDs.Base.Bands{ii,2});
end
