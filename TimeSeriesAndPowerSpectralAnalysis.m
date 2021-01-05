%%%%%%%%%%% Jo Suresh 2016
%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is written to detect spikes from time-series data and run power spectrum analysis using FFT
% V1 - basline detection based on bandpass filtering
% Input the filename which contains time series data at a sample rate of 25kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
sr=25000;
Nyq_freq=sr/2;

fnrawOrig=input('Filename :', 's');
%badChannels=input('Enter bad channels in an aray format eg., [1 2 3]: ');
t = input('Enter time in secs: ');

le=t*sr; % Length of the file
fN=sr/2;

fls=300;            % filter cutoff spikes
fhs=1500;

flf=1;              % filter cutoff local field potential (low freq)
fhf=5;

% integrator settings - Low frequency range
RC=50e-3;
f=1/(2*pi*RC);      % f(-3dB) for 50 ms integrator
UIR=(1/sr)*(1/RC);  % Unit Impulse response Amplitude

% filter coefficients
[b,a]=butter(2,[fls/fN fhs/fN]); % 300 - 1500 Hz
[g,d]=butter(1,f/fN,'low'); % 1 - 3 Hz
[e,f]=butter(2,[flf/fN fhf/fN]); % not used in this prog

% Calculating Instantaneous firing rate - averaged across all channels.
ch = 1:1:60;
fnraw = [fnrawOrig ];
fid=fopen(fnraw);
clear RAW filtLF filtHF FA SPIKE spikeCountPerBin BinnedLFP LFP;
RAW=fread(fid,[length(ch),le],'double');

LFP =  zeros(length(ch),le); % le - no of sample points, for2 mins data it is 3000000
HFP = zeros(length(ch),le);
SPIKES = zeros(length(ch),le);
for c=ch
   
    display(['Working on channel# ' num2str(c)]);
    %Step 1 - demean signal
    RAW(c,:) = RAW(c,:)- mean(RAW(c,:));
    %Step 2 - Extract Low freq
    LFP(c,:) = filter(g,d,RAW(c,:));
    %Step 3 - Extract High freq
    HFP(c,:) = filtfilt(b,a,RAW(c,:));
    FA =  HFP(c,:);
    
    th_neg=std(FA).*-5;       % threshold for spike detection - 5*Std deviation
    
    % Detecting spikes
    iprev=round(-sr*.002);
    tmp=-iprev;
    iprev_pos=round(sr*.002);
    tmp_pos = iprev_pos;
    for i = 1:length(FA )
      if ((FA(i) <th_neg)&&(i>iprev+tmp)&&((i+tmp)<length(FA ))&&(i-round(tmp/4)>0));
            SPIKES(c,i)=1;
            iprev=i;
        end;
    end;
end; % end of #channels


%% Calculating power spectrum of averaged raw signal
raw = mean(RAW,1);
Fs = 25000;
x = raw;
N=length(x); %get the number of points
k=0:N-1;     %create a vector from 0 to N-1
T=N/Fs;      %get the frequency interval
freq=k/T;    %create the frequency range
X=fft(x)/N; % normalize the data

%only want the first half of the FFT, since it is redundant
cutOff = ceil(N/2);

%take only the first half of the spectrum
X = X(2:cutOff);
freq = freq(2:cutOff);
Power = X.*conj(X);

%% Plot data

figure
subplot(211);
plot(1/sr:1/sr:t,mean(RAW,1));
subplot(212);
plot(1/sr:1/sr:t,mean(LFP,1));

%Low freq spectrum
figure;
hold all
subplot(211);
plot(freq(1:50) ,Power(1:50));
xlabel('Freq in Hz');
ylabel('Power');
title('Power Spectrum using FFT');

%High freq spectrum
subplot(212);
plot(freq(100:15000) ,Power(100:15000));
xlabel('Freq in Hz');
ylabel('Power');
title('Power Spectrum using FFT');


%%
fclose('all');




