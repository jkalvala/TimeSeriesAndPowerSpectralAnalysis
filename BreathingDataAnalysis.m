%%%%%%%%%%% Jo Suresh 2016
%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is written to analze breathing data with the following functionalities.
% V1 - basline detection based on bandpass filtering
% V2 - allows user to specify the start and stop times nedded for your analysis, 
%      from a long recording 
% V3 - added metrics like 20-80% riseTime, fallTime, width at half
%      amplitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
fileName = input('Enter the filename to analyze (program will add .abf extension as defualt):  ','s');
sr = input('Enter the sampling frequency in Hz:   ');
tStartTimeSecs = input('Enter START TIME in secs:   ');
tEndTimeSecs = input('Enter END TIME in secs:   ');

% Read the abf file
fileNameOpen = [fileName '.abf'];
[abfData,si,h]=abfload(fileNameOpen);
x_orig=abfData(:,1)';
t_orig=1/sr:1/sr:length(x_orig)/sr;
x_WithinTimeSlots = abfData(tStartTimeSecs*sr:tEndTimeSecs*sr,1)';
t = tStartTimeSecs+(1/sr):1/sr:tEndTimeSecs+(1/sr);

% Plot the signal
figure(1);
subplot(211)
plot(t_orig,abfData(:,1));
title('Entire Original Signal');
subplot(212);
plot(t,x_WithinTimeSlots);
xlim([tStartTimeSecs tEndTimeSecs]);
title(['Part of signal being analyzed (' num2str(tStartTimeSecs) 's - ' num2str(tEndTimeSecs) 's)']);
xlabel('Time in secs');
ylabel('Amplitude mV');

%% filter the signal
fN = sr/2; % sr - sampling frequency. fN - Nyquist freq.
%*************************************************
fls=0.4; % This is high pass filter setting in Hz
fhs = 3; % This is low pass filter setting in Hz
%*************************************************


[g,d]=butter(2,[fls/fN fhs/fN],'bandpass');
x_filt =filter(g,d,x_WithinTimeSlots);

figure(2)
hold all;
plot(t,x_WithinTimeSlots,'b');
plot(t,x_filt,'r');
grid on;
legend('OriginalSignal','FilteredSignal');
xlabel('Time in secs');
ylabel('Amplitude AU');

%% Find maxima/peaks
offset_time_in_sec = 2; % We need to ignore the first few secs due to filter artifacts. Hence the offset (in sec).
x_filt_offset = x_filt(offset_time_in_sec*sr:length(x_filt)); % getting rid of first 2 secs due to initial filter artifact
t_offset = t(offset_time_in_sec*sr:length(x_filt));

[~,locs_maxima] = findpeaks(x_filt_offset,1:1:length(x_filt_offset),'MinPeakDistance',0.3*sr,'MinPeakHeight',0); %Minimum distance between peaks is 0.2 secs

x_filt_offset_inverted = -x_filt_offset;
[~,locs_minima] = findpeaks(x_filt_offset_inverted,1:1:length(x_filt_offset),'MinPeakDistance',0.3*sr,'MinPeakHeight',0); %Minimum distance between peaks is 0.2 secs


% figure
% hold on
% plot(t_offset,x_filt_offset);
% plot((locs_maxima/sr)+offset_time_in_sec+tStartTimeSecs,x_filt_offset(locs_maxima),'rv','MarkerFaceColor','r');
% plot((locs_minima/sr)+offset_time_in_sec+tStartTimeSecs,x_filt_offset(locs_minima),'rs','MarkerFaceColor','b');
% grid on;
% legend('Filtered Signal','maxima','minima');


%% Pull out all the single waves
% We are only going to use all those maxima, that are enclosed
% between the first and the last minima.
first_index_of_maxima = find(locs_maxima>locs_minima(1),1); % find the first occurance of maxima that occurs after the first minima. This is where the first wave starts.
last_index_of_maxima = find(locs_maxima<locs_minima(length(locs_minima)),1,'last'); % find the first occurance of maxima that occurs after the first minima. This is where the first wave starts.
final_locs_maxima = locs_maxima(first_index_of_maxima:last_index_of_maxima);

figure(3)
hold on
plot(t_offset,x_filt_offset);
plot((final_locs_maxima/sr)+offset_time_in_sec+tStartTimeSecs,x_filt_offset(final_locs_maxima),'rv','MarkerFaceColor','r');
plot((locs_minima/sr)+offset_time_in_sec+tStartTimeSecs,x_filt_offset(locs_minima),'rs','MarkerFaceColor','b');
grid on;
legend('Filtered Signal','maxima','minima');
xlabel('Time in secs');
ylabel('Amplitude AU');
title('Bandpass Filtered Signal after applying 2 sec offset');



for i=1:length(final_locs_maxima)
    % capture the wave
    myWave{i}= x_filt_offset(locs_minima(i):locs_minima(i+1));
    myWaveTemp = x_filt_offset(locs_minima(i):locs_minima(i+1));
    myWaveRiseTime(i)= risetime(myWaveTemp,sr,'PercentReferenceLevels',[20 80],'StateLevels',[x_filt_offset(locs_minima(i)) x_filt_offset(final_locs_maxima(i))]);
    myWaveFallTime(i)= falltime(myWaveTemp,sr,'PercentReferenceLevels',[20 80],'StateLevels',[x_filt_offset(locs_minima(i+1)) x_filt_offset(final_locs_maxima(i))]);
    %myWaveWidthAtHalfAmp(i)= pulsewidth(myWaveTemp,sr,'StateLevels',[x_filt_offset(locs_minima(i+1)) x_filt_offset(final_locs_maxima(i))]);
  
end;


%% Periods
timeOfMaximaInSecs = final_locs_maxima/sr;
Periods = diff(timeOfMaximaInSecs);
Freq = 1./Periods;


%% Irregularity index
if length(Periods) < 3
    display('*************************************************************************************')
    display('ERROR - There are not enough maxima to calculate scores. Minimum required is 3 maxima');
    display('*************************************************************************************')
    
end
for i=1:length(Periods)-1
    Score(i) = (abs(Periods(i+1)- Periods(i))/Periods(i));
end;

% %% Extracting each wave separately to study the separate riseTime, fallTime and widthAtHalfAmp
% 
% % Step1: Invert the signal and find the maxima again. This effectively
% % captures the minima from the original signal.
% x_filt_offset_inverted = -x_filt_offset;
% [valleys,valleyLocs] = findpeaks(x_filt_offset_inverted,t_offset,'MinPeakDistance',0.2,'MinPeakHeight',0); %Minimum distance between peaks is 0.2 secs
% text(valleyLocs+.02,valleys,num2str((1:numel(valleys))'))


%% Writing the results to an excel sheet
filename = [fileName num2str(tStartTimeSecs) 's_' num2str(tEndTimeSecs) 's_Results.xlsx'];
header={'TimeOfPeak','Period','Freq','IrregularityScor','RiseTimeInSecs','FallTimeInSecs','WidthAtHalfAmpl'};
sheet = 2;
xlswrite(filename,header,sheet);
xlRange = 'A2';
timeCorrectionFactor = offset_time_in_sec+tStartTimeSecs;
xlswrite(filename,((timeOfMaximaInSecs)+timeCorrectionFactor)',sheet,xlRange);
xlRange = 'B2';
xlswrite(filename,Periods',sheet,xlRange);
xlRange = 'C2';
xlswrite(filename,Freq',sheet,xlRange);
xlRange = 'D2';
xlswrite(filename,Score',sheet,xlRange);
xlRange = 'E2';
xlswrite(filename,myWaveRiseTime',sheet,xlRange);
xlRange = 'F2';
xlswrite(filename,myWaveFallTime',sheet,xlRange);
% xlRange = 'G2';
% WidthAtHalfAmpl = cell2mat(myWaveWidthAtHalfAmp)
% xlswrite(filename,WidthAtHalfAmpl',sheet,xlRange);
%% THE END !
display('Data written to excel sheet, Complete !')
