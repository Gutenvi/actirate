% This the code for first data recorded with the ECG 
%% start with clean workspace
close all;                  
clear all;
%% loading files and extracting ECG and PPG signals
%importing the data into matlab
[hdr, record] = edfread('16-30-15.EDF');
ecg=flip(record(1,:));
ecg=ecg - mean(ecg);
%sampling frequency
ecg_fs=125;
%time scaling
%ecg_time=[1:length(ecg)]./ecg_fs;#
ecg_time=linspace(0,length(ecg)/ecg_fs,length(ecg)) ;
%% filtering
% Bandpass filter
%[b,a]=butter(5,[1 120]/125/2,'bandpass');
%bp_ecg_sig=filtfilt(b,a,ecg);

%Highpass filter
[b,a]=butter(5,1/125/2,'high');
hp_ecg_sig=filtfilt(b,a,ecg);

% Lowpass filter
[b,a]=butter(5,120/125/2,'low');
lp_ecg_sig=filtfilt(b,a,hp_ecg_sig);

% removing very low frequency noise that stoping signal to remain on zero line
% using wavelet transform
% Logical array for selecting reconstruction elements
%for this process, i used signal multiresolution app of MATLAB. 
%This app of MATLAB, divide and display different frequencies of signal.
%in this part, i divided the ECG frequencies in 8 different parts. 
%then i stopped 2 lowest frequencies by saying 
%"false false" which were creating deviation in ecg signal,
%and passing rest of the frequencies by saying "pass pass pass pass pass pass"
%in this command i am saying that i want 6 frequencies out of 8

levelForReconstruction = [true, true, true, true, true, true, false, false];

%in this command i am commanding the MATLAB to break the signals into 7 parts
% Perform the decomposition using modwt. 
wt = modwt(lp_ecg_sig, 'sym4', 6);
% Construct MRA matrix using modwtmra
mra = modwtmra(wt, 'sym4');
% Sum along selected multiresolution signals. In this command i am combining the true frequencies
filtered_ecg_sig = sum(mra(levelForReconstruction,:),1);
%% plots of the raw and filtered data
figure
plot(ecg_time,ecg)
hold on
plot(ecg_time,filtered_ecg_sig)
xlabel('ecg_time')
ylabel('Amplitude')
title 'Raw ECG and Filtered ECG'
legend('Raw ECG','Filtered ECG')
%% Periodogram
[Pxx,Freq] = periodogram(ecg,flattopwin(length(ecg)),length(ecg),125);
figure
subplot(121)
   plot(Freq,sqrt(Pxx))
   grid on;
   xlabel('frequency (Hz)'); ylabel('power/frequency (dB/Hz)');
   title 'Periodogram of Raw ECG'

 [Pxx,Freq] = periodogram(filtered_ecg_sig,flattopwin(length(filtered_ecg_sig)),length(filtered_ecg_sig),125);
 subplot(122)
   plot(Freq,sqrt(Pxx))
   grid on;
   xlabel('frequency (Hz)'); ylabel('power/frequency (dB/Hz)');
   
   title (['Periodogram of Filtered ECG'])
%% FFT of filtered HS signal and the maximum frequency
x = filtered_ecg_sig;
x = x - mean(x);                                            
nfft = 2^nextpow2(length(x)); % next larger power of 2
%y = fft(x,nfft); % Fast Fourier Transform
y = fft(x); % Fast Fourier Transform
y = abs(y.^2); % raw power spectrum density
y_hs = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y_hs); % find maximum
f_scale_hs = (0:nfft/2)* 125/nfft; % frequency scale
f_dominant_hs = f_scale_hs(k);
figure
plot(f_scale_hs, y_hs)
xlim([0 10]);
grid('on')
title(['Dominant Frequency ', num2str(f_dominant_hs), ' Hz'])
xline(f_dominant_hs,'--r','linewidth',3)
figure
plot(f_scale_hs, y_hs)
xlim([0 10]);
title (['FFT of Filtered ECG'])



%% Windowing

w=500;
%f_dominant_hs=zeros(1, length(ppg_shifted)-(w+0));
f_dominant_hs=[];
r=length(ecg)-w;
for i = 1:r
    window = zeros(1,length(ecg));
    window(i:i+w) = 1 ;       
windowed = filtered_ecg_sig.*window;
%taking the mean considering the zero
windowed(i:w+1) = windowed(i:w+1) - mean(windowed(i:w+1)); 
%high pass 
[b,a]=butter(5,1/125/2,'high');

hp_win_sig=filtfilt(b,a,windowed);

% Lowpass filter
[b,a]=butter(5,120/125/2,'low');
x=filtfilt(b,a,hp_win_sig);  

nfft = 2^nextpow2(length(x)); % next larger power of 2
%y = fft(x,nfft); % Fast Fourier Transform
y = fft(x); % Fast Fourier Transform

y = abs(y.^2); % raw power spectrum density
y= y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
f_scale_hs = (0:nfft/2)*125/nfft; % frequency scale
f_dominant_hs(i) = f_scale_hs(k);

end
MaxiF=max(f_dominant_hs);
t=length(f_dominant_hs);
y1=f_dominant_hs*60;
figure
plot((1:t)/100,y1)