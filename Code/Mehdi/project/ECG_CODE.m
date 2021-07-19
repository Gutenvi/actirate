%% start with clean workspace
close all;                  
clear all;
%% loading files and extracting ECG and PPG signals
r=readtable('0707-1630-A.csv');
vec=table2array(r(2:end,2));

ppg=flip(vec);

[hdr, record] = edfread('16-30-15.EDF');
ecg=flip(record(1,:));

ecg_fs=125;
ppg_fs=100;

ecg_time=[1:length(ecg)]./ecg_fs; 
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
wt = modwt(lp_ecg_sig, 'sym4', 7);
% Construct MRA matrix using modwtmra
mra = modwtmra(wt, 'sym4');
% Sum along selected multiresolution signals. In this command i am combining the true frequencies
filtered_ecg_sig = sum(mra(levelForReconstruction,:),1);
%% plotting signals over ecg_time
figure
plot(ecg_time,ecg)
hold on
plot(ecg_time,filtered_ecg_sig)
xlabel('ecg_time')
ylabel('Amplitude')
title 'Raw ECG and Filtered ECG'
legend('Raw ECG','Filtered ECG')
%% Periodogram
[Pxx,Freq] = periodogram(ecg,rectwin(length(ecg)),length(ecg),125);
figure
subplot(121)
   plot(Freq,10*log10(Pxx))
   grid on;
   xlabel('frequency (Hz)'); ylabel('power/frequency (dB/Hz)');
   title 'Periodogram of Raw ECG'

 [Pxx,Freq] = periodogram(filtered_ecg_sig,rectwin(length(filtered_ecg_sig)),length(filtered_ecg_sig),125);
 subplot(122)
   plot(Freq,10*log10(Pxx))
   grid on;
   xlabel('frequency (Hz)'); ylabel('power/frequency (dB/Hz)');
   ylim ([-20,70])
   title (['Periodogram of Filtered ECG'])
%% FFT of filtered ecg signal
l=length(filtered_ecg_sig);
nfft=2^nextpow2(l);
freq=ecg_fs/2*linspace(0,1,nfft/2+1);
f_xn=abs(fft(filtered_ecg_sig,nfft));
y = fft(filtered_ecg_sig,nfft); % Fast Fourier Transform
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
f_scale = (0:nfft/2)* 125/nfft; % frequency scale
f_dominant = f_scale(k);
figure

plot(freq , f_xn(1:length(freq)))
title 'FFT of Filtered ECG'
xline(f_dominant,'--r','linewidth',3)
title(['Dominant Frequency ', num2str(f_dominant), ' Hz'])

%% hanning window
Wn=0.5/ppg_fs;
n = 499; % No. of order
n_over = 60; 

b1_ecg = fir1(n, Wn, hann(n+1));
hnn_sig=filtfilt(b1_ecg,1,filtered_ecg_sig);
figure
plot(ecg_time,hnn_sig)
grid('on')
title(['hanning Window '])
%% hamming window
b2_ecg = fir1(n, Wn, hann(n+1));
hnn_sig=filtfilt(b2_ecg,1,filtered_ecg_sig);
figure
plot(ecg_time,hnn_sig)
grid('on')
title(['Hamming Window '])
%% Blackman window
b3_ecg = fir1(n, Wn, blackman(n+1));
blk_sig=filtfilt(b3_ecg,1,filtered_ecg_sig);
figure
plot(ecg_time,blk_sig)
grid('on')
title(['BlackMann Window '])
%% Extra


frequencies = [];
w=500;

for i = 1:(length(ecg)-w)
    window = zeros(length(ecg));
    window(i:i+w) = 1 ;       

windowed = filtered_ecg_sig*window;
[b,a]=butter(5,0.7/100/2,'high');
hp_win_sig=filtfilt(b,a,windowed);

% Lowpass filter
[b,a]=butter(5,4/100/2,'low');
windowed=filtfilt(b,a,hp_win_sig);
    
% FFT of filtered ecg signal
l=length(windowed);
nfft=2^nextpow2(l);
freq=ecg_fs/2*linspace(0,1,nfft/2+1);
f_xn=abs(fft(windowed,nfft));
% dominant frequency
x = windowed;
x = x - mean(x);                                            
nfft = 2^nextpow2(length(x)); % next larger power of 2
y = fft(x,nfft); % Fast Fourier Transform
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
f_scale = (0:nfft/2)* 125/nfft; % frequency scale
f_dominant_w = f_scale(k);
mp= (k-1)*125/nfft;
end
figure
plot(freq , f_xn(1:length(freq)))
title (['FFT of Filtered ECG'])
xline(f_dominant_w,'--r','linewidth',3)
title(['Dominant Frequency ', num2str(f_dominant_w), ' Hz'])


