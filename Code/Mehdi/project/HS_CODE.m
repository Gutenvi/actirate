%% start with clean workspace
close all;                  
clear all;
clc
%% loading files and extracting ECG and PPG signals
r=readtable('0707-1630-A.csv');
vec=table2array(r(5:end,2));

ppg_shifted=flip(vec)';

ppg=[];
for i=1:length(ppg_shifted)
    sub=ppg_shifted(i)-ppg_shifted(1);
    ppg=[ppg sub];
end

ppg_fs=100;

ppg_time=[1:length(ppg)]./ppg_fs;

%% filtering
% Highpass filter
[b,a]=butter(5,0.5/100/2,'high');
hp_ppg_sig=filtfilt(b,a,ppg);

% Lowpass filter
[b,a]=butter(5,4/100/2,'low');
lp_ppg_sig=filtfilt(b,a,hp_ppg_sig);

% removing very low frequency noise that stoping signal to remain on zero line
% using wavelet transform

% removing very low frequency noise that stoping signal to remain on zero line
% using wavelet transform
% Logical array for selecting reconstruction elements
%for this process, i used signal multiresolution app of MATLAB. 
%This app of MATLAB, divide and display different frequencies of signal.
%in this part, i divided the ppg frequencies in 8 different parts. 
%then i stopped 2 lowest frequencies by saying 
%"false false" which were creating deviation in ecg signal,
%and passing rest of the frequencies by saying "pass pass pass pass pass pass"
%in this command i am saying that i want 6 frequencies out of 8
levelForReconstruction = [true, true, true, true, true, true, true, false];
%in this command i am commanding the MATLAB to break the signals into 7 parts
% Perform the decomposition using modwt. 
wt = modwt(lp_ppg_sig, 'sym4', 7);
% Construct MRA matrix using modwtmra
mra = modwtmra(wt, 'sym4');
% Sum along selected multiresolution signals
filtered_ppg_sig = sum(mra(levelForReconstruction,:),1);


%% plotting signals over time
figure
plot (ppg_time,ppg)
hold on
%subplot(212)
plot(ppg_time,filtered_ppg_sig)
xlabel('time')
ylabel('Amplitude')
title 'Raw PPG and Filtered PPG'
legend('Raw PPG','Filtered PPG')

%% Periodogram
[Pxx,Freq] = periodogram(ppg,rectwin(length(ppg)),length(ppg),100);
figure
  semilogy(Freq,10*log10(Pxx))
   grid on;
   xlabel('frequency (Hz)'); ylabel('power/frequency (dB/Hz)');
   title 'Periodogram of Raw PPG'

 [Pxx,Freq] = periodogram(filtered_ppg_sig,rectwin(length(filtered_ppg_sig)),length(filtered_ppg_sig),100);
figure
  semilogy(Freq,10*log10(Pxx))
   grid on;
   xlabel('frequency (Hz)'); ylabel('power/frequency (dB/Hz)');
   
   title 'Periodogram of Filtered PPG'
   
%% FFT of filtered ecg signal
l=length(filtered_ppg_sig);
nfft=2^nextpow2(l);
freq=ppg_fs/2*linspace(0,1,nfft/2+1)
f_xn=abs(fft(filtered_ppg_sig,nfft));
y = fft(filtered_ppg_sig,nfft); % Fast Fourier Transform
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
f_scale = (0:nfft/2)* 125/nfft; % frequency scale
f_dominant = f_scale(k);
figure

plot(freq , f_xn(1:length(freq)))
title 'FFT of Filtered PPG'
xline(f_dominant,'--r','linewidth',3)
title(['Dominant Frequency ', num2str(f_dominant), ' Hz'])


%% Dominant frequency (optional)
%{
x = filtered_ppg_sig;
x = x - mean(x);                                            
nfft = 2^nextpow2(length(x)); % next larger power of 2
y = fft(x,nfft); % Fast Fourier Transform
y = abs(y.^2); % raw power spectrum density
y = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
f_scale = (0:nfft/2)* 125/nfft; % frequency scale
f_dominant = f_scale(k);

figure

plot(f_scale, y)

grid('on')
title(['Dominant Frequency ', num2str(f_dominant), ' Hz'])
xline(f_dominant,'--r','linewidth',3)
%}
%% hanning window
Wn=0.5/ppg_fs;
n = 499; % No. of order
n_over = 60; 

b1 = fir1(n, Wn, hann(n+1));
hnn_sig=filtfilt(b1,1,filtered_ppg_sig);
figure
plot(ppg_time,hnn_sig)
grid('on')
title(['hanning Window '])


%% hamming window
b2 = fir1(n, Wn, hamming(n+1)); 
hmm_sig=filtfilt(b2,1,filtered_ppg_sig);
figure
plot(ppg_time,hmm_sig)
grid('on')
title(['Hamming Window '])
%% Blackman window
b3 = fir1(n, Wn, blackman(n+1));
blk_sig=filtfilt(b3,1,filtered_ppg_sig);
figure
plot(ppg_time,blk_sig)
grid('on')
title(['Hamming Window '])
%% Extra Window (not done)
frequencies = [];
w=500;

for i = 1:(length(ppg)-w)
    window = zeros(length(ppg));
    window(i:i+w) = 1 ;       

windowed = filtered_ppg_sig*window;
[b,a]=butter(5,0.7/100/2,'high');
hp_win_sig=filtfilt(b,a,windowed);

% Lowpass filter
[b,a]=butter(5,4/100/2,'low');
windowed=filtfilt(b,a,hp_win_sig);
    
% FFT of filtered ecg signal
l=length(windowed);
nfft=2^nextpow2(l);
freq=ppg_fs/2*linspace(0,1,nfft/2+1);
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
title (['FFT of Filtered PPG'])
xline(f_dominant_w,'--r','linewidth',3)
title(['Dominant Frequency ', num2str(f_dominant_w), ' Hz'])

