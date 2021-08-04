% This the code for first data recorded with the HRS on the Heart while
% doing Sport
close all;                  
clear all;
clc
%% loading files and extracting HRS and PPG signals
%r=readtable('0707-1630-A.csv');
r=readtable('data.csv');
ppg_shifted=table2array(r(7:end,2));

%ppg_shifted=flip(vec)';
%{
ppg=[];
for i=1:length(ppg_shifted)
    sub=ppg_shifted(i)-ppg_shifted(1);
    ppg=[ppg sub];
end
%}
ppg=ppg_shifted - mean(ppg_shifted);
ppg_fs=100;
%x2=[1:length(ppg)]./ppg_fs
%ppg_time=[1:length(ppg)]./ppg_fs;
ppg_time = linspace(0,length(ppg_shifted)/ppg_fs,length(ppg_shifted)) ;

%% filtering
% Highpass filter

[b,a]=butter(5,0.4/100/2,'high');
hp_ppg_sig=filtfilt(b,a,ppg);

% Lowpass filter

[b,a]=butter(5,5/100/2,'low');
filtered_ppg_sig=filtfilt(b,a,hp_ppg_sig);
% removing very low frequency noise that stoping signal to remain on zero line
% using wavelet transform

% removing very low frequency noise that stoping signal to remain on zero line
% using wavelet transform
% Logical array for selecting reconstruction elements
%for this process, i used signal multiresolution app of MATLAB. 
%This app of MATLAB, divide and display different frequencies of signal.
%in this part, i divided the ppg frequencies in 8 different parts. 
%then i stopped 2 lowest frequencies by saying 
%"false false" which were creating deviation in HRS signal,
%and passing rest of the frequencies by saying "pass pass pass pass pass pass"
%in this command i am saying that i want 6 frequencies out of 8
levelForReconstruction = [true, true, true, true, true, true, false, false];
%in this command i am commanding the MATLAB to break the signals into 7 parts
% Perform the decomposition using modwt. 
wt = modwt(filtered_ppg_sig, 'sym4', 7);
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
title (['Raw and Filtered HBS Signal'])
legend('Raw HBSS','Filtered HBSS')
print(gcf,'Raw HRS and Filtered HRS - Heart - Sport','-depsc');
saveas(gcf,'Raw HRS and Filtered HRS - Heart - Sport.png')
%% Periodogram

[Pxx,Freq] = periodogram(ppg,flattopwin(length(ppg)),length(ppg),100);
figure
subplot(121)
   semilogy(Freq,sqrt(Pxx))
   
   grid on;
   xlabel('frequency (Hz)'); ylabel('power/frequency (dB/Hz)');
   title 'Periodogram of Raw PPG'

 [Pxx,Freq] = periodogram(filtered_ppg_sig,flattopwin(length(ppg)),length(ppg),100);
 subplot(122)
   semilogy(Freq,sqrt(Pxx))
    
   grid on;
   xlabel('frequency (Hz)'); ylabel('power/frequency (dB/Hz)');
  
   title 'Periodogram of Filtered PPG'
   
print(gcf,'Periodogram of Filtered HRS - Heart - Sport','-depsc');
saveas(gcf,'Periodogram of Filtered HRS - Heart - Sport.png')
%% FFT of filtered HS signal and the maximum frequency
x = filtered_ppg_sig;
x = x - mean(x);                                            
nfft = 2^nextpow2(length(x)); % next larger power of 2
%y = fft(x,nfft); % Fast Fourier Transform
y = fft(x); % Fast Fourier Transform
y = abs(y.^2); % raw power spectrum density
y_hs = y(1:1+nfft/2); % half-spectrum
[v,k] = max(y_hs); % find maximum
f_scale_hs = (0:nfft/2)* 100/nfft; % frequency scale
f_dominant_hs = f_scale_hs(k);
figure
plot(f_scale_hs, y_hs)
xlim([0 10]);
grid('on')
title(['Dominant Frequency ', num2str(f_dominant_hs), ' Hz'])
print(gcf,'FFT of Filtered HRS - Heart - Sport','-depsc');
saveas(gcf,'FFT of Filtered HRS - Heart - Sport.png')
xline(f_dominant_hs,'--r','linewidth',3)
figure
plot(f_scale_hs, y_hs)
xlim([0 10]);
title (['FFT of Filtered HBSS'])
%% Extra Window

w=500;
%f_dominant_hs=zeros(1, length(ppg_shifted)-(w+0));
f_dominant_hs=[];
r=length(ppg_shifted)-w;
for i = 1:r
    window = zeros(1,length(ppg_shifted));
    window(i:i+w) = 1 ;       
windowed = filtered_ppg_sig.*window;
%taking the mean considering the zero
windowed(i:w+1) = windowed(i:w+1) - mean(windowed(i:w+1)); 
%high pass 
[b,a]=butter(5,0.7/100/2,'high');

hp_win_sig=filtfilt(b,a,windowed);

% Lowpass filter
[b,a]=butter(5,13/100/2,'low');
x=filtfilt(b,a,hp_win_sig);  

nfft = 2^nextpow2(length(x)); % next larger power of 2
%y = fft(x,nfft); % Fast Fourier Transform
y = fft(x); % Fast Fourier Transform

y = abs(y.^2); % raw power spectrum density
y= y(1:1+nfft/2); % half-spectrum
[v,k] = max(y); % find maximum
f_scale_hs = (0:nfft/2)*100/nfft; % frequency scale
f_dominant_hs(i) = f_scale_hs(k);

end
MaxiF=max(f_dominant_hs);

t=length(f_dominant_hs);

y1=f_dominant_hs*60;

figure
plot((1:t)/100,y1)
print(gcf,'Windowed HRS - Heart - Sport','-depsc');
saveas(gcf,'Windowed HRS - Heart - Sport - Heart - Sport.png')
