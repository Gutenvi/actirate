%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CaCom
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          % ALWAYS start with clean workspace
close all;                  
clear all;
fontSize=8;
%%
data = readtable('1906-09-S70.csv');  % skips the first 5 rows of data
x= data{:,2};
x(1:1:6)=[];
x = x-mean(x);
%%
order = 4; 
sampling_freq = 100;
cutoff_freq_up = 4;
cutoff_freq_down = 0.5;
number_of_samples = length(x);
z=(1/sampling_freq)*number_of_samples;

time=linspace(0, z, number_of_samples);
z_f=1./time;
z_f(1)=[];
%% Filtrage
[b,a]=butter(order,[cutoff_freq_down,cutoff_freq_up]/(sampling_freq/2),'bandpass');
filtsig=filter(b,a,x);  %filtered signal
figure(1)
plot(time,filtsig);
hold on
plot(time,x);
hold off
xlabel('Time ');
ylabel('Amplitude ');
%% FFT
FFT = fft(filtsig)/length(filtsig);
n = length(filtsig)/2;
spectrum_signal = FFT(1:floor(n)+1);
spectrum_signal(2:end-1) = spectrum_signal(2:end-1)*2;
k=(0:floor(n))/n * sampling_freq/2;
v = abs(spectrum_signal);
figure(2)
stem(k,v);
xlabel('Frequency ');
ylabel('Amplitude ');
%print(gcf,'FFT','-depsc');
%% Periodogram
FFT_x = fft(x)/length(x);
n = length(x)/2;
spectrum_signal = FFT_x(1:floor(n)+1);
spectrum_signal=20*log10(spectrum_signal);
spectrum_signal(2:end-1)   = spectrum_signal(2:end-1);
f=(0:floor(n))/n * sampling_freq/2;
figure(4)
semilogx(f,spectrum_signal);
xlim([1 50])
xlabel('Frequency [Hz]');
ylabel('Linear spectrum [V RMS]');
%%
FFT_x = fft(filtsig)/length(filtsig);
n = length(filtsig)/2;
spectrum_signal = FFT_x(1:floor(n)+1);
spectrum_signal=20*log10(spectrum_signal);
spectrum_signal(2:end-1)   = spectrum_signal(2:end-1);
f=(0:floor(n))/n * sampling_freq/2;
figure(6)
semilogx(f,spectrum_signal);
xlim([1 50])
xlabel('Frequency [Hz]');
ylabel('Linear spectrum [V RMS]');
%%
