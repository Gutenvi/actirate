% This function filters the ECG signals to denoise
% Input parameter 1: ecg (raw ECG data)
% Input parameter 2: ecg_time (time points for the plot)

% Output parameter: ecg_filt (filtered ECG)
function filtered_ecg_sig = filterECG(ecg, ecg_time)

Fs = 125; %Sampling Frequency


%%%%% Filter Design Parameters Defined for ECG %%%%%
fcomb = [[0.5 1.0], [45 48 52 55]]; % Assuming the data is acquired in Europe
mags = [[0 1], [0 1]];
dev = [[0.5 0.1], [0.1 0.5]];
%%%%% Filter Design Parameters Defined for ECG %%%%%


%%%% Design kaiser window filter %%%%
[n,Wn,beta,ftype] = kaiserord(fcomb,mags,dev,Fs);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
%%%% Design kaiser window filter %%%%

filtered_ecg_sig = filtfilt(hh, 1, ecg); % Apply the designed filter on the input data

%%%% PLOT %%%%%
figure
plot(ecg_time, ecg)
grid
title('Original Signal')
hold on
plot(ecg_time, filtered_ecg_sig)
grid
title('Filtered Signal')
title 'Raw ECG and Filtered ECG'
legend('Raw ECG','Filtered ECG')
hold off

end
%%%% PLOT %%%%%


