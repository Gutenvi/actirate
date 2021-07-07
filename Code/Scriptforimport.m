v=importdata("rct3configurator\0506-16-2-S_ref.csv");
ans= edfread("G:\DATA\20210705\16-07-17.EDF")
ECG=0;
for i=1:16
ECG=[ECG;ans.ECG{i,1}]
end
ECG=ECG(3:2001);
ECG=ECG-mean(ECG);
v=v(2:2000);
v=v-mean(v);
plot(ECG)
hold
%plot(v)
legend
hold off