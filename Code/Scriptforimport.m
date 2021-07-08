<<<<<<< HEAD
v=importdata("C:\Users\simon\OneDrive\Dokumente\GitHub\actirate\Data\ECG\1\0707-1630-A.csv");
ans= edfread("C:\Users\simon\OneDrive\Dokumente\GitHub\actirate\Data\ECG\1\16-30-15.EDF")
=======
v=importdata("..\Data\ECG\1\0707-1630-A.csv");
ans= edfread("..\Data\ECG\1\16-30-15.EDF")
>>>>>>> e78110fdadf879db7938ddcca23370e0cc15c98d
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
plot(v)
legend
hold off