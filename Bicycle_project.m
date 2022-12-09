clc
clear all
close all

%% Upload the data


file = readmatrix('minute 1.txt','TrimNonNumeric',true); %here it takes the non numeric values from the file

%Assigning the data
muscles = zeros(length(file),5);
muscles = file(:,3:8);
time = zeros(length(file),1);



%% Transfer function for emg
%the transfer function was gotten from the data sheet.
%EMG = (ADC/2^n)-0.5)*Vcc/Gain

time = file(:,1)./1000; %setting the time into seconds

for i = 1:6
    data(:,i) = ((((muscles(:,i)/2^16)-0.5)*3)/1000) *1000; %transfer function
end


%% fast fourier transform (FFT)

fs = 1000;          %sampling frequency
L = length(data);
f = fs*(0:(L/2))/L; %FFT freq range is defined as frequency reselution from 0 to 1/2 data lenght

for i = 1:6
    p1 = fft(data(:,i));      % taking the fast fouier transform
    p1 = abs(p1/L);     % divided by L to normalize

    p1 = p1(1:L/2+1);
    p1(2:end-1) = 2*p1(2:end-1);

end 


%% 4th order butterworth filter

fnyq = fs/2; %Nyquist frequency
fcuthigh = 15; % Gotten from FFT
fcutlow = 250; %Gotten from the FFT

[b,a] = butter(4,[fcuthigh,fcutlow]/fnyq,'bandpass'); % 4th Butterworth filter - is this the best filter??

for i = 1:6
    muscles(:,i) = filtfilt(b,a,muscles(:,i));
end 


%% Full wave rectification (but not for the tracker)

rec_signal = zeros(length(muscles),5); 
 

for i=1:5
    rec_signal(:,i) = abs(muscles(:,i)); 
end 

rec_signal(:,6) = muscles(:,6);

%% Root mean square (not for the tracker)

for i =1:5
    rec_signal(:,i) = sqrt(movmean(rec_signal(:,i).^2, 10));
end


%% Taking the standard deviation from the signal to get rid of certain frequencys

s = std(rec_signal);

for i=1:5
    rec_signal(:,i)=rec_signal(:,i)-s(i);
end


% Making the minus values zero

for i = 1:length(rec_signal)
    for j = 1:5
        if (rec_signal(i,j) < 0)
            rec_signal(i,j) = 0;
            
        end
    end
end


%% Getting exactly one loop

% Finding the locations of the peaks (loops)
[pks,locs]= findpeaks(rec_signal(:,6),time);
avgx = max(pks)-std(pks);
[pks,locs] = findpeaks(rec_signal(:,6),time,'MinPeakProminence',avgx);

%Extract the data from one loop to another

for i=1:length(locs)-1
    loop{i}=rec_signal(locs(i)*1000:locs(i+1)*1000,:);
end 



%% Final plot - both activation and golden loop


oneloop = loop{5}; 
muscle1 = muscle1(:,1);
th1 =linspace(0,2*pi,length(muscle1));
subplot(2,3,1)
h = polarscatter(th1,muscle1,'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [0,pi/2], 'BinCounts',[max(muscle1)], 'FaceColor','red')
title('Glute')

subplot(2,3,2)
muscle2 = oneloop(:,2);
th2 =linspace(0,2*pi,length(muscle2));
h = polarscatter(th2,muscle2, 'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [deg2rad(85),deg2rad(175)], 'BinCounts',[max(muscle2)], 'FaceColor','magenta')
title('Medial quadriceps')

subplot(2,3,3)
muscle3 = oneloop(:,3);
th3 =linspace(0,2*pi,length(muscle3));
h = polarscatter(th3,muscle3, 'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [deg2rad(85),deg2rad(175)], 'BinCounts',[max(muscle3)], 'FaceColor','blue')
title('Lateral quadriceps')

subplot(2,3,4)
muscle4 = oneloop(:,4);
th4 =linspace(0,2*pi,length(muscle4));
h = polarscatter(th4,muscle4, 'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [deg2rad(225),deg2rad(270)], 'BinCounts',[max(muscle4)], 'FaceColor','green')
title('Hamstring')

subplot(2,3,5)
muscle5 = oneloop(:,5);
th5 =linspace(0,2*pi,length(muscle5));
h = polarscatter(th5,muscle5, 'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [deg2rad(160),deg2rad(180)], 'BinCounts',[max(muscle5)], 'FaceColor','cyan')
title('Calf')














