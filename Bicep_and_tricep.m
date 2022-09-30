clc
clear all
close all

%% Upload the data

file = readmatrix('bicep+tricep_6var.txt','TrimNonNumeric',true);

time = zeros(length(file),2);
muscle = zeros(length(file),2);

muscle = file(:,3:4);
%muscle1 = file(:,3);
%muscle2 = file(:,4);
%muscle3 = file(:,5);
%muscle4 = file(:,6);


%% Plotting the raw data

% figure(1)
% plot(gogn(:,1),gogn(:,3))
% xlim([0,1000])

%% Transfer function for emg
%the transfer function was gotten from the data sheat.
%EMG = (ADC/2^n)-0.5)*Vcc/Gain

time = file(:,1)./1000; %setting the time into seconds

for i=1:2
data(:,i) = ((((muscle(:,i)/2^16)-0.5)*3)/1000)*1000; %transfer function
end

for i=1:2
figure(1)
subplot(2,1,i)
plot(time,data(:,i));
xlabel('Time (s)');
ylabel('Voltage (mV)');
end

%% fast fourier transform (FFT)

fs = 1000;          %sampling frequency
L = length(data);   %Length of the matrix
f = fs*(0:(L/2))/L; %FFT freq range is defined as frequency reselution from 0 to 1/2 data lenght


for i=1:2
p1 = fft(data(:,i));     % taking the fast fouier transform
p1 = abs(p1/L);     % divided by L to normalize

p1 = p1(1:L/2+1);

p1(2:end-1) = 2*p1(2:end-1); %??? why (compute a single sided spectrum by taking the positive part of double sided spectrum and multiply by 2

figure(3)
subplot(2,1,i)
plot(f,p1)
xlabel('Frequency')
ylabel('Intensity')
title('FFT')
end

%% Band pass filter

fnyq = fs/2; %Nyquist frequency
fcuthigh = 15; %This was dicided manually
fcutlow = 300; %This was dicided manually

[b,a] = butter(4,[fcuthigh,fcutlow]/fnyq, 'bandpass'); % 4th Butterworth filter
for i=1:2
    data(:,i) = filtfilt(b,a,data(:,i));
    figure(5)
    subplot(2,1,i)
    plot(time,data(:,i))
    xlabel('Time [s]')
    ylabel('Voltage [mV]')
end


%% Full wave rectification

rec_signal = zeros(length(data),1); %making a matix with zeros
for i=1:2
    rec_signal(:,i) = abs(data(:,i)); %putting the abs value for the data 
end

% Plotting
for i= 1:2
    figure(5)
    subplot(2,1,i)
    plot(time,rec_signal(:,i))
    xlabel('Time [s]')
    ylabel('Voltage [mV]')
%xlim([0,5])
end