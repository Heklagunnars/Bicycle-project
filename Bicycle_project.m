clc
clear all
close all

%% Upload the data

file = readmatrix('bicep.txt','TrimNonNumeric',true); %here it takes the non numeric values from the file

time = zeros(length(file),1);
muscle = zeros(length(file),4);
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

L = length(muscle)
time = file(:,1)./1000; %setting the time into seconds
 
for i = 1:L
    muslce(i:1) = ((((file(i,3)/2^16)-0.5)*3)/1000) *1000; %transfer function
    muslce(i:2) = ((((file(i,4)/2^16)-0.5)*3)/1000) *1000;
    muslce(i:3) = ((((file(i,5)/2^16)-0.5)*3)/1000) *1000;
    muslce(i:4) = ((((file(i,6)/2^16)-0.5)*3)/1000) *1000;
end
%% Plotting the raw EMG data


%% fast fourier transform (FFT)

fs = 1000;          %sampling frequency

f = fs*(0:(L/2))/L; %FFT freq range is defined as frequency reselution from 0 to 1/2 data lenght

figure;
for i = 1:4
    p1 = fft(muscle(:,i));      % taking the fast fouier transform
    p1 = abs(p1/L);     % divided by L to normalize

    p1 = p1(1:L/2+1);
    p1(2:end-1) = 2*p1(2:end-1);

    subplot(4,1,i)
    plot(f,p1)
    xlabel('Frequency(Hz)');
    ylabel('Intensity')
    label('Vöðvi 1')
    grid
    
    if(i==1)
        title('Vöðvi 1')
    elseif (i==2)
        title('Vöðvi 2')
    elseif (i==3)
        title('Vöðvi 3')
    else 
        title('Vöðvi 4')
    end

    hold on

end 


%% 4th order butterworth filter

fnyq = fs/2; %Nyquist frequency
fcuthigh = 15; %This was dicided manually
fcutlow = 300; %This was dicided manually

[b,a] = butter(4,[fcuthigh,fcutlow]/fnyq,'bandpass'); % 4th Butterworth filter - is this the best filter??

for i = 1:4
    muscle(:,i) = filtfilt(b,a,muscle(:,i));
end 


%% Full wave rectification

rec_signal = zeros(length(muscle),4); %making a matix with zeros
 

for i=1:4
    rec_signal(:,i) = abs(muscle(:,i)); %putting the abs value for the data
end 
%% Plotting final signal

figure;
for i = 1:4
    subplot(4,1,i)
    plot(time, rec_signal(:,i));
    xlabel('Time[s]');
    ylabel('Voltage [mV]')
    label('Vöðvi 1')
    grid
    
    if(i==1)
        title('Vöðvi 1')
    elseif (i==2)
        title('Vöðvi 2')
    elseif (i==3)
        title('Vöðvi 3')
    else 
        title('Vöðvi 4')
    end
end
sgtitle('Final plot')
hold off


