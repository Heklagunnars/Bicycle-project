clc
clear all
close all

%% Upload the data


file = readmatrix('cord_v4.txt','TrimNonNumeric',true); %here it takes the non numeric values from the file

%Assigning the data
muscle = zeros(length(file),5);
muscle = file(:,3:7);
time = zeros(length(file),1);



%% Transfer function for emg
%the transfer function was gotten from the data sheet.
%EMG = (ADC/2^n)-0.5)*Vcc/Gain


time = file(:,1)./1000; %setting the time into seconds
for i = 1:5
    data(:,i) = ((((muscle(:,i)/2^16)-0.5)*3)/1000) *1000; %transfer function
    %muslce(i:2) = ((((file(i,4)/2^16)-0.5)*3)/1000) *1000;
    %muslce(i:3) = ((((file(i,5)/2^16)-0.5)*3)/1000) *1000;
    %muslce(i:4) = ((((file(i,6)/2^16)-0.5)*3)/1000) *1000;
end

for i=1:5
    figure(1)
    subplot(5,1,i)
    plot(time,data(:,i));
    xlabel('Time (s)');
    ylabel('Voltage (mV)');
end

%% fast fourier transform (FFT)

fs = 1000;          %sampling frequency
L = length(data);
f = fs*(0:(L/2))/L; %FFT freq range is defined as frequency reselution from 0 to 1/2 data lenght

figure;
for i = 1:5
    p1 = fft(data(:,i));      % taking the fast fouier transform
    p1 = abs(p1/L);     % divided by L to normalize

    p1 = p1(1:L/2+1);
    p1(2:end-1) = 2*p1(2:end-1);

    subplot(5,1,i)
    plot(f,p1);
    xlabel('Frequency(Hz)');
    ylabel('Intensity')
    zlabel('Glute')
    grid
    
    if(i==1)
        title('Glutes')
    elseif (i==2)
        title('Inner quad')
    elseif (i==3)
        title('Outer quad')
    elseif (i==4)
        title('Back thigh')
    else 
        title('Calve')
    end

    hold on

end 


%% 4th order butterworth filter

fnyq = fs/2; %Nyquist frequency
fcuthigh = 15; %This was dicided manually
fcutlow = 250; %This was dicided manually

[b,a] = butter(4,[fcuthigh,fcutlow]/fnyq,'bandpass'); % 4th Butterworth filter - is this the best filter??

for i = 1:5
    muscle(:,i) = filtfilt(b,a,muscle(:,i));
end 


%% Full wave rectification

rec_signal = zeros(length(muscle),5); %making a matix with zeros
 

for i=1:5
    rec_signal(:,i) = abs(muscle(:,i)); %putting the abs value for the data
end 
%% Plotting final signal

figure;
for i = 1:5
    subplot(5,1,i)
    plot(time, rec_signal(:,i));
    xlabel('Time[s]');
    ylabel('Voltage [mV]')
    title('Vöðvi 1')
    grid
    
    if(i==1)
        title('Glutes')
    elseif (i==2)
        title('Inner quad')
    elseif (i==3)
        title('Outer quad')
    elseif (i==4)
        title('Back thigh')
    else 
        title('Calve')
    end
    xlim([0,40])
end
sgtitle('Final plot')
hold off


