clc
clear all
close all

%% Upload the data


file = readmatrix('circuit_bicycle.txt','TrimNonNumeric',true); %here it takes the non numeric values from the file

%Assigning the data
muscle = zeros(length(file),5);
muscle = file(:,3:8);
time = zeros(length(file),1);



%% Transfer function for emg
%the transfer function was gotten from the data sheet.
%EMG = (ADC/2^n)-0.5)*Vcc/Gain

time = file(:,1)./1000; %setting the time into seconds

for i = 1:6
    data(:,i) = ((((muscle(:,i)/2^16)-0.5)*3)/1000) *1000; %transfer function
end


%% fast fourier transform (FFT)

fs = 1000;          %sampling frequency
L = length(data);
f = fs*(0:(L/2))/L; %FFT freq range is defined as frequency reselution from 0 to 1/2 data lenght

figure;
for i = 1:6
    p1 = fft(data(:,i));      % taking the fast fouier transform
    p1 = abs(p1/L);     % divided by L to normalize

    p1 = p1(1:L/2+1);
    p1(2:end-1) = 2*p1(2:end-1);

    subplot(6,1,i)
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
    elseif (i==5) 
        title('Calve')
    else 
        title('Tracker')
    end

    hold on

end 


%% 4th order butterworth filter

fnyq = fs/2; %Nyquist frequency
fcuthigh = 15; %This was dicided manually
fcutlow = 250; %This was dicided manually

[b,a] = butter(4,[fcuthigh,fcutlow]/fnyq,'bandpass'); % 4th Butterworth filter - is this the best filter??

for i = 1:6
    muscle(:,i) = filtfilt(b,a,muscle(:,i));
end 


%% Taking abselute value (but not for the tracker)

rec_signal = zeros(length(muscle),5); 
 

for i=1:5
    rec_signal(:,i) = abs(muscle(:,i)); 
end 

rec_signal(:,6) = muscle(:,6);


%% Taking the standard deviation from the signal to get rid of certain frequencys

s = std(rec_signal);

for i=1:5
    rec_signal(:,i)=rec_signal(:,i)-s(i);
end

% Making the minus values zero
for i = 1:length(rec_signal)
    for j = 1:5
        if (rec_signal(i,j) < 0)
            rec_signal(i,j)=0;
        end
    end
end

%% Plotting the final signal after std and filtering
figure;
for i = 1:6
    subplot(6,1,i)

    if(i==1)
        plot(time, rec_signal(:,1),'blue');
        title('Glute')
    elseif (i==2)
        plot(time, rec_signal(:,2),'yellow');
        title('Medial Quadriceps')
    elseif (i==3)
        plot(time, rec_signal(:,3),'red');
        title('Lateral Quadriceps')
    elseif (i==4)
        plot(time, rec_signal(:,4),'magenta');
        title('Hamstring')
    elseif (i==5)
        plot(time, rec_signal(:,5),'cyan');
        title('Calve')
    else 
        plot(time,rec_signal(:,6))
        title('Tracker')
    end

    %xlim([9,14])
    xlabel('Time[s]');
    ylabel('Voltage [mV]')
end 
hold off

%% Getting exactly one loop

% Finding the locations of the peaks (loops)
[pks,locs]= findpeaks(rec_signal(:,6),time);
avgx = max(pks)-std(pks);
[pks,locs] = findpeaks(rec_signal(:,6),time,'MinPeakProminence',avgx);

%Extract the data from one loop to another

for i=1:length(locs)-1
    acc{i}=rec_signal(locs(i)*1000:locs(i+1)*1000,:);
end 


%% Plotting the final signal in a polarplot



% x1=rec_signal(1000:1600,1);
% th1 = linspace(0,2*pi,length(x1));
% x2=rec_signal(1000:1600,2);
% th2 = linspace(0,2*pi,length(x2));
% x3=rec_signal(1000:1600,3);
% th3 = linspace(0,2*pi,length(x3));
% x4=rec_signal(1000:1600,4);
% th4 = linspace(0,2*pi,length(x4));
% x5=rec_signal(1000:1600,5);
% th5 = linspace(0,2*pi,length(x5));
% 
% 
% 
% figure(7)
% subplot(2,3,1)
% h1 = polarscatter(th1,x1,'b');
% ax= ancestor(h1,'polaraxes');
% ax.ThetaZeroLocation='right';
% ax.ThetaDir="clockwise";
% title('Glute')
% subplot(2,3,2)
% h2 = polarscatter(th2,x2,'y')
% ax= ancestor(h2,'polaraxes');
% ax.ThetaZeroLocation='right';
% ax.ThetaDir="clockwise";
% title('Medial Quadriceps')
% subplot(2,3,3)
% h3 = polarscatter(th3,x3,'r')
% ax= ancestor(h3,'polaraxes');
% ax.ThetaZeroLocation='right';
% ax.ThetaDir="clockwise";
% title('Lateral Quadriceps')
% subplot(2,3,4)
% h4 = polarscatter(th4,x4,'m')
% ax= ancestor(h4,'polaraxes');
% ax.ThetaZeroLocation='right';
% ax.ThetaDir="clockwise";
% title('Hamstring')
% subplot(2,3,5)
% h5 = polarscatter(th5,x5,'c')
% ax= ancestor(h5,'polaraxes');
% ax.ThetaZeroLocation='right';
% ax.ThetaDir="clockwise";
% title('Calve')
