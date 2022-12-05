clc
clear all
close all

%% Upload the data


file = readmatrix('minute 1.txt','TrimNonNumeric',true); %here it takes the non numeric values from the file

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

%% Root mean square

for i =1:5
    rec_signal(:,i) = sqrt(movmean(rec_signal(:,i).^2, 10));
end


%% Taking the standard deviation from the signal to get rid of certain frequencys

s = std(rec_signal);

for i=1:5
    rec_signal(:,i)=rec_signal(:,i)-s(i);
end


% Making the minus values zero
%rec_signal_positive = []
for i = 1:length(rec_signal)
    for j = 1:5
        if (rec_signal(i,j) < 0)
            %rec_signal_positive(end+1,j) = rec_signal(i,j);
            rec_signal(i,j) = 0;
            
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
    lopp{i}=rec_signal(locs(i)*1000:locs(i+1)*1000,:);
end 

%% Plotting the final signal in a polarplot

% loop = lopp{17};
% x = {zeros(length(loop),5)};
% th = {zeros(length(loop),5)};
% 
% for i=1:5
%     x{i}=loop(:,i);
%     th{i}=linspace(0,2*pi,length(x{i}));
% end
% 
% 
% for i=1:5
%     newcolors = {'red','magenta','blue','green','cyan'};
%     subplot(2,3,i)
%     h = polarscatter(th{i},x{i},newcolors{i});
%     ax = ancestor(h,'polaraxes');
%     ax.ThetaZeroLocation = 'top';
%     ax.ThetaDir = "clockwise";
%      
%    
%     if i == 1
%         hold on
%         title('Glute')
%     end
%     if i == 2
%         title('Medial quadriceps')
%     end
%     if i == 3
%         
%         title('Lateral quadriceps')
%     end
%     if i == 4
%         title('Hamstring')
%     end
%     if i == 5
%         title('Calf')
%     end
% end

%% Polarplots and polarhistogram togthere
% Still have to make the code more beautiful and make it into a for loop

loop = lopp{17};
x1 = loop(:,1);
th1 =linspace(0,2*pi,length(x1));
subplot(2,3,1)
h = polarscatter(th1,x1, 'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [0,pi/2], 'BinCounts',[max(x1)], 'FaceColor','red')
title('Glute')

subplot(2,3,2)
x2 = loop(:,2);
th2 =linspace(0,2*pi,length(x2));
h = polarscatter(th2,x2, 'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [deg2rad(85),deg2rad(175)], 'BinCounts',[max(x2)], 'FaceColor','magenta')
title('Medial quadriceps')

subplot(2,3,3)
x3 = loop(:,3);
th3 =linspace(0,2*pi,length(x3));
h = polarscatter(th3,x3, 'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [deg2rad(85),deg2rad(175)], 'BinCounts',[max(x3)], 'FaceColor','blue')
title('Lateral quadriceps')

subplot(2,3,4)
x4 = loop(:,4);
th4 =linspace(0,2*pi,length(x4));
h = polarscatter(th4,x4, 'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [deg2rad(225),deg2rad(270)], 'BinCounts',[max(x4)], 'FaceColor','green')
title('Hamstring')

subplot(2,3,5)
x5 = loop(:,5);
th5 =linspace(0,2*pi,length(x5));
h = polarscatter(th5,x5, 'black');
ax = ancestor(h,'polaraxes');
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = "clockwise";
hold on 
polarhistogram('BinEdges', [deg2rad(160),deg2rad(180)], 'BinCounts',[max(x5)], 'FaceColor','cyan')
title('Calf')













