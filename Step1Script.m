% Houston Data

clear;
data = readtable("C:\Users\brucevang\Documents\MATLAB\ME5312-Project-master\HOU_Data.csv");

%to calculate monthly average daily solar irradiation,
%GHI for each day of each month must be totaled
%to do this we need to add up GHI over every day of each month
month = data.Month;                 %array of months
day = data.Day;                     %array of days
GHI = data.GHI;                     %array of GHI data
lastD = 1;                          %updating value to keep track of last day
D = 1;                              %updating value for current day
num_meas = 0;                       %# of measurements so far
meas_sum = 0;                       %sum of measurements so far
dailyAvgs = [];                     %array of daily average GHI
L = length(day);                    %length of csv file

DPM = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; %days per month

%iterate down every row of file
for row_num = 1:L
    %get day of month from data
    D = day(row_num);
    %is this row from the same day as the last one
    if D ~= lastD %if not, average GHI measurements, add to array, reset measurement sums
        avg = meas_sum/num_meas;
        meas_sum = 0;
        num_meas = 0;
        dailyAvgs = [dailyAvgs, avg];
    end
    %if it is the same day still
    if GHI(row_num) > 0 %if GHI data is 0 (nighttime) skip measurement
        num_meas = num_meas+1; %keep track of # of measurements per day
        meas_sum = meas_sum+GHI(row_num); %add to a running GHI todal
    end
    lastD = D; %set last day to current day

end
avg = meas_sum/num_meas;
dailyAvgs = [dailyAvgs, avg];

Hbar = []; %%initialize array for montly daily average
    
        
DOY = 1;
lastAvg = 0;
monthTracker = 1;
dayTracker = 0;
runningAvg = 0;
for avg_num = 1:length(dailyAvgs)    
        runningAvg = runningAvg+dailyAvgs(avg_num);
        if dayTracker == DPM(monthTracker)
            monthlyDailyAvg = runningAvg/DPM(monthTracker);
            runningAvg = 0;
            dayTracker = 0;
            Hbar = [Hbar, monthlyDailyAvg];
            monthTracker = monthTracker+1;
        end
        dayTracker = dayTracker+1;
end
monthlyDailyAvg = runningAvg/DPM(monthTracker);
Hbar = [Hbar, monthlyDailyAvg];
% Convert Hbar data from W/m^2 to MJ/m^2

% Gather average n and declination to calculate Hbar_o and clearness index
n_avg = [17,47,75,105,135,162,198,228,258,288,318,344];
dec = [-20.9,-13,-2.4,9.4,18.8,23.1,21.2,13.5,2.2,-9.6,-18.9,-23]; % From tables
sunsetangle = zeros(1,length(dec));
Hbar_o = zeros(1,length(dec));
Kbar_T = zeros(1,length(dec));
N = zeros(1,length(dec));
lat = 29.76;
for i = 1:length(n_avg)
    sunsetangle(i) = acosd(-tand(dec(i))*tand(lat));
    N(i) = (2/15)*sunsetangle(i);
    Hbar_o(i) = 10^(-6)*((24*3600*1367)/pi)*(1 + 0.033*cosd((360*n_avg(i))/365))*(cosd(lat)*cosd(dec(i))*sind(sunsetangle(i)) + ((pi*sunsetangle(i))/180)*sind(lat)*sind(dec(i)));
    
end

for i = 1:length(Hbar)
    Hbar(i) = .0036*Hbar(i)*N(i); % Convert values for Hbar to MJ/m^2
    Kbar_T(i) = Hbar(i)/Hbar_o(i); % Get clearness index
end

% Calculate Hbar_d and Hbar_b
Hbar_d = zeros(1,length(Hbar));
Hbar_b = zeros(1,length(Hbar));
for i = 1:length(Hbar)
    if sunsetangle(i) <= 81.4 % Using the formula for Hbar_d/Hbar with a known sunset angle
        Hbar_d(i) = Hbar(i)*(1.391 - 3.560*Kbar_T(i) + 4.189*Kbar_T(i)^2 - 2.137*Kbar_T(i)^3);
    else 
        Hbar_d(i) = Hbar(i)*(1.311 - 3.022*Kbar_T(i) + 3.427*Kbar_T(i)^2 - 1.821*Kbar_T(i)^3);
    end
    Hbar_b(i) = Hbar(i) - Hbar_d(i); % Calculate Hbar_b
end

slope = [0:10:90];
reflectance = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2];

Hbar_T = zeros(length(slope),length(Hbar));
Hbar_bT = Hbar_T;
Hbar_dT = Hbar_T;
R_bar_b = Hbar_T;
for i = 1:length(slope)
    for j = 1:length(Hbar) % Calculate Hbar_T with the Rbar_b value in the book. 
        R_bar_b(i,j) = (cosd(lat-slope(i))*cosd(dec(j))*sind(sunsetangle(j)) + (pi/180)*sunsetangle(j)*sind(lat - slope(i))*sind(dec(j)))/(cosd(lat)*cosd(dec(j))*sind(sunsetangle(j)) + (pi/180)*sunsetangle(j)*sind(lat)*sind(dec(j)));
        Hbar_T(i,j) = Hbar_b(j)*R_bar_b(i,j) + Hbar_d(j)*((1 + cosd(slope(i)))/2) + Hbar(j)*reflectance(j)*((1 - cosd(slope(i)))/2);
        Hbar_bT(i,j) = Hbar_b(j)*R_bar_b(i,j);
        Hbar_dT(i,j) = Hbar_d(j)*((1 + cosd(slope(i)))/2);
    end
end

%% Equinoxes and solstices: 20 March, 21 June, 23 September, 21 December

March20 = find(data.Month == 3 & data.Day == 20);
IM = [];
June21 = find(data.Month == 6 & data.Day == 21);
IJ = [];
September23 = find(data.Month == 9 & data.Day == 23);
IS = [];
December21 = find(data.Month == 12 & data.Day == 21);
ID = [];

hourtracker = 0;
tracker = 0;
append = [];
for i = 1:length(March20)
    if data.GHI(March20(i)) > 0 || data.GHI(March20(i) -1) > 0 
        append = [append, March20(i)];
        if March20(i) == 3758
            IM = [IM,.0036*(data.GHI(March20(i))/2)];
        else
            tracker = GHI(March20(i)) + tracker;
            if hourtracker == 1
               hourtracker = 0;
               tracker = (tracker/2)*.0036;
               IM = [IM, tracker];
               tracker = 0;
            else
            hourtracker = hourtracker + 1;
            end
         end
    end
end

hourtracker = 0;
tracker = 0;
append = [];
for i = 1:length(June21)
    if data.GHI(June21(i)) > 0 || data.GHI(June21(i) -1) > 0 
        append = [append, June21(i)];
        if June21(i) == 8220
            IJ = [IJ,.0036*(data.GHI(June21(i))/2)];
        else
            tracker = GHI(June21(i)) + tracker;
            if hourtracker == 1
               hourtracker = 0;
               tracker = (tracker/2)*.0036;
               IJ = [IJ, tracker];
               tracker = 0;
            else
            hourtracker = hourtracker + 1;
            end
         end
    end
end

hourtracker = 0;
tracker = 0;
append = [];
for i = 1:length(September23)
    if data.GHI(September23(i)) > 0 || data.GHI(September23(i) -1) > 0 
        append = [append, September23(i)];
        if September23(i) == 12734
            IS = [IS,.0036*(data.GHI(September23(i))/2)];
        else
            tracker = GHI(September23(i)) + tracker;
            if hourtracker == 1
               hourtracker = 0;
               tracker = (tracker/2)*.0036;
               IS = [IS, tracker];
               tracker = 0;
            else
            hourtracker = hourtracker + 1;
            end
         end
    end
end

hourtracker = 0;
tracker = 0;
append = [];
for i = 1:length(December21)
    if data.GHI(December21(i)) > 0 || data.GHI(December21(i) -1) > 0 
        append = [append, December21(i)];
        if December21(i) == 17008 
            ID = [ID,.0036*(data.GHI(December21(i))/2)];
        else
            tracker = GHI(December21(i)) + tracker;
            if hourtracker == 1
               hourtracker = 0;
               tracker = (tracker/2)*.0036;
               ID = [ID, tracker];
               tracker = 0;
            else
            hourtracker = hourtracker + 1;
            end
         end
    end
end


%% Plots
months = [1,2,3,4,5,6,7,8,9,10,11,12];
figure(1);
plot(months,Hbar_T)
legend(['0' char(176)],['10' char(176)],['20' char(176)],['30' char(176)],['40' char(176)],['50' char(176)],['60' char(176)],['70' char(176)],['80' char(176)],['90' char(176)]);
legend('Location','Northwest')
set(gca,'XMinorTick','on','YMinorTick','on')
xlabel('\fontname{Times}Average Daily Monthly Radiation [MJ/m^2]','FontSize',12)
ylabel('\fontname{Times}Month','FontSize',12)


figure(2);
hold on;
hours = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19];
% Add zeros for the days that were shorter than the equinox
zero = [0];
IM = [0,IM,0];
IS = [0,IS,0];
ID = [0,0,ID,0,0];
plot(hours,ID)
plot(hours,IJ)
plot(hours,IM)
plot(hours,IS)
legend('Winter Solstice','Summer Solstice', 'March Equinox','September Equinox')
legend('Location','Northwest')
set(gca,'XMinorTick','on','YMinorTick','on')
xlabel('\fontname{Times}Average Hourly Radiation [MJ/m^2]','FontSize',12)
ylabel('\fontname{Times}Hour','FontSize',12)

%% Step 2 Question 4 - Finding Qu
%Finding Ibar
a = .409 + .5016*sind(sunsetangle - 60);
b = .6609 - .4767*sind(sunsetangle - 60);
hour_angle = -172.5:15:172.5; %from 12:30 AM to 11:30 PM
r_t = zeros(length(hour_angle), length(sunsetangle));
for i = 1:length(hour_angle)
    for j = 1:length(sunsetangle)
        r_t(i,j) = (pi/24)*(a(j)+b(j)*cosd(hour_angle(i))) * (cosd(hour_angle(i)) - cosd(sunsetangle(j)))...
            / (sind(sunsetangle(j)) - (pi*sunsetangle(j)*cosd(sunsetangle(j)))/180);
        if r_t(i,j) <= 0 %setting all negative values to 0
            r_t(i,j) = 0;
        else
            r_t(i,j) = r_t(i,j);
        end
    end
end
for i = 1:24
    for j = 1:12
        Ibar(i,j) = r_t(i,j)*Hbar(j); %Ibar for the midpoint of each hour in a day for each month
    end
end

%Finding kbar_T
omega1 = -180:15:165;
omega2 = -165:15:180;
for i = 1:length(omega1)
    for j = 1:length(n_avg)
        Ibar_o(i,j) = 10^(-6)*((12*3600*1367)/pi)*(1 + 0.033*cosd((360*n_avg(j))/365))*(cosd(lat)*cosd(dec(j))*(sind(omega2(i))-sind(omega1(i)))+...
            (pi*(omega2(i)-omega1(i))*sind(lat)*sind(dec(j))/180));
        if Ibar_o(i,j) <= 0 %setting all negative values to 0
            Ibar_o(i,j) = 0;
        else
            Ibar_o(i,j) = Ibar_o(i,j);
        end
    end
end
for i = 1:24
    for j = 1:12
        if Ibar_o(i,j) == 0
            kbar_T(i,j) = 0;
        else
            kbar_T(i,j) = Ibar(i,j)/Ibar_o(i,j);
        end
    end
end

%erbs correlation
for i = 1:24
    for j = 1:12
        if kbar_T(i,j) <= .22
            erbs(i,j) = 1 - .09*kbar_T(i,j);
            if erbs(i,j) == 1
                erbs(i,j) = 0;
            end
        elseif (kbar_T(i,j) <= .8) && (kbar_T(i,j) > .22)
            erbs(i,j) = .9511 - .1604*kbar_T(i,j) + 4.388*kbar_T(i,j)^2 - 16.638*kbar_T(i,j)^3 + 12.336*kbar_T(i,j)^4;
        elseif kbar_T(i,j) > .8
            erbs(i,j) = .165;
        end
    end
end
for i = 1:24
    for j = 1:12
        Ibar_d(i,j) = Ibar(i,j)*erbs(i,j);
        Ibar_b(i,j) = Ibar(i,j) - Ibar_d(i,j);
    end
end

%To further simplify, set slope = latitude -->
beta = lat; %Note "beta" is used instead of "slope" since the latter has been used earlier in this code
for i = 1:24
    for j = 1:12
        Rb(i,j) = cosd(lat-beta)*cosd(dec(j))*cosd(hour_angle(i))/(cosd(lat)*cosd(dec(j))*cosd(hour_angle(i))+sind(lat)*sind(dec(j))); 
        %eqn 1.8.2
        %i --> midpoint of hour, j given month
    end
end

%FR(taualpha)& FR_UL --> SRCC Data Sheet
FRtaualpha = .589;
FR_UL = -3.564;

%Calculate Ibar_T * Kbar_taualpha:
b_o = -.1881; %attained via linear regression of 1/cos(theta)-1 and K_taualpha (done via Python)
theta_e_g = 90 - .5788*beta + .002693*beta^2; %eqn 5.4.1
theta_e_d = 59.7 - .1388*beta + .001497*beta^2; %eqn 5.4.2
for i = 1:24
    for j = 1:12
        theta_e_b(i,j) = acosd(cosd(lat-beta)*cosd(dec(j))*cosd(hour_angle(i))); %eqn 1.6.7, note that the second term on the RHS of the equation cancels out since sin(lat-beta)=0
    end
end

%Find K_taualpha for these angles
K_taualpha_g = 1-.1881*(1/cosd(theta_e_g) - 1);
K_taualpha_d = 1-.1881*(1/cosd(theta_e_d) - 1);
for i = 1:24
    for j = 1:12
        K_taualpha_b(i,j) = 1-.1881*(1/cosd(theta_e_b(i,j)) - 1);
    end
end

for i = 1:24
    for j = 1:12 
        Ibar_T_and_Kbar_taualpha(i,j) = Ibar_b(i,j)*Rb(i,j)*K_taualpha_b(i,j) + Ibar_d(i,j)*K_taualpha_d*.5*(1+cosd(beta)) + Ibar(i,j)*.2*K_taualpha_d*.5*(1-cosd(beta));
        %Reliable values from 8:30 AM to 5:30 PM (or alternatively from 8 AM to 6 PM)
    end
end

%Calculate Qu
Ac = .944; %SRCC --> Gross Collector Area
%Need inlet and ambient temperature in order to calculate Qu
%assume Ti = 30, 40, 50, degrees F --> 4.44, 10, 15.56 degrees C
%Assume Ta varies monthly via TMY data. Do it for the avg day of each month (17th Jan, 16th Feb, 16th Mar, 15 Apr, 15 May, 11 June, 17 July,
%16 Aug, 15 Sept, 15 Oct, 14 Nov, 10 Dec)
%Ta will be for i & j indices (midpoint hour & avg day of month)
Ti = [4.44, 10, 15.56];
data2 = readtable("C:\Users\brucevang\Documents\MATLAB\ME5312-Project-master\HOU_Data_TMY.csv");
month2 = data2.Month;
Temp2 = data2.Temperature;
Month_Temp = [month2,Temp2];
%For downloaded TMY data with first two rows deleted (A1 = "Source", A2 = "NSRDB")
%The average day for each month will correspond to the following rows:

for i = 1:12 %Month (On Avg Day)
    for j = 1:24 %Midpoint Hour
        if i == 1
            Ta1(i,j) = Month_Temp(384+j,2);
        elseif i == 2 
            Ta1(i,j) = Month_Temp(1104+j,2);
        elseif i == 3
            Ta1(i,j) = Month_Temp(1776+j,2);
        elseif i == 4
            Ta1(i,j) = Month_Temp(2496+j,2);
        elseif i == 5
            Ta1(i,j) = Month_Temp(3216+j,2);
        elseif i == 6
            Ta1(i,j) = Month_Temp(3864+j,2);
        elseif i == 7
            Ta1(i,j) = Month_Temp(4728+j,2);
        elseif i == 8
            Ta1(i,j) = Month_Temp(5448+j,2);
        elseif i == 9
            Ta1(i,j) = Month_Temp(6168+j,2);
        elseif i == 10
            Ta1(i,j) = Month_Temp(6888+j,2);
        elseif i == 11
            Ta1(i,j) = Month_Temp(7609+j,2);
        elseif i == 12
            Ta1(i,j) = Month_Temp(8232+j,2);
        end
    end
end

%Transpose to rearrange into Midpoint Hour rows, Average Day in Month columns
Ta2 = transpose(Ta1);

%Qu Calculation
for i = 1:24 %Midpoint Hour
    for j = 1:12 %Month (On Avg Day)
        for k = 1:3 %Varying Ti
            Qu(i,j,k) = Ac*(Ibar_T_and_Kbar_taualpha(i,j)*FRtaualpha - FR_UL*3600*(Ti(k) - Ta2(i,j))/(1e+6));
            if Qu(i,j,k) < 0 %Set negative Qu values equal to 0 since it's equivalent to having no useful energy
                Qu(i,j,k) = 0;
            end
            %In MJ
            %Note conversion factors 3600 and 1e+6 are necessary
        end
    end
end

%Qu values are negative for i = [1,7] and i = [18,24] (i.e. from 12 AM to 8 AM and 6 PM to 12 PM)
%This means that there is no useful energy for these hours
%Most reliable values from Qu is from i = [8, 17]
%Therefore, graph values from 8 AM to 7 PM
%3 Graphs, 1 for each different Ti
%X-axis Month, Y-Axis Qu, 10 lines for each hour interval (8-9 AM, 9-10 AM, ..., 5-6 PM)

figure(3);
hold on;
plot(months, Qu(8,1:12,1))
plot(months, Qu(9,1:12,1))
plot(months, Qu(10,1:12,1))
plot(months, Qu(11,1:12,1))
plot(months, Qu(12,1:12,1))
plot(months, Qu(13,1:12,1))
plot(months, Qu(14,1:12,1))
plot(months, Qu(15,1:12,1))
plot(months, Qu(16,1:12,1))
plot(months, Qu(17,1:12,1))
legend('8-9 AM','9-10 AM', '10-11 AM','11 AM-12 PM','12-1 PM','1-2 PM','2-3 PM','3-4 PM','4-5 PM','5-6 PM')
set(gca,'YMinorTick','on')
xlim([1,12])
ylim([0,2])
xlabel('\fontname{Times}Month','FontSize',12)
ylabel('\fontname{Times}Useful Energy, Qu [MJ]','FontSize',12)
title('Useful Energy Throughout the Year for Ti = 30 F = 4.44 C')

figure(4);
hold on;
plot(months, Qu(8,1:12,2))
plot(months, Qu(9,1:12,2))
plot(months, Qu(10,1:12,2))
plot(months, Qu(11,1:12,2))
plot(months, Qu(12,1:12,2))
plot(months, Qu(13,1:12,2))
plot(months, Qu(14,1:12,2))
plot(months, Qu(15,1:12,2))
plot(months, Qu(16,1:12,2))
plot(months, Qu(17,1:12,2))
legend('8-9 AM','9-10 AM', '10-11 AM','11 AM-12 PM','12-1 PM','1-2 PM','2-3 PM','3-4 PM','4-5 PM','5-6 PM')
set(gca,'YMinorTick','on')
xlim([1,12])
ylim([0,2])
xlabel('\fontname{Times}Month','FontSize',12)
ylabel('\fontname{Times}Useful Energy, Qu [MJ]','FontSize',12)
title('Useful Energy Throughout the Year for Ti = 40 F = 10 C')

figure(5);
hold on;
plot(months, Qu(8,1:12,3))
plot(months, Qu(9,1:12,3))
plot(months, Qu(10,1:12,3))
plot(months, Qu(11,1:12,3))
plot(months, Qu(12,1:12,3))
plot(months, Qu(13,1:12,3))
plot(months, Qu(14,1:12,3))
plot(months, Qu(15,1:12,3))
plot(months, Qu(16,1:12,3))
plot(months, Qu(17,1:12,3))
legend('8-9 AM','9-10 AM', '10-11 AM','11 AM-12 PM','12-1 PM','1-2 PM','2-3 PM','3-4 PM','4-5 PM','5-6 PM')
set(gca,'YMinorTick','on')
xlim([1,12])
ylim([0,2])
xlabel('\fontname{Times}Month','FontSize',12)
ylabel('\fontname{Times}Useful Energy, Qu [MJ]','FontSize',12)
title('Useful Energy Throughout the Year for Ti = 50 F = 15.56 C')
