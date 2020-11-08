clear;
data = readtable("E:\2020\ME5312\debugData.csv");

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

DPM = [31, 30, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; %days per month

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
dayTracker = 1;
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
lat = 44.99;
for i = 1:length(n_avg)
    sunsetangle(i) = acosd(-tand(dec(i))*tand(lat));
    N(i) = (2/15)*sunsetangle(i);
    Hbar_o(i) = 10^(-6)*((24*3600*1367)/pi)*(1 + 0.033*cosd((360*n_avg(i))/365))*(cosd(lat)*cosd(dec(i))*sind(sunsetangle(i)) + ((pi*sunsetangle(i))/180)*sind(lat)*sind(dec(i)));
    
end

for i = 1:length(Hbar)
    Hbar(i) = .0036*Hbar(i)*N(i);
    Kbar_T(i) = Hbar(i)/Hbar_o(i);
end

% Calculate Hbar_d and Hbar_b
Hbar_d = zeros(1,length(Hbar));
Hbar_b = zeros(1,length(Hbar));
for i = 1:length(Hbar)
    if sunsetangle(i) <= 81.4
        Hbar_d(i) = Hbar(i)*(1.391 - 3.560*Kbar_T(i) + 4.189*Kbar_T(i)^2 - 2.137*Kbar_T(i)^3);
    else 
        Hbar_d(i) = Hbar(i)*(1.311 - 3.022*Kbar_T(i) + 3.427*Kbar_T(i)^2 - 1.821*Kbar_T(i)^3);
    end
    Hbar_b(i) = Hbar(i) - Hbar_d(i);
end

slope = [0,10,20,30,40,45];
reflectance = [0.7,0.7,0.5,0.4,0.2,0.2,0.2,0.2,0.3,0.4,0.5,0.7];

Hbar_T = zeros(length(slope),length(Hbar));
Hbar_bT = Hbar_T;
Hbar_dT = Hbar_T;
R_bar_b = Hbar_T;
for i = 1:length(slope)
    for j = 1:length(Hbar) % Calculate Hbar_T with the Rbar_b value in the book. 
        if acosd(-tand(lat - slope(i))*tand(dec(j))) < sunsetangle(j)
            sunsetprime = acosd(-tand(lat - slope(i))*tand(dec(j)));
        else
            sunsetprime = sunsetangle(j);
        end
        R_bar_b(i,j) = (cosd(lat-slope(i))*cosd(dec(j))*sind(sunsetprime) + (pi/180)*sunsetprime*sind(lat - slope(i))*sind(dec(j)))/(cosd(lat)*cosd(dec(j))*sind(sunsetangle(j)) + (pi/180)*sunsetangle(j)*sind(lat)*sind(dec(j)));
        Hbar_T(i,j) = Hbar_b(j)*R_bar_b(i,j) + Hbar_d(j)*((1 + cosd(slope(i)))/2) + Hbar(j)*reflectance(j)*((1 - cosd(slope(i)))/2);
        Hbar_bT(i,j) = Hbar_b(j)*R_bar_b(i,j);
        Hbar_dT(i,j) = Hbar_d(j)*((1 + cosd(slope(i)))/2);
    end
end

%% Equinoxes and solstices: 20 March, 21 June, 23 September, 21 December

HMarch20 = dailyAvgs(79);
HJune21 = dailyAvgs(172);
HSeptember23 = dailyAvgs(266);
HDecember21 = dailyAvgs(355);
n_day = [79,172,266,355];
for i = 1:length(n_day)
    dec(i) = 23.45*sind(360*((284+ n_day(i))/365));
    sunset(i) = acosd(-tand(dec(i))*tand(lat));
    leng(i) = (2/15)*sunset(i);
end
marchhour = [-89.193,-82.5,-67.5,-52.5,-37.5,-22.5,-7.5,7.5,22.5,37.5,52.5,67.5,82.5,89.193];
junehour = [-115.7,-112.5,-97.5,-82.5,-67.5,-52.5,-37.5,-22.5,-7.5,7.5,22.5,37.5,52.5,67.5,82.5,97.5,112.5,115.7];
septemberhour = [-88.99,-82.5,-67.5,-52.5,-37.5,-22.5,-7.5,7.5,22.5,37.5,52.5,67.5,82.5,88.99];
decemberhour = [-64.3,-52.5,-37.5,-22.5,-7.5,7.5,22.5,37.5,52.5,64.3];

for i = 1:length(marchhour)
    r_t_march(i) = (pi/24)*((0.409 + 0.5016*sind(sunset(1) - 60)) + ((0.6609 - 0.4767*sind(sunset(1) - 60))*cosd(marchhour(i))))*((cosd(marchhour(i)) - cosd(sunset(1)))/(sind(sunset(1)) - ((pi*sunset(1))/180)*cosd(sunset(1))));
    I_march(i) = (HMarch20*.0036*leng(1))*r_t_march(i);
    % Take K_bar_T values from each month and say assume they are same for
    % specific day
    H_d_march =  Kbar_T(3)*HMarch20*.0036*leng(1);
    I_d_march(i) = ((pi/24)*((cosd(marchhour(i)) - cosd(sunset(1)))/(sind(sunset(1)) - ((pi*sunset(1))/180)*cosd(sunset(1)))))*H_d_march;
    I_b_march(i) = I_march(i) - I_d_march(i);
end

for i = 1:length(septemberhour)
    r_t_sept(i) = (pi/24)*((0.409 + 0.5016*sind(sunset(3) - 60)) + ((0.6609 - 0.4767*sind(sunset(3) - 60))*cosd(septemberhour(i))))*((cosd(septemberhour(i)) - cosd(sunset(3)))/(sind(sunset(3)) - ((pi*sunset(3))/180)*cosd(sunset(3))));
    I_sept(i) = (HSeptember23*.0036*leng(3))*r_t_sept(i);
    H_d_september =  Kbar_T(9)*HSeptember23*.0036*leng(3);
    I_d_sept(i) = ((pi/24)*((cosd(septemberhour(i)) - cosd(sunset(3)))/(sind(sunset(3)) - ((pi*sunset(3))/180)*cosd(sunset(3)))))*H_d_september;
    I_b_sept(i) = I_sept(i) - I_d_sept(i);
end
for i = 1:length(junehour)
    r_t_june(i) = (pi/24)*((0.409 + 0.5016*sind(sunset(2) - 60)) + ((0.6609 - 0.4767*sind(sunset(2) - 60))*cosd(junehour(i))))*((cosd(junehour(i)) - cosd(sunset(2)))/(sind(sunset(2)) - ((pi*sunset(3))/180)*cosd(sunset(2))));
    I_june(i) = (HJune21*.0036*leng(2))*r_t_june(i);
    H_d_june =  Kbar_T(6)*HJune21*.0036*leng(2);
    I_d_june(i) = ((pi/24)*((cosd(junehour(i)) - cosd(sunset(2)))/(sind(sunset(2)) - ((pi*sunset(2))/180)*cosd(sunset(2)))))*H_d_june;
    I_b_june(i) = I_june(i) - I_d_june(i);
end
for i = 1:length(decemberhour)
    r_t_dec(i) = (pi/24)*((0.409 + 0.5016*sind(sunset(4) - 60)) + ((0.6609 - 0.4767*sind(sunset(4) - 60))*cosd(decemberhour(i))))*((cosd(decemberhour(i)) - cosd(sunset(4)))/(sind(sunset(4)) - ((pi*sunset(4))/180)*cosd(sunset(4))));
    I_dec(i) = (HDecember21*.0036*leng(4))*r_t_dec(i);
    H_d_dec =  Kbar_T(12)*HDecember21*.0036*leng(4);
    I_d_dec(i) = ((pi/24)*((cosd(decemberhour(i)) - cosd(sunset(4)))/(sind(sunset(4)) - ((pi*sunset(4))/180)*cosd(sunset(4)))))*H_d_dec;
    I_b_dec(i) = I_dec(i) - I_d_dec(i);
end

      
%% Plots
months = [1,2,3,4,5,6,7,8,9,10,11,12];
figure(1);
plot(months,Hbar_T)
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('\fontname{Times}Average Daily Monthly Radiation [MJ/m^2]','FontSize',12)
xlabel('\fontname{Times}Month','FontSize',12);
xlim([1 12])
title('\fontname{Times}Average Daily Monthly Radiation in Minneapolis')
legend(['0' char(176)],['10' char(176)],['20' char(176)],['30' char(176)],['40' char(176)],['45' char(176)]);

figure(2);
hold on;
plot(marchhour,I_march)
plot(marchhour,I_b_march)
plot(marchhour,I_d_march)
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('\fontname{Times}Hourly Radiation [MJ/m^2]','FontSize',12)
xlabel('\fontname{Times}Hour angle, \omega [deg]','FontSize',12);
title('\fontname{Times}March Equinox in Minneapolis')
xlim([-sunset(1) sunset(1)])
legend('I','I_b','I_d')
hold off;

figure(3);
hold on;
plot(junehour,I_june)
plot(junehour,I_b_june)
plot(junehour,I_d_june)
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('\fontname{Times}Hourly Radiation [MJ/m^2]','FontSize',12)
xlabel('\fontname{Times}Hour angle, \omega [deg]','FontSize',12);
title('\fontname{Times}Summer Solstice in Minneapolis')
xlim([-sunset(2) sunset(2)])
legend('I','I_b','I_d')
hold off;

figure(4)
hold on;
plot(septemberhour,I_sept)
plot(septemberhour,I_b_sept)
plot(septemberhour,I_d_sept)
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('\fontname{Times}Hourly Radiation [MJ/m^2]','FontSize',12)
xlabel('\fontname{Times}Hour angle, \omega [deg]','FontSize',12);
title('\fontname{Times}September Equinox in Minneapolis')
xlim([-sunset(3) sunset(3)])
legend('I','I_b','I_d')
hold off;

figure(5)
hold on;
plot(decemberhour,I_dec)
plot(decemberhour,I_b_dec)
plot(decemberhour,I_d_dec)
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('\fontname{Times}Hourly Radiation [MJ/m^2]','FontSize',12)
xlabel('\fontname{Times}Hour angle, \omega [deg]','FontSize',12);
title('\fontname{Times}Winter Solstice in Minneapolis')
xlim([-sunset(4) sunset(4)])
ylim([0 1.4])
legend('I','I_b','I_d')
hold off;
