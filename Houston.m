clear;
data = readtable("C:\Users\nickd\OneDrive\Documents\MATLAB\ME5312\HOU_data.csv");

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