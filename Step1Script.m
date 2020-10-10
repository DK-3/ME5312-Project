clear;
data = readtable("E:\2020\ME5312\testData.csv");

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
