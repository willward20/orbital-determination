%% Estimate the position of the sun and the moon at all times.  
% Initial time (UTC).
year = 2018;
month = 3;
day = 23;
hour = 8;
minute = 55;
second = 3;

sec_per_day = 24*60*60;

% Define the time span over seven days.
tspan = 518460:60:7*sec_per_day;
disp(size(tspan))

% Precompute Sun and Moon positions
init_date = datetime(year,month,day,hour,minute,second);
dates = juliandate(init_date + seconds(tspan));
rsunMat = zeros(3, length(tspan));
rmoonMat = zeros(3, length(tspan));

% Get the Sun and Moon positions [km] in 
% Earth-Centered Inertial (ECI) frame.
for i = 1:length(tspan)
    disp(i)
    [rsunMat(:,i), ~] = planetEphemeris(dates(i), 'Earth', 'Sun');
    [rmoonMat(:,i), ~] = planetEphemeris(dates(i), 'Earth', 'Moon');
end

% Combine the data into a single matrix where:
% - First column is time (Julian date)
% - Next three columns are the Sun's position
% - Last three columns are the Moon's position
time_sun_moon = [tspan', rsunMat', rmoonMat'];

% Save the data to a MAT file
save('sun_moon_positions_last_day.mat', 'time_sun_moon');