clear all
%clc
format long

% Initial time (UTC).
year = 2018;
month = 3;
day = 23;
hour = 8;
minute = 55;
second = 3;
UTC = [year month day hour minute second];
UTCdatetime = datetime(UTC);

% Initial Satellite State in ECI
R0 = [ 6985.31847092135;
       1621.22693026875;
       15.7331358497728]; % [km] 
V0 = [ -1.67429646939288;
        7.25923523079474; 
        0.261383935118773]; % [km/s]

% Spacecraft Properties
m = 2000; % [kg] mass
AreaX = 6/(1e6); % [km^2] +X/-X Face Area
AreaY = 8/(1e6); % [km^2] +Y/-Y Face Area
AreaZ = 12/(1e6); % [km^2] +Z/-Z Face Area
A_panel = 15/(1e6); % [km^2] Solar Panel Area

% Celestial constants
mu = 398600.4415; % [km^3/s^2] Earth gravitation
RE = 6378.1363; % [km] Earth radius
mu_Sun = 132712440018; % [km^3/s^2] Sun gravitation
AU = 149597870.7; % [km] 1 Astronomical Unity
mu_Moon = 4902.800066; % [km^3/s^2] Moon's gravitation
e_E = 0.081819221456; % [no units] Earth's eccentricity
om_E = 7.292115146706979e-5; % [rad/s] Earth's rot vel

% % JGM-3 Earth Gravity Model (Tapley Table D.1).
% J2 = -(-0.48416954845647e-03) * sqrt(2*2 + 1);
% J3 = -(0.95717059088800e-06) * sqrt(2*3 + 1);
% J4 = -(0.53977706835730e-06) * sqrt(2*4 + 1);
% J5 = -(0.68658987986543e-07) * sqrt(2*5 + 1);
% J6 = -(-0.14967156178604e-06) * sqrt(2*6 + 1);
% 
% % Vallado 1st Edition Table D-4. (Uses JGM-2).
% J_2 = 0.0010826269;
% J_3 = -0.0000025323;
% J_4 = -0.0000016204;

% EGM-08 (Vallado 4th Edition.
J2 = 1.08262617385222e-3;
J3 = -2.53241051856772e-6;
J4 = -1.61989759991697e-6;
J5 = -2.27753590730836e-7;
J6 = 5.40666576283813e-7;
J7 = -3.50551795713742e-7;
J8 = -2.03993125929884e-7;
J9 = -1.22127958919496e-7;
J10 = -2.44390769772693e-7;


% Drag coefficients
C_D = 1.88;
rho0 = (3.614e-13)*1e9; % [kg/km^3]
r0 = 700000/1000 + RE; % [km]
H = 88667.0/1000; % [km]

% Solar radiation pressure constants
c = 299792.458;  % speed of light [km/s]
p_sr = (1353/1e6)/c; % rought estimate of solar pressure (Vallado p. 518)
c_r = 1; % solar reflectivity (very rough estimate). 


% Load the position of the sun and mooon at Julian dates.
sun_moon_pos1 = load('sun_moon_positions_first_3_days.mat');
sun_moon_pos2 = load('sun_moon_positions_second_3_days.mat');
sun_moon_pos3 = load('sun_moon_positions_last_day.mat');
sun_moon_pos1 = cell2mat(struct2cell(sun_moon_pos1));
sun_moon_pos2 = cell2mat(struct2cell(sun_moon_pos2));
sun_moon_pos3 = cell2mat(struct2cell(sun_moon_pos3));
sun_moon_pos = [sun_moon_pos1; sun_moon_pos2; sun_moon_pos3];

% Extract time and sun/moon position data.
rsunMat = sun_moon_pos(:,2:4);
rmoonMat = sun_moon_pos(:,5:7);



% Get symbolic equations for dU/dx.
[dUdx, dUdy, dUdz] = getGradU();

% Substitute numerical values into symbolic expressions
dUdx_num = subs(dUdx, {'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'muu', 'RE', 'C_D', 'A', 'm', 'rho0', 'r0', ...
    'H', 'om_E', 'mu_Sun', 'mu_Moon', 'p_sr', 'c_r', 'A_panel', 'A_Y', 'A_Z'}, {J2, J3, J4, J5, J6, J7, J8, J9, J10, mu, RE, C_D, AreaX, m, rho0, ...
    r0, H, om_E, mu_Sun, mu_Moon, p_sr, c_r, A_panel, AreaY, AreaZ});
dUdy_num = subs(dUdy, {'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'muu', 'RE', 'C_D', 'A', 'm', 'rho0', 'r0', ...
    'H', 'om_E', 'mu_Sun', 'mu_Moon', 'p_sr', 'c_r', 'A_panel', 'A_Y', 'A_Z'}, {J2, J3, J4, J5, J6, J7, J8, J9, J10, mu, RE, C_D, AreaX, m, rho0, ...
    r0, H, om_E, mu_Sun, mu_Moon, p_sr, c_r, A_panel, AreaY, AreaZ});
dUdz_num = subs(dUdz, {'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'muu', 'RE', 'C_D', 'A', 'm', 'rho0', 'r0', ...
    'H', 'om_E', 'mu_Sun', 'mu_Moon', 'p_sr', 'c_r', 'A_panel', 'A_Y', 'A_Z'}, {J2, J3, J4, J5, J6, J7, J8, J9, J10, mu, RE, C_D, AreaX, m, rho0, ...
    r0, H, om_E, mu_Sun, mu_Moon, p_sr, c_r, A_panel, AreaY, AreaZ});

% Convert symbolic expressions to numerical functions
dUdx_func = matlabFunction(dUdx_num, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz'), sym('xMoon'), ...
    sym('yMoon'), sym('zMoon'), sym('xSun'), sym('ySun'), sym('zSun')]);
dUdy_func = matlabFunction(dUdy_num, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz'), sym('xMoon'), ...
    sym('yMoon'), sym('zMoon'), sym('xSun'), sym('ySun'), sym('zSun')]);
dUdz_func = matlabFunction(dUdz_num, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz'), sym('xMoon'), ...
    sym('yMoon'), sym('zMoon'), sym('xSun'), sym('ySun'), sym('zSun')]);



% Get the approximate sun and moon positions.
rmoon = rmoonMat(1, :);
rsun = rsunMat(1, :);

% Compare with given accelerations.
% Get the acceleration of the satellite. 
dVx = dUdx_func(R0(1), R0(2), R0(3), V0(1), V0(2), V0(3), ...
    rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));
dVy = dUdy_func(R0(1), R0(2), R0(3), V0(1), V0(2), V0(3), ...
    rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));
dVz = dUdz_func(R0(1), R0(2), R0(3), V0(1), V0(2), V0(3), ...
    rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));
acc = [dVx; dVy; dVz];


truth_acc_2BJ2J3J4 = [ -0.0075601602480519 ;
 -0.00175474460199184 ;
 -1.71462729023098e-05 ];

truth_acc_drag = [1.28565567565462e-11 ;
 -5.57684843178588e-11 ;
 -2.15959896971973e-12 
];

truth_acc_3rd = [1.03767688038019e-10 ;
 6.44383775293363e-10 ;
 2.98396598163909e-10 
];

truth_acc_srp = [ -2.68799817041286e-11 ;
 -2.9509733168165e-12;
 -4.95076071757261e-13 
];

relErr = abs((acc - truth_acc_3rd)./truth_acc_3rd)