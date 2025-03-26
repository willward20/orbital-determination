clear all
clc
format long

A_t0 = load('A_t0.mat');
A_t0 = cell2mat(struct2cell(A_t0));
H_Tilde_t0 = load('H_Tilde_t0.mat');
H_Tilde_t0 = cell2mat(struct2cell(H_Tilde_t0));
STM_21600_0 = load('Phi_21600_0.mat');
STM_21600_0 = cell2mat(struct2cell(STM_21600_0));


% Celestial constants
mu = 398600.4415; % [km^3/s^2] Earth gravitation
RE = 6378.1363; % [km] Earth radius
mu_Sun = 132712440018; % [km^3/s^2] Sun gravitation
AU = 149597870.7; % [km] 1 Astronomical Unity
mu_Moon = 4902.800066; % [km^3/s^2] Moon's gravitation
e_E = 0.081819221456; % [no units] Earth's eccentricity
om_E = 7.292115146706979e-5; % [rad/s] Earth's rot vel

J2 = 0.00108248; % [no units] COPIED FROM HW 2

% Drag coefficients
C_D = 1.88;
rho0 = (3.614e-13)*1e9; % [kg/km^3]
r0 = 700000/1000 + RE; % [km]
H = 88667.0/1000; % [km]

% Spacecraft Properties
m = 2000; % [kg] mass
AreaX = 6/(1e6); % [km^2] +X/-X Face Area
AreaY = 8/(1e6); % [km^2] +Y/-Y Face Area
AreaZ = 12/(1e6); % [km^2] +Z/-Z Face Area
AreaPanel = 15/(1e6); % [km^2] Solar Panel Area

% Initial Satellite State
R0 = [6990077.798814194;
      1617465.311978378;
      22679.810569245355]/1000; % [km] ECI FRAME??
V0 = [-1675.13972506056
       7273.72441330686
       252.688512916741]/1000; % [km/s] ECI FRAME??

% Tracking Statation Coordinates (ECEF)
atoll = [-6143584; 1364250; 1033743]/1000; % [km]
diego = [1907295; 6030810; -817119]/1000; % [km]
arecibo = [2390310; -5564341; 1994578]/1000; % [km]

% Initial time (UTC).
year = 2018;
month = 2;
day = 1;
hour = 5;
minute = 0;
second = 0;
UTC = [year month day hour minute second];

% Do I need leap seconds? 
delAT = 32; % [seconds]
% Hardcoded EOP values from Vallado Ex 3-15.
dUT1 = -0.4399619; % [seconds]
% Polar motion values hardcoded from Vallado. 
pm = [-0.140682 0.333309]./3600; % [arcseconds]

% Convert initial satellite position from ECI to ECEF
R0ECEF = eci2ecef(UTC, R0, 'dAT', delAT, 'dUT1', dUT1, 'pm', pm);
V0ECEF = eci2ecef(UTC, V0, 'dAT', delAT, 'dUT1', dUT1, 'pm', pm);

%% Estimate the position of the sun and the moon at all times.  
% Define the time span
tspan = linspace(0, 21600, 60); % Choose appropriate range
% Precompute Sun and Moon positions
init_date = datetime(year,month,day,hour,minute,second);
dates = juliandate(init_date + seconds(tspan));
rsunMat = zeros(3, length(tspan));
rmoonMat = zeros(3, length(tspan));

% Get the Sun and Moon positions [km] in 
% Earth-Centered Inertial (ECI) frame.
% for i = 1:length(tspan)
%     [rsunMat(:,i), ~] = planetEphemeris(dates(i), 'Earth', 'Sun');
%     [rmoonMat(:,i), ~] = planetEphemeris(dates(i), 'Earth', 'Moon');
% end

% Calculate the sun and moon positions relative to the satellite at epoch.
% rsun0 = rsunMat(:,1);
% rmoon0 = rmoonMat(:,1);
[rsun0, ~] = planetEphemeris(dates(1), 'Earth', 'Sun');
[rmoon0, ~] = planetEphemeris(dates(1), 'Earth', 'Moon');
R0Sat2Sun = rsun0 - R0;
R0Sat2Moon = rmoon0 - R0;

% Calculate the acceleration due to sun and moon (Vallado Accelerations Due to Third Body (p. 515 book). 
aSun = mu_Sun*((R0Sat2Sun./norm(R0Sat2Sun)^3) - rsun0./norm(rsun0)^3);
aMoon = mu_Moon*((R0Sat2Moon./norm(R0Sat2Moon)^3) - rmoon0./norm(rmoon0)^3);

% Interpolation functions
% interp_rsun = @(t) interp1(tspan, rsunMat', t, 'linear')';
% interp_rmoon = @(t) interp1(tspan, rmoonMat', t, 'linear')';


%% Sanity Check: acceleration due to 2b+J2. 
% Get symbolic equations for dU/dx.
[dUdx, dUdy, dUdz] = getGradU();

% Substitute numerical values into symbolic expressions
dUdx_num = subs(dUdx, {'J2', 'muu', 'RE'}, {J2, mu, RE});
dUdy_num = subs(dUdy, {'J2', 'muu', 'RE'}, {J2, mu, RE});
dUdz_num = subs(dUdz, {'J2', 'muu', 'RE'}, {J2, mu, RE});

% Convert symbolic expressions to numerical functions
dUdx_func = matlabFunction(dUdx_num, 'Vars', [sym('x'), sym('y'), sym('z')]);
dUdy_func = matlabFunction(dUdy_num, 'Vars', [sym('x'), sym('y'), sym('z')]);
dUdz_func = matlabFunction(dUdz_num, 'Vars', [sym('x'), sym('y'), sym('z')]);

% Evaluate for R0.
acc2BJ2 = [dUdx_func(R0(1), R0(2), R0(3));
           dUdy_func(R0(1), R0(2), R0(3));
           dUdz_func(R0(1), R0(2), R0(3))];

% Ben's 2B + J2 acceleration.
ben2BJ2 = [-7.55345e-3; -1.74783e-3; -2.45705e-5];

% Check the error.
twoBJ2Err = acc2BJ2 - ben2BJ2


%% Sanity Check: acceleration due to drag.
x = R0(1);
y = R0(2);
z = R0(3);
vx = V0(1);
vy = V0(2);
vz = V0(3);
r = sqrt(x^2 + y^2 + z^2); % [km]

% Acceleration due to drag.
rhoA = rho0*exp(-(r - r0)/H);
VAvec = [vx + om_E*y; ...
         vy - om_E*x; ...
         vz];
VA = sqrt((vx + om_E*y)^2 + (vy - om_E*x)^2 + vz^2);
accDrag = -0.5*C_D*(AreaX/m)*rhoA*VA*VAvec;

% Compare against true values. Error should be on the order of 1e-16. 
dragErr = accDrag - [3.94136e-12; -1.71202e-11; -6.39572e-13]


%% Get symbolic function for A matrix.
% NOTE: Provided A mat is in km and km/s
AMat = getA();

% Substitute numerical values into A
AMat_num = subs(AMat, {'J2', 'muu', 'RE', 'C_D', 'A', 'm', 'rho0', 'r0', ...
    'H', 'om_E', 'mu_Sun', 'mu_Moon'}, {J2, mu, RE, C_D, AreaX, m, rho0, ...
    r0, H, om_E, mu_Sun, mu_Moon});

% Convert symbolic expressions to numerical functions
AMat_func = matlabFunction(AMat_num, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz'), sym('xMoon'), ...
    sym('yMoon'), sym('zMoon'), sym('xSun'), sym('ySun'), sym('zSun')]);

% Calculate A0.
AMat0 = AMat_func(R0(1), R0(2), R0(3), V0(1), V0(2), V0(3), ...
    rmoon0(1), rmoon0(2), rmoon0(3), rsun0(1), rsun0(2), rsun0(3));

% Calculate the relative error in A0. 
relDiffA = (AMat0 - A_t0)./A_t0

histogram(reshape(log10(abs(relDiffA)),49,1), 10)
xlabel("Log of the Relative Difference")
ylabel("Frequency")
title("Histogram of the Log of Relative Differences in A Matrix Components")


%% Get symbolic functions for HTilde for each station.
HtildeFunc = getHtilde();

% Substitute numerical values for Atoll station. 
HtildeAtoll = subs(HtildeFunc, {'C_D', 'xI', 'yI', 'zI', 'vxI', 'vyI', 'vzI'}, ...
    {C_D, atoll(1), atoll(2), atoll(3), 0, 0, 0});
% Substitute numerical values for Diego Garcia station. 
HtildeDiego = subs(HtildeFunc, {'C_D', 'xI', 'yI', 'zI', 'vxI', 'vyI', 'vzI'}, ...
    {C_D, diego(1), diego(2), diego(3), 0, 0, 0});
% Substitute numerical values for Arecibo station. 
HtildeArecibo = subs(HtildeFunc, {'C_D', 'xI', 'yI', 'zI', 'vxI', 'vyI', 'vzI'}, ...
    {C_D, arecibo(1), arecibo(2), arecibo(3), 0, 0, 0});

% Convert symbolic expressions to numerical functions
Ht_Atoll_func = matlabFunction(HtildeAtoll, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz')]);
Ht_Diego_func = matlabFunction(HtildeDiego, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz')]);
Ht_Arecibo_func = matlabFunction(HtildeArecibo, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz')]);


%% Calculate Htilde0 (for Tracking Station 1 - Atoll)
% Calculate Htile0 for each station. These should use ECI. 
Htilde0 = Ht_Atoll_func(R0(1), R0(2), R0(3), ...
    V0(1), V0(2), V0(3));

% Calculate error with provided Htilde0.
Htilde0 - H_Tilde_t0
relDiff = abs((Htilde0-H_Tilde_t0)./H_Tilde_t0)


%% Propagate the STM and the state forward.
% Initialize ode45.
myoptions = odeset('RelTol',3e-14,'AbsTol',1e-16);
times = 0:60:21600; % [s]

% Create an initial state vector.
XSTM0 = [R0; V0; C_D; reshape(eye(7), 7^2, 1)];

% Propagate the kinematics forward in time.
[T,Y] = ode45(@propStateAndSTM, times, XSTM0, myoptions, AMat_func, ...
    interp_rsun, interp_rmoon); 

% Extract the components of Ystar
Xt = Y(:,1:7);
STMt = reshape(Y(:,8:56)', 7, 7, 361);
if STMt(:,:,1) ~= eye(7)
    disp('ERROR: STM(1) not equal to identity!')
end

% Get the error in the STM mapping from 0 to 21600.
STMerr = STMt(:,:,end) - STM_21600_0


%% Estimate the light time. 
c = 299792458/1000; % speed of light [km/s]
tol = 1e-6; % tolerance (1 mm) [km]
del = [1; 1; 1]; % initialize delta [km]

% while del(1) > tol & del(2) > tol & del(3) > tol
%     t_corr = 
% end