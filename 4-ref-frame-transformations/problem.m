%% IAU-76/FK5 reduction of position from ECEF (ITRF) to ECI (GCRF)
% Define position and velocity in ECEF (ITRF).
rECEF = [-1033.4793830; 7901.2952754; 6380.3565958]; % [km]
vECEF = [-3.225636520; -2.872451450; 5.531924446]; % [km/s]

% Define position and velocity in ECI (GCRF).
% Earth Centered Inertial (ECI). 
rGCRF = [5102.508958; 6123.011401; 6378.136928]; % [km]
vGCRF = [-4.74322016; 0.79053650; 5.53375528]; % [km/s]

% Final time (UTC).
year = 2004;
month = 4;
day = 6;
hour = 7;
minute = 51;
second = 28.386009;

% Create a datetime variable.
dateTimeUTC = datetime(year,month,day,hour,minute,second);

% DUT1 value for April 6, 2004 reported in "finalls.all.csv" by the IERS. 
% DUT1 = -0.4399498; 

% Convert UTC to UT1.
dateTimeUT1 = dateTimeUTC;
dateTimeUT1.Second = dateTimeUTC + DUT1; 

% NOTE: you need to do more conversions here (see Vallado ex. 3-15 p, 230). 

% Convert Calendart Date + UT1 time into Julian date using this online
% convertor to start with https://aa.usno.navy.mil/calculated/juliandate?ID=AA&date=2004-04-06&era=AD&time=07%3A51%3A27.900&submit=Get+Date 
JDUT1 = 2453101.827406; % Julian Date in UT1 time.
% Here's an equation from Vallado Algorithm 14
% JDUT1 = 367*year - (7*(year + ((month+9)/12))/4) ...
%         + (275*month/9) + day + 1721013.5 ...
%         + (((UT1second/60) + minute)/60 + hour)/24;

% Convert the Julian date (JDUT1) to Greenwich Mean Sidereal Time (GMST)
% using a function provided by the TA. 
GMST = JD2GMST(JDUT1);

% Define the nutation in longitude.
del_nu = 0.0;



%%
clear all 
clc

% Define position and velocity in ECEF (ITRF).
rECEF = [-1033.4793830; 7901.2952754; 6380.3565958]; % [km]
vECEF = [-3.225636520; -2.872451450; 5.531924446]; % [km/s]

% Define position and velocity in ECI (GCRF). 
rGCRF = [5102.508958; 6123.011401; 6378.136928]; % [km]
vGCRF = [-4.74322016; 0.79053650; 5.53375528]; % [km/s]

% Load the nutation date (1980) and extract values.
nu_data = readmatrix('nut80.dat');
a_n1 = nu_data(:,1);
a_n2 = nu_data(:,2);
a_n3 = nu_data(:,3);
a_n4 = nu_data(:,4);
a_n5 = nu_data(:,5);
A = nu_data(:,6); % [0.0001 arcseconds]
B = nu_data(:,7); % [0.0001 arcseconds]
C = nu_data(:,8); % [0.0001 arcseconds]
D = nu_data(:,9); % [0.0001 arcseconds]

% Final time (UTC).
year = 2004;
month = 4;
day = 6;
hour = 7;
minute = 51;
second = 28.386009;

% Create a datetime variable.
dateTimeUTC = datetime(year,month,day,hour,minute,second);

% Hardcoded EOP values from Vallado Ex 3-15.
dUT1 = -0.4399619; % [seconds]
xp = -0.140682; % [arcseconds]
yp = 0.333309; % [arcseconds]
LOD = 0.0015563; % [seconds]
% Leap seconds.
delAT = 32; % [seconds]
% Correction for nutation and obliquity. 
deps1980 = -0.003875; % [arcseconds]
dpsi1980 = -0.052195; % [arcseconds]

% Convert UTC to UT1.
dateTimeUT1 =  dateTimeUTC + seconds(dUT1);

% Get the Julian date for UT1. 
JDUT1 = juliandate(dateTimeUT1);

% Get fractional century for JDUT1.
T_UT1 = (JDUT1 - 2451545.0) / 36525.0; % [unitless]

% Convert UTC to TAI.
dateTimeTAI = dateTimeUTC + seconds(delAT);

% Get TT.
dateTimeTT = dateTimeTAI + seconds(32.184);

% Get Julian date from TT.
JDTT = juliandate(dateTimeTT);

% Check for errors in JDTT.
disp('Check for errors in JDTT ==========================================')
JDTTerr = 2453101.828154745 - JDTT

%% Calculate T_TT (equation from Vallado example 3-15). 
% T_TT is fractional century. 
T_TT = (JDTT - 2451545.0)/36525; % [unitless]

% Check for T_TT errors against Vallado.
disp('Check for errors in T_TT ==========================================')
T_TTErr = 0.0426236319 - T_TT

%% Convert ECEF to PEF
% Convert xp and yp to radians.
xp_rad = xp/3600*pi/180; % [rad]
yp_rad = yp/3600*pi/180; % [rad]

% Use the small angle approximation to define the transformation 
% matrix from ITRF to PEF (Eq 3-78 from Vallado).
W = [1,            0,   -xp_rad
     0,            1,    yp_rad
     xp_rad  -yp_rad,    1];
% For reference, the full equation (Eq 3-77)
% W = [cos(xp),          0,       -sin(xp)
%      sin(xp)*sin(yp), cos(yp), cos(xp)*sin(yp)
%      sin(xp)*cos(yp), -sin(yp), cos(xp)*cos(yp)];

rPEF = orthodcm(W)*rECEF;
vPEF = orthodcm(W)*vECEF;

% Check for errors.
disp("Check for errors in rPEF and VPEF =================================")
rPEFerr = [-1033.4750313; 7901.3055856; 6380.3445328] - rPEF
vPEFerr = [-3.225632747; -2.872442511; + 5.531931288] - vPEF


%% Get Fundamental Nutation Arguments.
% Convert 360 degrees to arcseconds
r = 360; % [degrees]

% Calculate the fundamental nutation arguments. 
M_moon  = 134.96298139 + (1325*r + 198.8673981)*T_TT + 0.0086972*T_TT^2 + (1.78e-5)*T_TT^3;
M_sun  = 357.52772333 + (99*r + 359.0503400)*T_TT   - 0.0001603*T_TT^2 + (3.3e-8)*T_TT^3;
u_Mmoon =  93.27191028 + (1342*r + 82.0175381)*T_TT  - 0.0036825*T_TT^2 - (3.1e-7)*T_TT^3;
D_sun  = 297.85036306 + (1236*r + 307.1114800)*T_TT - 0.0019142*T_TT^2 + (5.3e-6)*T_TT^3;
Om_moon = 125.04452222 - (5*r + 134.1362608)*T_TT    + 0.0020708*T_TT^2 + (2.2e-6)*T_TT^3;

% Convert these angles back to values between 0 and 360 degrees.
M_moon = mod(M_moon, 360); % [degrees]
M_sun = mod(M_sun, 360); % [degrees]
u_Mmoon = mod(u_Mmoon, 360); % [degrees]
D_sun = mod(D_sun, 360); % [degrees]
Om_moon = mod(Om_moon, 360); % [degrees]

% Check for errors in the angles of nutation.
disp("Check for errors in M_moon, M_circ, u_Mmoon, D_circ, Om_moon ======")
M_moonErr = 314.9118590 - M_moon
M_circErr = 91.9379931 - M_sun
u_MmoonErr = 169.0968272 - u_Mmoon
D_circErr = 196.7518116 - D_sun
Om_moon_Err = 42.6046140 - Om_moon


%% Get the Change in Nutation Arguments.
% Calculate the nutation in longitude (1980) and the 
% nutation in obliquity (1980) using nut80.dat. 
% Note that we divide by 10000 to get A, B, C, and D
% into arcseconds. 
delPsi1980 = 0.0;
deleps1980 = 0.0;
for ii = 1:106
    a_p_i =   a_n1(ii)*M_moon ...
           + a_n2(ii)*M_sun  ...
           + a_n3(ii)*u_Mmoon ...
           + a_n4(ii)*D_sun ...
           + a_n5(ii)*Om_moon; % [degrees]
    delPsi1980 = delPsi1980 + (A(ii) + B(ii)*T_TT)*sin(a_p_i*pi/180)/10000;
    deleps1980 = deleps1980 + (C(ii) + D(ii)*T_TT)*cos(a_p_i*pi/180)/10000;
end

% Convert the nutation from arcseconds to degrees. 
delPsi1980 = delPsi1980/3600; % [degrees]
deleps1980 = deleps1980/3600; % [degrees]

% Check for errors in delNuLong and delNuObli.
disp('Check for errors in delNuLong and delNuObli =======================')
delPsiErr = -0.0034108 - delPsi1980
delepsErr = 0.0020316 - deleps1980


%% Get the Mean Obliquity. 
% Add the nutation and obliquity corrections from 
% Eq 3-84 in Vallado and convert to degrees. 
delPsi1980_cor = delPsi1980 + dpsi1980/3600; % [degrees]
deleps1980 = deleps1980 + deps1980/3600; % [degrees]

% Compute mean obliquity of the ecliptic (epsilon_bar_1980) in arcseconds
epsBar1980 = (84381.448 - 46.8150*T_TT - 0.00059*T_TT^2 + 0.001813*T_TT^3);
% Convert to degrees. 
epsBar1980deg = epsBar1980 / 3600;

% Check for errors in mean obliquity.
disp('Check for errors in epsBar1980 ====================================')
epsBarErr = 23.4387368 - epsBar1980deg


%% Calculate eps in degrees.
eps = epsBar1980deg + deleps1980; % [degrees]

% Check for errors in eps.
disp('Check for errors in eps (w/o corrections) =========================')
epsErr = 23.4407685 - (eps - deps1980/3600)


%% Get GMST1982.
% Calculate GMST1982 (https://astrogreg.com/risesetalgorithm.html)
thetaGMST = mod(280.46061837 + 360.98564736629 * (JDUT1 - 2451545.0) ...
                + 0.000387933 * T_UT1^2 - (T_UT1^3) / 38710000.0, 360); % [degrees]

% Check for errors in GMST1982.
disp('Check for errors in GMST1982 ======================================')
tGMSTerr = 312.8098943 - thetaGMST


%% Calcualte thetaGAST1982.
% Convert from degrees to radians.
epsBar1980rad = epsBar1980deg*pi/180; % [radians]
Om_moon_rad = Om_moon*pi/180; % [radians]

% Get the equation of equinoxes (Eq 3-79 from Vallado).
equniox1982 = delPsi1980*cos(epsBar1980rad) + 0.00264*sin(Om_moon_rad)/3600 ...
              + 0.000063*sin(2*Om_moon_rad)/3600;
thetaGAST = thetaGMST + equniox1982; % [degrees]

% Check for errors in thetaGast.
disp('Check for errors in thetaGast =====================================')
thetaGASTerr = 312.8067654 - thetaGAST % [degrees]


%% Get the rotation rate of the Earth. 
omega = [0; 0; (7.29211514670697910e-5)*(1 - LOD/86400)];
% Get the correction (omega cross rPEF).
omXrPEF = cross(omega, rPEF); % [km/s]

% Check errors in omega cross rPEF.
disp('Check for errors in omega cross rPEF ==============================')
crossErr = [-0.57617229; -0.07536219; 0.000]' - omXrPEF'

%% Convert PEF to TOD 
% Convert thetaGAST to radians. 
tGASTrad = thetaGAST*pi/180;
% Calculate rotation matrix for thetaGAST about z axis. 
% Note that I do NOT do this rotation about -thetaGAST. 
Rot3GAST = [cos(tGASTrad),  -sin(tGASTrad),  0
            sin(tGASTrad),   cos(tGASTrad),  0
            0             ,   0             ,  1];

% Convert PEF to TOD.
rTOD = orthodcm(Rot3GAST)*rPEF; % [km]
vTOD = orthodcm(Rot3GAST)*(vPEF + omXrPEF); % [km/s]

% Check for errors in rTOD and vTOD.
disp('Check for errors in rTOD and vTOD =================================')
rTODerr = [5094.5147804; 6127.3664612; 6380.3445328] - rTOD
vTODerr = [-4.746088567; 0.786077222; 5.531931288] - vTOD


%% TODO: Convert TOD to MOD. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert delPsi and eps to radians. 
dNL = delPsi1980_cor*pi/180; % [rad]
epsrad = eps*pi/180; % [rad]

% For siplicity, change the notation.
e = epsrad; % [rad]
eB = epsBar1980rad; % [rad]

% Equation after Eq. 3-72 from Vallado.
% THIS STILL DOES NOT WORK!!!!!!
R131 = [ cos(dNL),              sin(dNL)*cos(e) ,                                    sin(e)*sin(dNL)
        -sin(dNL)*cos(eB),  cos(e)*cos(dNL)*cos(eB) + sin(e)*sin(eB),  sin(e)*cos(dNL)*cos(eB) - sin(eB)*cos(e)
        -sin(dNL)*sin(eB),  sin(eB)*cos(e)*cos(dNL) - sin(e)*cos(eB),  sin(e)*sin(eB)*cos(dNL) + cos(e)*cos(eB)];

rMOD = orthodcm(R131)*rTOD;
vMOD = orthodcm(R131)*vTOD;

% COPIED FROM VALLADO
rmod = [5094.0283745; 6127.8708164; 6380.2485164];
vmod = [-4.746263052; 0.786014045; 5.531790562];

% Check for errors in rMOD and vMOD.
disp("Check for errors in rMOD and vMOD =================================")
rMODerr = rmod - rMOD
vMODerr = vmod - vMOD


%% Determine roation angles for precession from Eq. 3-88. 
% Note we convert from arcseconds to degrees by dividing by 3600. 
eta = 2306.2181*T_TT/3600 + 0.30188*T_TT^2/3600 + 0.017998*T_TT^3/3600;
theta = 2004.3109*T_TT/3600 - 0.42665*T_TT^2/3600 - 0.041833*T_TT^3/3600;
z = 2306.2181*T_TT/3600 + 1.09468*T_TT^2/3600 + 0.018203*T_TT^3/3600;

% Check for angle errors against Vallado values.
disp("Check for errors in eta, theta, z angles ==========================")
etaErr = 0.0273055 - eta
thetaErr = 0.0237306 - theta
zErr = 0.0273059 - z


%% Convert MOD to GCRF.
% Convert degrees to radians.
eta = eta*pi/180;
theta = theta*pi/180;
z = z*pi/180;

% Rotate to the final ECI (GCRF) coordinates.
R323 = [ cos(theta)*cos(z)*cos(eta) - sin(z)*sin(eta),  sin(z)*cos(theta)*cos(eta) + sin(eta)*cos(z),  sin(theta)*cos(eta)
        -sin(eta)*cos(theta)*cos(z) - sin(z)*cos(eta), -sin(z)*sin(eta)*cos(theta) + cos(z)*cos(eta), -sin(theta)*sin(eta)
        -sin(theta)*cos(z),                            -sin(theta)*sin(z)                              cos(theta)];

rGCRFcalc = R323*rmod;
vGCRFcalc = R323*vmod;

% Check for errors in final calcualted GCRF values.
disp("Check for final errors ============================================")
rErr = (rGCRF - rGCRFcalc)'
vErr = (vGCRF - vGCRFcalc)'
