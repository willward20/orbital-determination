clear all
clc
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
R0 = [6984.45711518852;
      1612.2547582643;
      13.0925904314402]; % [km] 
V0 = [-1.67667852227336;
       7.26143715396544;
       0.259889857225218]; % [km/s]

% Tucker's Initial State in ECI
R0best = [6978642.93555757;
      1616591.37423624;
      19501.3039227613]./1000; % [km]
V0best = [-1663.03111545246;
       7260.77059418911;
       270.559997671147]./1000; % [km/s]

% Tracking Statation Coordinates (ECEF)
ts1 = [-6143584; 1364250; 1033743]/1000; % [km]
ts2 = [1907295; 6030810; -817119]/1000; % [km]
ts3 = [2390310; -5564341; 1994578]/1000; % [km]

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
mu_Moon = 4902.800066; % [km^3/s^2] Moon's gravitation
om_E = 7.292115146706979e-5; % [rad/s] Earth's rot vel

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

% Time Constants
sec_per_day = 24*60*60;

% Load the tracking data.
LEO_data_part1 = load('LEO_DATA_Apparent_3Days.mat');
LEO_data_part1 = cell2mat(struct2cell(LEO_data_part1));
LEO_data_part2 = load('LEO_DATA_Apparent_Days4-6.mat');
LEO_data_part2 = cell2mat(struct2cell(LEO_data_part2));
LEO_data_part2(:,2) =  3*sec_per_day + LEO_data_part2(:,2);
LEO_data = [LEO_data_part1; LEO_data_part2];

% Define the tracking data noise standard deviations.
ts1_sigmas = [(10/1000); 0.5*(1e-3)];
ts2_sigmas = [(5/1000); 1*(1e-3)];
ts3_sigmas = [(10/1000); 0.5*(1e-3)];


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


% REPLACE THESE  WITH THE RIGHT ONES
delAT = 37; % [seconds]
dUT1 = 0.1440613; %0.1964940; % [seconds]
% Polar motion values. 
% pm = [0.015447 0.288187]./3600; % [degrees]
pm = [0.020888 0.381012]./3600; % [degrees]


% Get symbolic function for A matrix.
AMat = getA();

% Substitute numerical values into A
AMat_num = subs(AMat, {'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8', 'J9', 'J10', 'muu', 'RE', 'C_D', 'A', 'm', 'rho0', 'r0', ...
    'H', 'om_E', 'mu_Sun', 'mu_Moon', 'p_sr', 'c_r', 'A_panel', 'A_Y', 'A_Z'}, ...
    {J2, J3, J4, J5, J6, J7, J8, J9, J10, mu, RE, C_D, AreaX, m, rho0, r0, H, om_E, mu_Sun, mu_Moon, p_sr, c_r, A_panel, AreaY, AreaZ});

% Convert symbolic expressions to numerical functions
AMat_func = matlabFunction(AMat_num, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz'), sym('xMoon'), ...
    sym('yMoon'), sym('zMoon'), sym('xSun'), sym('ySun'), sym('zSun')]);


% Get symbolic functions for HTilde for each station.
HtildeFunc = getHtilde();

% Substitute numerical values for tracking station 1. 
Htilde_ts1 = subs(HtildeFunc, {'C_D'}, {C_D});
% Substitute numerical values for tracking station 2. 
Htilde_ts2 = subs(HtildeFunc, {'C_D'}, {C_D});
% Substitute numerical values for tracking station 3. 
Htilde_ts3 = subs(HtildeFunc, {'C_D'}, {C_D});

% Convert symbolic expressions to numerical functions
Ht_ts1_func = matlabFunction(Htilde_ts1, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz'), ...
    sym('xI'), sym('yI'), sym('zI'), ...
    sym('vxI'), sym('vyI'), sym('vzI')]);
Ht_ts2_func = matlabFunction(Htilde_ts2, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz'), ...
    sym('xI'), sym('yI'), sym('zI'), ...
    sym('vxI'), sym('vyI'), sym('vzI')]);
Ht_ts3_func = matlabFunction(Htilde_ts3, 'Vars', [sym('x'), sym('y'), ...
    sym('z'), sym('vx'), sym('vy'), sym('vz'), ...
    sym('xI'), sym('yI'), sym('zI'), ...
    sym('vxI'), sym('vyI'), sym('vzI')]);


% Initialize ode45.
myoptions = odeset('RelTol',3e-14,'AbsTol',1e-14);


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


% Prepare for EKF
% Unpack the LEO Data
LEO_data = LEO_data(2134:end, :);
% LEO_data = LEO_data(LEO_data(:,1) == 1, :);
stationVec = LEO_data(:,1); % station ID number.
tVecM = LEO_data(:,2); % measurement times
zVec = LEO_data(:,3:4); % measurement data. 

% Measurement uncertainty.
R1 = diag(ts1_sigmas.^2);
R2 = diag(ts2_sigmas.^2);
R3 = diag(ts3_sigmas.^2);

% Process noise statistics (add later?)
vBar = zeros(6,1); % assume unbiased noise for now

% Define the process noise trans. mat (PNTM) for fixed 60 second intervals.
sigX = (1e-6)/5;
sigY = (1e-10)/0.1;
sigZ = (1e-6)/0.5;


% The prior state PDF (initial condition and covar)
xBar = [R0best; V0best];
PBar = diag([10 10 10 0.1 0.001 0.1]);

% Initialize the EKF
len = length(tVecM); % iterations
xkHat = xBar;
Pk = PBar;

% Mats for recording.
xHatMat = zeros(6, len+1); % 6 states (don't estimate C_D) 
PkMat = zeros(6,6,len+1); % covariance
delzks = zeros(2, len);
prefit_delzks = zeros(len, 2); %length(measurementTs),2);
postfit_delzks = zeros(len, 2); %length(measurementTs),2);
Pzzks = zeros(2,2,len); %length(measurementTs));

% Add the beginning time.
tVecM = [0; tVecM];
% Add the end time.
tVecM = [tVecM; 7*sec_per_day];

% EKF Algorithm (Lecture 16, slide 24)
for k = 1:len+1
    msg = "Iteration %d/%d";
    str = sprintf(msg,k,len);
    disp(str)

    % Time Update.
    if tVecM(k) == 0
        % Don't perform a time update because we're already at t = 0. 
        xkBar = xBar;
        PkBar = PBar;
        Fk = eye(6); % STM from t=0 to t=0. 
    else
        % Prepare ode45.
        Fk_padded = eye(7);
        xSTMk = [xkHat; C_D; reshape(Fk_padded, 7^2,1)];

        % Propagate the kinematics forward in time.
        if k == 1
            [T,Y] = ode45(@propStateAndSTM, [0:60:tVecM(k)]', xSTMk, myoptions, AMat_func, ...
                dUdx_func, dUdy_func, dUdz_func, rsunMat, rmoonMat); 
        else
            [T,Y] = ode45(@propStateAndSTM, tVecM(k-1:k), xSTMk, myoptions, AMat_func, ...
                dUdx_func, dUdy_func, dUdz_func, rsunMat, rmoonMat); 
        end
        
        % Extract the next state vector.
        xkBar = Y(end,1:6)';
        % Extract the STM.
        Fk_padded = reshape(Y(end,8:56)', 7, 7);
        Fk = Fk_padded(1:6, 1:6);
        % Update PkBar for each model.
        if k == 1
            PkBar = Fk*Pk*Fk';% + getQ(tVecM(k), sigX, sigY, sigZ);
        else
            PkBar = Fk*Pk*Fk' + getQ(tVecM(k) - tVecM(k-1), sigX, sigY, sigZ);
        end
    end

    if k == len+1 || k == 1
        % There is no measurement. 
        xkHat = xkBar;
        Pk = PkBar; 
    else
        % Perform a measurement update. 

        % Update the time.
        UTCplus = datevec(UTCdatetime + seconds(tVecM(k)));
    
        % Get the ground station location at the new time
        % and get the script Hk (Htilde) for that station. 
        % also get the R mat for that station. 
        if stationVec(k) == 1
            [gsR, gsV] = ecef2eci(UTCplus, ts1, [0;0;0], 'dAT', delAT, 'dUT1', dUT1, 'pm', pm);
            scriptHk = Ht_ts1_func(xkBar(1), xkBar(2), xkBar(3), ...
                xkBar(4), xkBar(5), xkBar(6), gsR(1), gsR(2), gsR(3), ...
                gsV(1), gsV(2), gsV(3));
            R = R1;
        elseif stationVec(k) == 2
            [gsR, gsV] = ecef2eci(UTCplus, ts2, [0;0;0], 'dAT', delAT, 'dUT1', dUT1, 'pm', pm);
            scriptHk = Ht_ts2_func(xkBar(1), xkBar(2), xkBar(3), ...
                xkBar(4), xkBar(5), xkBar(6), gsR(1), gsR(2), gsR(3), ...
                gsV(1), gsV(2), gsV(3));
            R = R2;
        else
            [gsR, gsV] = ecef2eci(UTCplus, ts3, [0;0;0], 'dAT', delAT, 'dUT1', dUT1, 'pm', pm);
            scriptHk = Ht_ts3_func(xkBar(1), xkBar(2), xkBar(3), ...
                xkBar(4), xkBar(5), xkBar(6), gsR(1), gsR(2), gsR(3), ...
                gsV(1), gsV(2), gsV(3));
            R = R3;
        end

        % Remove the C_D part of Htilde.
        scriptHk = scriptHk(:,1:6);
    
        % Simulate measurements. 
        range = norm(xkBar(1:3) - gsR);
        rangerate = (xkBar(1:3) - gsR).'*(xkBar(4:6) - gsV)/range;
    
        % Extract the satellite position and velocities.
        pos = xkBar(1:3);
        vel = xkBar(4:6);
    
        % Get the approximate sun and moon positions.
        rmoon = rmoonMat((tVecM(k)/60)+1, :);
        rsun = rsunMat((tVecM(k)/60)+1, :);
    
        % Get the acceleration of the satellite. 
        dVx = dUdx_func(pos(1), pos(2), pos(3), vel(1), vel(2), vel(3), ...
            rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));
        dVy = dUdy_func(pos(1), pos(2), pos(3), vel(1), vel(2), vel(3), ...
            rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));
        dVz = dUdz_func(pos(1), pos(2), pos(3), vel(1), vel(2), vel(3), ...
            rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));
        acc = [dVx; dVy; dVz];
    
        % Compute light time corrections
        if stationVec(k) == 3
            bias = 20/1000; % 20 meters
        else
            bias = 0.0; 
        end
        [range, rangerate] = getLightTimeCorrections( ...
            range, pos, vel, acc, gsR, gsV, bias);
        
        % Calculate measurement residuals.
        delzk = zVec(k,:)' - [range; rangerate]
    
        % Update estimate. 
        Kk = PkBar*scriptHk'*inv(scriptHk*PkBar*scriptHk' + R);
    
        xkHat = xkBar + Kk*delzk; % don't do measurement update for now
    
        I = eye(size(Kk*scriptHk));
        Pk = (I - Kk*scriptHk)*PkBar*(I - Kk*scriptHk)' + Kk*R*Kk';
    
        % Calc Residuals and Pzzk.
        Pzzk = scriptHk*PkBar*scriptHk' + R;
    
        % Simulate measurements again for post fit. 
        range = norm(xkHat(1:3) - gsR);
        rangerate = (xkHat(1:3) - gsR).'*(xkHat(4:6) - gsV)/range;
    
        % Calculate measurement residuals.
        postfit_delzk = zVec(k,:)' - [range; rangerate];
    
        % Record the measurements and measurement error. 
        prefit_delzks(k,:) = delzk';
        postfit_delzks(k,:) = postfit_delzk';
        Pzzks(:,:,k) = Pzzk;
    end

    % Record data.
    xHatMat(:,k) = xkHat;
    PkMat(:,:,k) = Pk;

    if isnan(eig(Pk))
        break
    end
end

Xt_final = xHatMat(1:6,end);
cov = PkMat(1:3,1:3,end);

%% Propagate to the end of the seventh day.
% next_times = tVecM(end):60:sec_per_day*7; % [s]
% XSTMinit = [xHatMat(:,end); C_D; reshape(eye(7), 7^2, 1)];
% 
% % Propagate the kinematics forward in time.
% [T2,Y2] = ode45(@propStateAndSTM, next_times, XSTMinit, myoptions, AMat_func, ...
%     dUdx_func, dUdy_func, dUdz_func, rsunMat, rmoonMat);
% 
% % Extract the components of Ystar
% Xt_final = Y2(:,1:7);
% STMt_final = reshape(Y2(:,8:56)', 7, 7, length(next_times));

% Compute relative error after 24 hours using truth. 
% rel_err = abs((Xt_final(end,1:6)' - truth)./truth)

%% Save the final data to a MAT file. 
% Get the final covariance.
% cov = STMt_final(1:6,1:6,end)*Pk*STMt_final(1:6,1:6,end)';
ward_pos_case = Xt_final(1:3)
ward_cov_case = cov(1:3,1:3,end)
% save('caseG.mat', "ward_pos_case", "ward_cov_case")

%% Plot the measurement residuals (post-fit)
% Plot MMSE Values
sigMat = zeros(len,2);
for k = 1:len
    sigMat(k,:) = sqrt(diag(Pzzks(:,:,k)))';
end
figure(1)
scatter(tVecM(1:len), postfit_delzks(:,1), 'filled')
hold on
grid on
plot(tVecM(1:len), postfit_delzks(:,1) - 3*sigMat(:,1), 'k', 'LineStyle','--','LineWidth',1.5)
plot(tVecM(1:len), postfit_delzks(:,1) + 3*sigMat(:,1), 'k', 'LineStyle','--','LineWidth',1.5)
title('Range Measurement Postfit Residuals over Time', 'FontSize', 16)
xlabel('Time [s]', 'FontSize', 14)
ylabel('Range Error [km]', 'FontSize', 14)
set(gca, 'FontSize', 12)

figure(2)
scatter(tVecM(1:len), postfit_delzks(:,2), 'filled')
hold on
grid on
plot(tVecM(1:len), postfit_delzks(:,2) - 3*sigMat(:,2), 'k', 'LineStyle','--','LineWidth',1.5)
plot(tVecM(1:len), postfit_delzks(:,2) + 3*sigMat(:,2), 'k', 'LineStyle','--','LineWidth',1.5)
title('Range-Rate Measurement Postfit Residuals over Time', 'FontSize', 16)
xlabel('Time [s]', 'FontSize', 14)
ylabel('Range-Rate Error [km/s]', 'FontSize', 14)
set(gca, 'FontSize', 12)



%% Plot the velocity over time.
% figure(3)
% plot(tVecM, xHatMat(1,:))
% hold on;
% grid on;
% plot(tVecM, xHatMat(2,:))
% plot(tVecM, xHatMat(3,:))
% legend('x','y','z')