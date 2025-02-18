%% Problem 2: Two-Body Kinematics + Drag
clear all
clc

% Set the initial conditions. 
rVec = [-2436.45, -2436.45, 6891.037]'; % [km]
vVec = [5.088611, -5.088611, 0.0]'; % [km/s]
X0 = [rVec; vVec];

% Set integrator check values for two body + J2 + drag.
rCheck = [-5751.50585441435, 4721.10759954747, 2046.09502951715];
vCheck = [-0.797610370476, -3.656553079577, 6.139595227675];

% Constants
J2 = 0.00108248; 
mu = 398600.4; % [km^3/s^2]
REarth = 6378.145; % [km]
secPerDay = 86400; % [s]
CD = 2.0;
A = 3.6 / 1e6; % [km^2]
m = 1350; % [kg]
rho0 = (4e-13)*(1e9); % [kg/km^3]
r0 = 7298.145; % [km]
H = 200.0; % [km]
thetaDot = 7.29211585530066e-5; % [rad/s]

% Get symbolic equations for dU/dx.
[dUdx, dUdy, dUdz] = getGradU();

% Substitute numerical values into symbolic expressions
dUdx_num = subs(dUdx, {'J2', 'muu', 'RE'}, {J2, mu, REarth});
dUdy_num = subs(dUdy, {'J2', 'muu', 'RE'}, {J2, mu, REarth});
dUdz_num = subs(dUdz, {'J2', 'muu', 'RE'}, {J2, mu, REarth});

% Convert symbolic expressions to numerical functions
dUdx_func = matlabFunction(dUdx_num, 'Vars', [sym('x'), sym('y'), sym('z')]);
dUdy_func = matlabFunction(dUdy_num, 'Vars', [sym('x'), sym('y'), sym('z')]);
dUdz_func = matlabFunction(dUdz_num, 'Vars', [sym('x'), sym('y'), sym('z')]);

% Initialize ode45.
myoptions = odeset('RelTol',3e-14,'AbsTol',1e-16);
times = 0:20:secPerDay; % [s]

% Propagate the kinematics forward in time using J2 + drag.
[T,Y] = ode45(@(t, X) orbitKinematicsJ2Drag(t, X, dUdx_func, dUdy_func, ...
    dUdz_func, rho0, r0, H, thetaDot, CD, A, m), times, X0, myoptions);

% Extract the position and velocity. 
R = Y(:,1:3);
V = Y(:,4:6);
% Get the euc-norms of each row of positions and velocities.
RNorms = vecnorm(R, 2, 2); 
VNorms = vecnorm(V, 2, 2); 

% Compare final value with the check values. 
R(end,:) - rCheck
V(end,:) - vCheck

%% Problem 2a: Compute the Specific Energy
% Calculate the kinetic, potential, and change in total energy.
Ek = VNorms.^2 ./ 2; % [km^2/s^2]
Ep = (mu./RNorms).*(1 - J2*(REarth./RNorms).^2.*(1.5*(R(:,3)./RNorms).^2 - 0.5)); % U
dEt = Ek - Ep;

% Plot the change in total energy over time. 
figure(1)
plot(T, dEt - dEt(1))
set(gca(), 'fontsize', 12)
title("Change in Energy over Time", fontsize=15)
xlabel("Time [s]", fontsize=14)
ylabel("Change in Energy [km^2/s^2]", fontsize=14)


%% Problem 2b: Plot orbital elements over one day. 
a_t = [];
e_t = [];
inc_t = [];
Om_t = [];
w_t = [];
Tp_t = [];
P_t = [];

% Get the orbital elements for 2B + J2 + drag
for ii = 1:length(R)
    [a,e,inc,Om,w,v,Tp,P] = cart2kepv2(R(ii,:)', V(ii,:)', mu);
    a_t = [a_t; a];
    e_t = [e_t; e];
    inc_t = [inc_t; inc];
    Om_t = [Om_t; Om];
    w_t = [w_t; w];
    Tp_t = [Tp_t; Tp];
    P_t = [P_t; P];
end

% Load the orbital elements for 2B + J2
kep2BJ2 = load('2BJ2elements.mat');

% Plot the results.
figure(2)
subplot(6,1,1)
sgtitle('Difference in 2B+J2 and 2B+J2+Drag Orbital Elements Over Time', fontweight='bold')
plot(T,a_t - kep2BJ2.a_t)
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('a_t [km]', fontsize=14, fontweight='bold')
grid on
subplot(6,1,2)
plot(T,e_t - kep2BJ2.e_t)
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('e', fontsize=14, fontweight='bold')
grid on
subplot(6,1,3)
plot(T,inc_t - kep2BJ2.inc_t)
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('inc [rad]', fontsize=14, fontweight='bold')
grid on
subplot(6,1,4)
plot(T,Om_t - kep2BJ2.Om_t)
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('Om [rad]', fontsize=14, fontweight='bold')
grid on
subplot(6,1,5)
plot(T,w_t - kep2BJ2.w_t)
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('w [rad]', fontsize=14, fontweight='bold')
grid on
subplot(6,1,6)
plot(T,Tp_t - kep2BJ2.Tp_t)
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('Time [s]', fontsize=14, fontweight='bold')
ylabel('Tp [s]', fontsize=14, fontweight='bold')
grid on


%% Extra Plots
% Plot the orbit in 3D.
% figure(4)
% hold on
% grid on
% plot3(R(:,1), R(:,2), R(:,3),'k')
% scatter3(0,0,0,'k','filled')
% scatter3(rVec(1), rVec(2), rVec(3), 'g', 'filled')
% scatter3(rCheck(1), rCheck(2), rCheck(3), 'r', 'filled')
% hold off

% Plot the orbital velocity in 3D.
% figure(5)
% hold on
% grid on
% plot3(V(:,1), V(:,2), V(:,3),'k')
% scatter3(0,0,0,'k','filled')
% scatter3(vVec(1), vVec(2), vVec(3), 'g', 'filled')
% scatter3(vCheck(1), vCheck(2), vCheck(3), 'r', 'filled')
% hold off