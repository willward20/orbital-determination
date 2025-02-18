%% Problem 1: J2 Kinematics
clear all
clc

% Set the initial conditions. 
rVec = [-2436.45, -2436.45, 6891.037]'; % [km]
vVec = [5.088611, -5.088611, 0.0]'; % [km/s]
X0 = [rVec; vVec];

% Set integrator check values for two body + J2.
rCheck = [-5751.49900721589, 4721.14371040552, 2046.03583664311];
vCheck = [-0.797658631074, -3.656513108387, 6.139612016678];

% Constants
J2 = 0.00108248; 
mu = 398600.4; % [km^3/s^2]
REarth = 6378.145; % [km]
secPerDay = 86400; % [s]

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

% Propagate the kinematics forward in time using J2 kinematics.
[T,Y] = ode45(@(t, X) orbitKinematicsJ2(t, X, dUdx_func, dUdy_func, dUdz_func), ...
              times, X0, myoptions);

% Extract the position and velocity. 
R = Y(:,1:3);
V = Y(:,4:6);
% Get the euc-norms of each row of positions and velocities.
RNorms = vecnorm(R, 2, 2); 
VNorms = vecnorm(V, 2, 2); 

% Compare final value with the check values. 
% R(end,:) - rCheck
% V(end,:) - vCheck

%% Problem 1b: Plot orbital elements over one day. 
a_t = [];
e_t = [];
inc_t = [];
Om_t = [];
w_t = [];
Tp_t = [];
P_t = [];
M_t = [];

for ii = 1:length(R)
    % Get the orbital elements
    [a,e,inc,Om,w,v,Tp,P,M] = cart2kepv2(R(ii,:)', V(ii,:)', mu);
    a_t = [a_t; a];
    e_t = [e_t; e];
    inc_t = [inc_t; inc];
    Om_t = [Om_t; Om];
    w_t = [w_t; w];
    Tp_t = [Tp_t; Tp];
    P_t = [P_t; P];
    M_t = [M_t; M];
end

% Plot the results.
figure(1)
subplot(6,1,1)
sgtitle('Orbital Elements Over Time (Two-Body + J2 Kinematics)', fontweight='bold')
plot(T,a_t)
hold on
xline(P_t(1)*[1:12])
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('a_t [km]', fontsize=14, fontweight='bold')
grid on
subplot(6,1,2)
plot(T,e_t)
hold on
xline(P_t(1)*[1:12])
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('e', fontsize=14, fontweight='bold')
grid on
subplot(6,1,3)
plot(T,inc_t)
hold on
xline(P_t(1)*[1:12])
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('inc [rad]', fontsize=14, fontweight='bold')
grid on
subplot(6,1,4)
plot(T,Om_t)
hold on
xline(P_t(1)*[1:12])
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('Om [rad]', fontsize=14, fontweight='bold')
grid on
subplot(6,1,5)
plot(T,w_t)
hold on
xline(P_t(1)*[1:12])
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ylabel('w [rad]', fontsize=14, fontweight='bold')
grid on
subplot(6,1,6)
plot(T,Tp_t)
hold on
xline(P_t(1)*[1:12])
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('Time [s]', fontsize=14, fontweight='bold')
ylabel('Tp [s]', fontsize=14, fontweight='bold')
grid on

% Save orbital elements to a MAT file to compare with 
% results from Problem 2 (2B + J2 + drag).
save 2BJ2elements.mat a_t e_t inc_t Om_t w_t Tp_t;


%% Problem 1c: Compute the Specific Energy
% Calculate the kinetic, potential, and change in total energy.
Ek = VNorms.^2 ./ 2; % [km^2/s^2]
Ep = (mu./RNorms).*(1 - J2*(REarth./RNorms).^2.*(1.5*(R(:,3)./RNorms).^2 - 0.5)); % U
dEt = Ek - Ep;

% Plot the change in total energy over time. 
figure(2)
plot(T, dEt - dEt(1))
set(gca(), 'fontsize', 12)
title("Change in Energy over Time", fontsize=15)
xlabel("Time [s]", fontsize=14)
ylabel("Change in Energy [km^2/s^2]", fontsize=14)


%% Problem 1d: Compure the Angular Momentum
% Calculate the angular momentum vector.
h = cross(R, V, 2);

% Plot the angular momentum as a function of time (3D scatter plot). 
figure(3)
plot(T, h(:,3)-h(1,3))
set(gca(), 'fontsize', 12)
title("Angular Momentum over Time", fontsize=15)
xlabel("Time [s]", fontsize=14)
ylabel("h_z [km^2/s]", fontsize=14)

% Extra Plots
% % Plot the orbit in 3D.
% figure(4)
% hold on
% grid on
% plot3(R(:,1), R(:,2), R(:,3),'k')
% scatter3(0,0,0,'k','filled')
% scatter3(rVec(1), rVec(2), rVec(3), 'g', 'filled')
% scatter3(rCheck(1), rCheck(2), rCheck(3), 'r', 'filled')
% hold off
% 
% % Plot the orbital velocity in 3D.
% figure(5)
% hold on
% grid on
% plot3(V(:,1), V(:,2), V(:,3),'k')
% scatter3(0,0,0,'k','filled')
% scatter3(vVec(1), vVec(2), vVec(3), 'g', 'filled')
% scatter3(vCheck(1), vCheck(2), vCheck(3), 'r', 'filled')
% hold off