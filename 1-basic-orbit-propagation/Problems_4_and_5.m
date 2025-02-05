% Integrates orbital kinematics forward for two orbits,
% given position & velocity initial conditions. 
clear all
clc

% Set the initial conditions. 
rVec = [-2436.45, -2436.45, 6891.037]'; % [km]
vVec = [5.088611, -5.088611, 0.0]'; % [km/s]
X0 = [rVec; vVec]; % state vector (6x1)
mu = 398600.5;  % [km^3/s^2]
a = 7.712184983762814e+03; % semi-major axis [km]
myoptions = odeset('RelTol',1e-12,'AbsTol',1e-20);

% Calculate the period and set time horizon.
P = sqrt(4*pi^2*a^3/mu); % [s]
times = 0:20:2*P; % [s]
    
% Propagate the kinematics forward in time.
[T,Y] = ode45(@orbitKinematics, times, X0, myoptions, mu); 

% Extract the position and velocity. 
R = Y(:,1:3);
V = Y(:,4:6);
% Get the euc-norms of each row of positions and velocities.
RNorms = vecnorm(R, 2, 2); 
VNorms = vecnorm(V, 2, 2); 
% Calculate acceleration over time and its norm. 
Acc = -mu*R./RNorms.^3; % divide each row of R by an element of RNorms.
AccNorms = vecnorm(Acc, 2, 2); 

% Calculate the angular momentum vector.
h = cross(R, V, 2);

% Calculate the kinetic, potential, and change in total energy.
Ek = VNorms.^2 ./ 2; % [km^2/s^2]
Ep = mu./RNorms;
dEt = Ek - Ep;

% Plot magnitudes of position, velocity, and acceleration over time. 
figure(1)
subplot(3,1,1)
sgtitle('Position, Velocity, and Acceleration Magnitudes Over Time', fontweight='bold')
plot(T,RNorms)
ax = gca;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
ylabel('P [km]', fontsize=14)
grid on
subplot(3,1,2)
plot(T,VNorms)
ax = gca;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
ylabel('V [km/s]', fontsize=14)
grid on
subplot(3,1,3)
plot(T,AccNorms)
ax = gca;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
xlabel('Time [s]', fontsize=14)
ylabel('A [km/s^2]', fontsize=14)
grid on

% Plot the angular momentum as a function of time (3D scatter plot). 
figure(2)
scatter3(h(:,1), h(:,2), h(:,3), 4, T, 'filled')
title("Angular Momentum over Time", fontsize=18)
xlabel("h_x [km^2/s]", fontsize=14, fontweight='bold')
ylabel("h_y [km^2/s]", fontsize=14, fontweight='bold')
zlabel("h_z [km^2/s]", fontsize=14, fontweight='bold')
set(gca(), 'fontsize', 12)

% Colormap and colorbar
colormap('winter'); 
c = colorbar('southoutside'); % Show color legend
c.Label.String = 'Time [s]';
c.Label.FontSize = 14;
clim([min(T) max(T)]);

% Plot the change in total energy over time. 
figure(3)
plot(T, dEt - dEt(1))
set(gca(), 'fontsize', 10)
title("Change in Energy over Time", fontsize=15)
xlabel("Time [s]", fontsize=14)
ylabel("Change in Energy [km^2/s^2]", fontsize=14)

% Plot the orbit in 3D.
% figure(4)
% hold on
% grid on
% plot3(R(:,1), R(:,2), R(:,3))
% scatter3(0,0,0,'k','filled')
% hold off


