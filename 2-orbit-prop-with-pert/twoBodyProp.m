% Integrates two-body orbital kinematics forward 
% for one day (86400 seconds) from initial 
% position & velocity. (Satellite orbiting Earth). 
clear all
clc

% Set the initial conditions. 
rVec = [-2436.45, -2436.45, 6891.037]'; % [km]
vVec = [5.088611, -5.088611, 0.0]'; % [km/s]
X0 = [rVec; vVec]; % state vector (6x1)

% Set integrator check values. 
rCheck = [-5971.19544867343, 3945.58315019255, 2864.53021742433];
vCheck = [0.049002818030, -4.185030861883, 5.848985672439];

% Constants
J2 = 0.00108248; 
mu = 398600.4; % [km^3/s^2]
REarth = 6378.145; % [km]
secPerDay = 86400; % [s]

% Initialize ode45.
myoptions = odeset('RelTol',3e-14,'AbsTol',1e-16);
times = 0:20:secPerDay; % [s]
    
% Propagate the kinematics forward in time.
[T,Y] = ode45(@orbitKinematics, times, X0, myoptions, mu); 

% Extract the position and velocity. 
R = Y(:,1:3);
V = Y(:,4:6);

% Compare final value with the check values. 
R(end,:) - rCheck
V(end,:) - vCheck

% Plot the orbit in 3D.
figure(4)
hold on
grid on
plot3(R(:,1), R(:,2), R(:,3),'k')
scatter3(0,0,0,'k','filled')
scatter3(rVec(1), rVec(2), rVec(3), 'g', 'filled')
scatter3(rCheck(1), rCheck(2), rCheck(3), 'r', 'filled')
hold off

% Plot the orbital velocity in 3D.
figure(5)
hold on
grid on
plot3(V(:,1), V(:,2), V(:,3),'k')
scatter3(0,0,0,'k','filled')
scatter3(vVec(1), vVec(2), vVec(3), 'g', 'filled')
scatter3(vCheck(1), vCheck(2), vCheck(3), 'r', 'filled')
hold off


