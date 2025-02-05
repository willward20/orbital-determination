function [T,Y] = propOrbit(rVec, vVec, mu, times)
% Integrates orbital kinematics forward for two orbits
% given position & velocity initial conditions. 
%
% INPUTS
%
% rVec = (3x1) position [km]
% vVec = (3x1) velocity [km/s]
% mu = double [km^3/s^2]
%
%+============================================================+
   % Initialize ode45.
    myoptions = odeset('RelTol',1e-12,'AbsTol',1e-20);
    X0 = [rVec; vVec]; % state vector (6x1)
    
    % Call ode45.
    [T,Y] = ode45(@orbitKinematics, times, X0, myoptions, mu); 
end