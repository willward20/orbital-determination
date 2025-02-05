function dX = orbitKinematics(t, X, mu)
%  Defines kinematic equations for an orbit. 
%
% INPUTS
%
% rVec = (3x1) position [km]
% vVec = (3x1) velocity [km/s]
%
%+============================================================+
    % Unpack the state vector.
    R = X(1:3); % position
    V = X(4:6); % velocity
    % Built the derivative vector.
    dX = zeros(6,1);
    dX(1:3) = V; % Rdot
    dX(4:6) = -mu*R/norm(R)^3; % Vdot
end