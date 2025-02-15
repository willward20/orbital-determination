function dX = orbitKinematicsJ2(t, X, dUdx_func, dUdy_func, dUdz_func)
% Defines kinematic equations for an orbit, including
% the Earth's oblateness (J2). 
%
% INPUTS
%
%+============================================================+
    % Unpack the state vector.
    R = X(1:3); % position
    V = X(4:6); % velocity

    % Built the derivative vector.
    dX = zeros(6,1);
    dX(1:3) = V; % Rdot
    dX(4) = dUdx_func(R(1), R(2), R(3));
    dX(5) = dUdy_func(R(1), R(2), R(3));
    dX(6) = dUdz_func(R(1), R(2), R(3));
end