function dX = twoBody2DKinematics(t, X)
%  Defines two-dimenstional two-body
%  kinematic equations for an orbit. 
%
% INPUTS
%
% X = (4x1) state vector
%   X(1) = x position [km]
%   X(2) = y position [km]
%   X(3) = x velocity [km/s]
%   X(4) = y velocity [km/s]
%
%+============================================================+
    % Unpack the state vector. 
    r = sqrt(X(1)^2 + X(2)^2);

    % Built the derivative vector.
    dX = zeros(4,1);
    dX(1:2) = X(3:4); % [xDot; yDot]
    dX(3) = -X(1) / r^3; % xDotDot
    dX(4) = -X(2) / r^3; % yDotDot
end