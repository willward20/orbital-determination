function dX = orbitKinematicsJ2Drag(t, X, dUdx_func, dUdy_func, dUdz_func, ...
    rho0, r0, H, thetaDot, CD, A, m)
% Defines kinematic equations for an orbit, including
% the Earth's oblateness (J2). 
%
% INPUTS
%
%+============================================================+
    % Unpack the state vector.
    R = X(1:3); % position
    V = X(4:6); % velocity
    RNorm = norm(R, 2); 

    % Built the derivative vector using two-body + J2 kinematics.
    dX = zeros(6,1);
    dX(1:3) = V; % Rdot
    dX(4) = dUdx_func(R(1), R(2), R(3));
    dX(5) = dUdy_func(R(1), R(2), R(3));
    dX(6) = dUdz_func(R(1), R(2), R(3));

    % Calculate acceleration due to drag.
    rhoA = rho0*exp(-(RNorm - r0)/H);
    VAvec = [V(1) + thetaDot*R(2); ...
             V(2) - thetaDot*R(1); ...
             V(3)];
    VA = sqrt((V(1) + thetaDot*R(2))^2 + (V(2) - thetaDot*R(1))^2 + V(3)^2);
    accDrag = -0.5*CD*(A/m)*rhoA*VA*VAvec;

    % Add drag to derivative vector.
    dX(4) = dX(4) + accDrag(1);
    dX(5) = dX(5) + accDrag(2);
    dX(6) = dX(6) + accDrag(3);
end