function dXdt = propSTM(t, X, AMat_func, dUdx_func, dUdy_func, ...
    dUdz_func, interp_rsun, interp_rmoon)
% Propagates the kinematics and STM of 2B+J2+Drag
% kinematic equations for a satellite. 
%
% INPUTS
%
% X = (56x1) state vector
%   X(1:3) = position in ECI [km]
%   X(4:6) = velocity in ECI [km]
%   X(7) = coefficient of drag C_D [unitless]
%   X(8:49) = vectorized STM (7 x 7)
%
%+============================================================+
    % Extract the inputs.
    R = X(1:3);
    V = X(4:6);
    C_D = X(7);
    STM = reshape(X(8:56), 7, 7);

    % Interpolate Sun and Moon positions instead of calling planetEphemeris
    rsun = interp_rsun(t);
    rmoon = interp_rmoon(t);

    % Calculate the A matrix.
    AMat = AMat_func(R(1), R(2), R(3), V(1), V(2), V(3), ...
        rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));

    % Compute the derivative for the states and STM. 
    % dStatesdt = AMat*[R;V;C_D];
    dR = V;
    dVx = dUdx_func(R(1), R(2), R(3), V(1), V(2), V(3), ...
        rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));
    dVy = dUdy_func(R(1), R(2), R(3), V(1), V(2), V(3), ...
        rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));
    dVz = dUdz_func(R(1), R(2), R(3), V(1), V(2), V(3), ...
        rmoon(1), rmoon(2), rmoon(3), rsun(1), rsun(2), rsun(3));
    dSTMdt = AMat*STM;

    % Reshape into the output derivative vector.
    dXdt = [dR; dVx; dVy; dVz; 0; reshape(dSTMdt, 49, 1)];
end