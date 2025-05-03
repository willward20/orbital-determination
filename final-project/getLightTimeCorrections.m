function [rho_new_biased, rhoDot_new] = getLightTimeCorrections( ...
    rho, R, V, A, gsR, gsV, b)
% 
% INPUTS
%
% rho = range measurement
% R = current position
% V = current velocity
% A = current acceleration
% gsR = ground station position
% gsV = ground station velocity
%
%+============================================================+
    % Set the initial light time estimate. 
    c = 299792.458;  % speed of light [km/s]
    lt = rho/c; % light time [s]

    % Set the tolerance for light time updates.
    tol = 1e-6; % [km] (converts to 1 mm)
    delta = 100; 
    R_old = R;

    iter = 0;

    while delta > tol
        % Propagate the state backwards in time using the constant
        % acceleration model. 
        dt = -lt;
        R_new = R + V*dt + 0.5*A*dt^2; 
        V_new = V + A*dt;

        % Back propagate the ground station position and velocity. 
        gsR_new = gsR + gsV*dt;
    
        % Simulate a new range measurment.
        rho_new = norm(R_new - gsR_new);
        rho_new_biased = norm(R_new - gsR_new) + b;
        rhoDot_new = (R_new - gsR_new).'*(V_new - gsV)/rho_new;
    
        % Get the new light time.
        lt = rho_new/c;
    
        % Compute the change.
        delta = norm(R_new - R_old);

        R_old = R_new;

        iter = iter + 1;
    end
end