function [dUdx, dUdy, dUdz] = getGradU()
% Returns symbolic equations for the gradient
% of U, which includes J2 (accounting for the
% Earth's oblateness), drag, and the sun/moon.
%
% INPUTS
% 
% None
%
% OUTPUTS
%
% dUdx = partial derivative of U w.r.t. x 
% dUdy = partial derivative of U w.r.t. y
% dUdz = partial derivative of U w.r.t. z 
%
% +============================================================+
    syms x y z vx vy vz J2 muu RE C_D A m rho0 r0 H om_E real
    syms xMoon yMoon zMoon xSun ySun zSun mu_Sun mu_Moon real

    r = sqrt(x^2 + y^2 + z^2);
    U = (muu/r)*(1 - J2*(RE/r)^2*(1.5*(z/r)^2 - 0.5));
    dUdx_2BJ2 = simplify(diff(U,x));
    dUdy_2BJ2 = simplify(diff(U,y));
    dUdz_2BJ2 = simplify(diff(U,z));

    % Acceleration due to drag.
    rhoA = rho0*exp(-(r - r0)/H);
    VAvec = [vx + om_E*y; ...
             vy - om_E*x; ...
             vz];
    VA = sqrt((vx + om_E*y)^2 + (vy - om_E*x)^2 + vz^2);
    accDrag = -0.5*C_D*(A/m)*rhoA*VA*VAvec;

    % Calculate the sun and moon positions relative to the satellite.
    RSat = [x; y; z];
    RSun = [xSun; ySun; zSun];
    RMoon = [xMoon; yMoon; zMoon];
    RSat2Sun = RSun - RSat;
    RSat2Moon = RMoon - RSat;
    
    % Calculate the acceleration due to sun and moon (Vallado Accelerations Due to Third Body (p. 515 book). 
    aSun = mu_Sun*((RSat2Sun./norm(RSat2Sun)^3) - RSun./norm(RSun)^3);
    aMoon = mu_Moon*((RSat2Moon./norm(RSat2Moon)^3) - RMoon./norm(RMoon)^3);

    % Total acceleration.
    accel = accDrag + [dUdx_2BJ2; dUdy_2BJ2; dUdz_2BJ2] + aSun + aMoon;

    % Return the components
    dUdx = accel(1);
    dUdy = accel(2);
    dUdz = accel(3);
end