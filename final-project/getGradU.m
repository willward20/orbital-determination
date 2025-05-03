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
    syms x y z vx vy vz J2 J3 J4 J5 J6 J7 J8 J9 J10 muu RE C_D A m real
    syms rho0 r0 H om_E p_sr c_r A_panel A_Y A_Z real
    syms xMoon yMoon zMoon xSun ySun zSun mu_Sun mu_Moon real

    r = sqrt(x^2 + y^2 + z^2);
    % U = (muu/r)*(1 - J2*(RE/r)^2*(1.5*(z/r)^2 - 0.5));
    % dUdx_2BJ2 = simplify(diff(U,x));
    % dUdy_2BJ2 = simplify(diff(U,y));
    % dUdz_2BJ2 = simplify(diff(U,z));

    % 2B + J2 + J3 + J4
    U = (muu/r) * (1 ...
        - J2 * (RE/r)^2 * (1/2)*(3*(z/r)^2 - 1) ...
        - J3 * (RE/r)^3 * (1/2)*(5*(z/r)^3 - 3*(z/r)) ...
        - J4 * (RE/r)^4 * (1/8)*(35*(z/r)^4 - 30*(z/r)^2 + 3) ...
        - J5 * (RE/r)^5 * (1/8)*(63*(z/r)^5 - 70*(z/r)^3 + 15*(z/r)) ...
        - J6 * (RE/r)^6 * (1/16)*(231*(z/r)^6 - 315*(z/r)^4 + 105*(z/r)^2 - 5) ...
        - J7 * (RE/r)^7 * (1/16)*(429*(z/r)^7 - 693*(z/r)^5 + 315*(z/r)^3 - 35*(z/r)) ...
        - J8 * (RE/r)^8 * (1/128)*(6435*(z/r)^8 - 12012*(z/r)^6 + 6930*(z/r)^4 - 1260*(z/r)^2 + 35) ...
        - J9 * (RE/r)^9 * (1/128)*(12155*(z/r)^9 - 25740*(z/r)^7 + 18018*(z/r)^5 - 4620*(z/r)^3 + 315*(z/r)) ...
        - J10 * (RE/r)^10 * (1/256)*(46189*(z/r)^10 - 109395*(z/r)^8 + 90090*(z/r)^6 - 30030*(z/r)^4 + 3465*(z/r)^2 - 63));
    dUdx_2BJ2J3J4J5J6 = simplify(diff(U,x));
    dUdy_2BJ2J3J4J5J6 = simplify(diff(U,y));
    dUdz_2BJ2J3J4J5J6 = simplify(diff(U,z));

    % Acceleration due to drag on the forward face.
    rhoA = rho0*exp(-(r - r0)/H);
    VAvec = [vx + om_E*y; ...
             vy - om_E*x; ...
             vz];
    VA = sqrt(VAvec(1)^2 + VAvec(2)^2 + VAvec(3)^2);
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

    % Acceleration due to drag from solar panel. 
    theta = acos(dot(VAvec, RSat2Sun)/norm(VAvec)/norm(RSat2Sun));
    C_Dtheta = C_D * sin(theta)^2;
    A_eff = A_panel * abs(cos(theta));
    aDragPanel = -0.5*C_Dtheta*(A_eff/m)*rhoA*VA*VAvec;

    % Acceleration due to drag from the Z face.
    theta = acos(dot(VAvec, -RSat)/norm(VAvec)/norm(RSat));
    C_Dtheta = C_D * sin(theta)^2;
    A_eff = A_Z * abs(cos(theta));
    aDragZ = -0.5*C_Dtheta*(A_eff/m)*rhoA*VA*VAvec;

    % Assume that RSat, VAvec, and xHat are co-planar.
    % Then, approximate yHat.
    yVec = cross(-RSat, VAvec);
    theta = acos(dot(VAvec, yVec)/norm(VAvec)/norm(yVec));
    C_Dtheta = C_D * sin(theta)^2;
    A_eff = A_Y * abs(cos(theta));
    aDragY = -0.5*C_Dtheta*(A_eff/m)*rhoA*VA*VAvec;

    % Calculate the acceleration due to solar radiation. 
    aSRP = -p_sr * c_r * A_panel * RSat2Sun ./ norm(RSat2Sun) / m; 

    % Total acceleration.
    accel = [dUdx_2BJ2J3J4J5J6; dUdy_2BJ2J3J4J5J6; dUdz_2BJ2J3J4J5J6] + ...
        accDrag + aDragPanel + aDragZ + aSun + aMoon;% + aSRP;
    % accel = accDrag + aDragZ + aDragPanel; % + aDragY;
    % accel = aSun + aMoon;
    % accel = aSRP;

    % Return the components
    dUdx = accel(1);
    dUdy = accel(2);
    dUdz = accel(3);
end