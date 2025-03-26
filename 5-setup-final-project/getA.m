function A = getA()
% Returns 
%
% INPUTS
% 
% None
%
% OUTPUTS
%
% +============================================================+
    syms x y z vx vy vz J2 muu RE C_D A m rho0 r0 H om_E real
    syms xMoon yMoon zMoon xSun ySun zSun mu_Sun mu_Moon real
    r = sqrt(x^2 + y^2 + z^2);

    % Potential function for 2B and J2.
    U = (muu/r)*(1 - J2*(RE/r)^2*(1.5*(z/r)^2 - 0.5));

    % Get components of acceleration due to 2B and J2. 
    dUdx = simplify(diff(U,x));
    dUdy = simplify(diff(U,y));
    dUdz = simplify(diff(U,z));

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
    accel = accDrag + [dUdx; dUdy; dUdz] + aSun + aMoon;

    % Calculate the partial derivative components of A. 
    dvdr = 0 * eye(3);
    dvdv = eye(3);
    dvdCD = [0; 0; 0];
    dadr = [diff(accel(1),x)   diff(accel(1),y)   diff(accel(1), z);
            diff(accel(2),x)   diff(accel(2),y)   diff(accel(2), z);
            diff(accel(3),x)   diff(accel(3),y)   diff(accel(3), z)];
    dadv = [diff(accel(1),vx)   diff(accel(1),vy)   diff(accel(1), vz);
            diff(accel(2),vx)   diff(accel(2),vy)   diff(accel(2), vz);
            diff(accel(3),vx)   diff(accel(3),vy)   diff(accel(3), vz)];
    dadCD = [diff(accel(1),C_D);
            diff(accel(2),C_D);
            diff(accel(3),C_D)];
    dCDdr = [0 0 0];
    dCDdv = [0 0 0];
    dCDdCD = 0; 

    A = [dvdr   dvdv   dvdCD;
         dadr   dadv   dadCD;
         dCDdr  dCDdv  dCDdCD];    
end