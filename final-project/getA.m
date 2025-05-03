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
    % Define symbolic constants.
    syms J2 J3 J4 J5 J6 J7 J8 J9 J10 muu RE C_D A m rho0 r0 H om_E real
    syms mu_Sun mu_Moon p_sr c_r A_panel real %delAT dUT1 pm real
    % Define symbolic variables.
    syms x y z vx vy vz real %UTC real
    syms xMoon yMoon zMoon xSun ySun zSun  real
    
    % Get the magnitude of the satellite position in ECI. 
    r = sqrt(x^2 + y^2 + z^2);

    % Get the acceleration.
    [dUdx, dUdy, dUdz] = getGradU();
    accel = [dUdx; dUdy; dUdz];

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