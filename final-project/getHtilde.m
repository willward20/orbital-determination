function Htilde = getHtilde()
% Returns an Htilde matrix for one ground station. 
%
% INPUTS
% 
% None
%
% OUTPUTS
%
% +============================================================+
    % Symbolic variables for the position and velocities
    % of the satellite and the ground station (I for instrument).
    syms x y z vx vy vz C_D xI yI zI vxI vyI vzI
    R = [x; y; z]; % satellite position
    V = [vx; vy; vz]; % satellite velocity
    RI = [xI; yI; zI]; % ground station position
    VI = [vxI; vyI; vzI]; % ground station velocity

    range = sqrt((R - RI).'*(R - RI));
    rangerate = (R - RI).'*(V - VI)/range;

    % Get the components of Htilde for range measurements. 
    dpdx = simplify(diff(range,x));
    dpdy = simplify(diff(range,y));
    dpdz = simplify(diff(range,z));
    dpdvx = simplify(diff(range,vx));
    dpdvy = simplify(diff(range,vy));
    dpdvz = simplify(diff(range,vz));
    dpdCD = simplify(diff(range,C_D));

    % Get the components of Htilde for range-rate measurements. 
    dpDotdx = simplify(diff(rangerate,x));
    dpDotdy = simplify(diff(rangerate,y));
    dpDotdz = simplify(diff(rangerate,z));
    dpDotdvx = simplify(diff(rangerate,vx));
    dpDotdvy = simplify(diff(rangerate,vy));
    dpDotdvz = simplify(diff(rangerate,vz));
    dpDotdCD = simplify(diff(rangerate,C_D));

    % Build Htilde.
    Htilde = [dpdx dpdy dpdz dpdvx dpdvy dpdvz dpdCD;
              dpDotdx dpDotdy dpDotdz dpDotdvx dpDotdvy dpDotdvz dpDotdCD];
end