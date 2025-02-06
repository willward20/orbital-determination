function [rVec, vVec] = kep2cart(a,e,inc,Om,w,v)
% Keplerian orbital elements to cartesian coordinates.  
%
% INPUTS
%
% a = semi-major axis [km]
% e = eccentricity [1]
% inc = inclination [rad]
% Om = longitude of ascending node [rad]
% w = argument of periapsis [rad]
% v = true anomoly [rad]
%
% OUTPUTS
%
% rVec = (3x1) position [km]
% vVec = (3x1) velocity [km/s]
%
% REFERENCES
%
% Code sourced from: 
% https://orbital-mechanics.space/classical-orbital-elements/orbital-elements-and-the-state-vector.html#orbital-elements-state-vector
% +============================================================+
    mu = 398600.5; % [km^3/s^2]

    % Find angular momentum.
    h = sqrt(mu*a*(1-e^2));

    % Transform to perifocal frame.
    r_w = h^2 / mu / (1 + e * cos(v)) .* [cos(v) sin(v) 0];
    v_w = mu / h .* [-sin(v) e + cos(v) 0];

    % Rotate to perifocal frame.
    R1 = [cos(-w) -sin(-w) 0; sin(-w) cos(-w) 0; 0 0 1];
    R2 = [1 0 0; 0 cos(-inc) -sin(-inc); 0 sin(-inc) cos(-inc)];
    R3 = [cos(-Om) -sin(-Om) 0; sin(-Om) cos(-Om) 0; 0 0 1];
    rVec = r_w * R1 * R2 * R3;
    vVec = v_w * R1 * R2 * R3;
end