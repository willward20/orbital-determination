% Test cart2kep vs ijk2keplerian.
clear all
clc

% Initial conditions.
rVec = [-2436.45, -2436.45, 6891.037]'; % position [km]
vVec = [5.088611, -5.088611, 0.0]'; % velocity [km/2]

[a,e,inc,Om,w,v] = cart2kep(rVec, vVec); % a [km], e [1], else [rad]

fprintf('    Semi-major axis (a): %d [km]\n', a)
fprintf('       Eccentricity (e): %d \n', e)
fprintf('        Inclination (i): %d [rad] \n', inc)
fprintf('              RAAN (Om): %d [rad] \n', Om)
fprintf('Argument of Perigee (w): %d [rad] \n', w)
fprintf('       True Anomaly (v): %d [rad] \n', v)