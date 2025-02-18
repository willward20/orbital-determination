function [a,e,inc,Om,w,v,Tp,P,M] = cart2kepv2(rVec, vVec, mu)
% Converts Cartesian coordinates to Keplerian orbital elements. 
% 
% Version 2 includes the time sine perigee and period.
%
% INPUTS
%
% rVec = (3x1) position [km]
% vVec = (3x1) velocity [km/s]
% mu = [km^3/s^2]
%
% OUTPUTS
%
% a = semi-major axis [km]
% e = eccentricity [1]
% inc = inclination [rad]
% Om = longitude of ascending node [rad]
% w = argument of periapsis [rad]
% v = true anomoly [rad]
% Tp = time of perigee passage [sec]
% P = period [sec]
% M = mean anomaly [rad]
%
% REFERENCES
%
% Cartesian State Vectors -> Keplerian Orbit Elements (Rene Schwarz)
% Statistical Orbit Determination by Tapley, Schutz, Born. 
% https://orbital-mechanics.space/time-since-periapsis-and-keplers-equation/elliptical-orbit-example.html
% +============================================================+
    r = norm(rVec);

    % Calculate orbital momentum vector and its norm.
    hVec = cross(rVec, vVec); % [km^2/s]
    h = norm(hVec); % [km^2/s]

    % Calculate the eccentricity.
    AVec = crossProdEquiv(vVec)*hVec - mu*rVec/r; % [km^3/s^2]
    eVec = AVec / mu; % [unitless]
    e = norm(eVec); % [unitless]

    % Calculate the semi-major axis.
    a = h^2 / (mu*(1 - e^2)); % [km]

    % Calculate the inclination.
    inc = acos(hVec(3) / h); % [rad]

    % Calculate the true anomaly. 
    if dot(rVec,vVec) >= 0
        v = acos(dot(eVec,rVec)/(e*r)); % [rad]
    else
        v = 2*pi - acos(dot(eVec,rVec)/(e*r)); % [rad]
    end

    % Calculate the vector pointing towards the ascending node.
    jVec = [0,0,1]';
    nVec = crossProdEquiv(jVec)*hVec;

    % Calculate the longitude of ascending node.
    n = norm(nVec);
    if nVec(2) >= 0
        Om = acos(nVec(1)/n); % [rad]
    else
        Om = 2*pi - acos(nVec(1)/n); % [rad]
    end

    % Calculate the argument of periapsis.
    if eVec(3) >= 0
        w = acos(dot(nVec, eVec) / (n*e)); % [rad]
    else
        w = 2*pi - acos(dot(nVec, eVec) / (n*e)); % [rad]
    end

    % Calculate the eccentricity anomaly (not returned). 
    E = 2*atan2(tan(v/2), sqrt((1+e)/(1-e))); % [rad]

    % Calculate the mean anomaly. 
    M = E - e*sin(E); % [rad]

    % Calculate the period. 
    P = 2*pi*sqrt(a^3/mu); % [sec]

    % Calculate the time of perigee passage.
    Tp = M*P/2/pi; % [sec]
end