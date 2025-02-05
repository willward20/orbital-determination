clear all;
close all;
clear clf;
clc

% Define constants.
A = 1.34;
kmratio = 1;
phi = pi/3.0;
times = linspace(0, 20, 101);

%% Problem 1: Evaluate the analytic solution to the harmonic oscillator.
% Calculate the harmonic oscillator solution.
x = harmOscillataorSolution(A, kmratio, phi, times);

% Plot the results. 
figure(1)
plot(times, x)
grid on;
title('Displacement over Time', fontsize=18)
xlabel('t [s]', fontsize=16)
ylabel('x(t) [m]', fontsize=16)
set(gca(), 'fontsize', 12)

%% Problem 2: Solve the Diff Eq using ode45.
% Initialize ode45.
myoptions = odeset('RelTol',1e-12,'AbsTol',1e-20);
y0 = zeros(2,1);
y0(1) = A*cos(phi);
y0(2) = -A*sqrt(kmratio)*sin(phi);

% Call ode45.
[T,Y] = ode45(@harmOscillatorDiffEq, times, y0, myoptions, kmratio);

% Calculate state error. 
error = Y(:,1)' - x;

% Plot state errors. 
figure(2)
plot(times, error)
grid on;
title('Displacement Errors over Time', fontsize=18)
xlabel('t [s]', fontsize=16)
ylabel('Displacement Errors [m]', fontsize=16)
set(gca(), 'fontsize', 12)

%% Function Definitions. 
function traj = harmOscillataorSolution(A, kmratio, phi, times)
    % Return dispalcement over time. 
    traj = A*cos(sqrt(kmratio)*times + phi);
end

function dx = harmOscillatorDiffEq(t, x, kmratio)
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -kmratio*x(1);
end