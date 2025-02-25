%% Problem 1: 
clear all
clc

% Set the initial conditions. 
X0 = [1; 0; 0; 1]; % [x0, y0, xDot0, yDot0]

% Set constants.
mu = 1;

% Initialize ode45.
myoptions = odeset('RelTol',3e-14,'AbsTol',1e-16);
times = 0:1:10; % [s]

%% Part (a): Integrate the Equations of Motion.
% Set checker values.
X10chk = [-0.839071529; -0.544021111; 0.544021111; -0.839071529];

% Propagate the kinematics forward in time.
[T,X] = ode45(@twoBody2DKinematics, times, X0, myoptions); 

% Check that integrated values pass the first check. 
disp('Part (a): Integration Errors: ')
err = X(end,:)' - X10chk

%% Part (b): Integrate the Perturbed Equations.
% Define the perturbation.
delX0 = [1e-6; -1e-6; 1e-6; 1e-6];
X0star = X0 - delX0;

% Define the full state vector (including the STM).
XSTM0star = [X0star; reshape(eye(4), 4^2, 1)];

% Set the checker values.
X10starchk = [-0.839031098; -0.544071486; 0.544076120; -0.839041244];
STMchk = [-19.2963174705 -1.0005919528 -1.5446240948 -20.5922746780
               24.5395368984  2.5430400375  3.3820224390  24.9959638293
              -26.6284485803 -1.2470410802 -2.0860289935 -27.5413748340
              -15.0754226454 -1.4570972848 -2.0011442064 -14.6674122500];

% Propagate the kinematics forward in time.
[Tstar,Ystar] = ode45(@prop2D2BKinSTM, times, XSTM0star, myoptions); 

Xstar = Ystar(:,1:4);
STM = reshape(Ystar(:,5:20)', 4, 4, 11);

% Check that integrated values pass the first check. 
disp('Part (b): Perturbed State Integration Errors: ')
err = Xstar(end,:)' - X10starchk
disp('Part (b): Perturbed STM Integration Errors: ')
err = STM(:,:,end) - STMchk

%% Part (c):


%% Part (d): Calculating the Perturbation Vector
% Set the checker values.
delX10_m1chk = [-0.000040431037; 0.000050375590; -0.000055009526; -0.000030284890];
delX10_m2chk = [-0.000040432624; 0.000050374483; -0.000055008811; -0.000030286882];

% Calculate delX using two methods.
delX_m1 = X - Xstar;
delX_m2 = pagemtimes(STM, delX0);

% Check the methods.
disp('Part (d): Perturbation Vector Error (m1): ')
err = delX_m1(end,:)' - delX10_m1chk
disp('Part (d): Perturbation Vector Error (m2): ')
err = delX_m2(:,:,end) - delX10_m2chk
disp('Part (d): Perturbation Vector Error (m1 - m2): ')
err = delX_m1(end,:)' - delX_m2(:,:,end)