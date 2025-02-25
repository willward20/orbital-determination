function dXdt = prop2D2BKinSTM(t, X)
% Propagates the kinematics and STM of 
% two-dimenstional two-body kinematic equations 
% for an orbit. 
%
% INPUTS
%
% X = (20x1) state vector
%   X(1) = x position [km]
%   X(2) = y position [km]
%   X(3) = x velocity [km/s]
%   X(4) = y velocity [km/s]
%   X(5:20) = vectorized STM
%
%+============================================================+
    % X = [x1, x2, x1', x2', F(t,t0)^T]^T
    r = sqrt(X(1)^2 + X(2)^2);
    A = [0  0  3*X(1)^2*r^-5 - r^-3  3*X(1)*X(2)*r^-5     ;
         0  0  3*X(1)*X(2)*r^-5      3*X(2)^2*r^-5 - r^-3 ;
         1  0              0                         0    ;
         0  1              0                         0   ];
    F = X(5:20);
    dFdt = A'*reshape(F, 4, 4);

    dXdt = [ X(3);
             X(4);
            -X(1)/r^3;
            -X(2)/r^3
            reshape(dFdt, 16, 1)];
end