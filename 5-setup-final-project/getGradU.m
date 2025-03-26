function [dUdx, dUdy, dUdz] = getGradU()
% Returns symbolic equations for the gradient
% of U, which includes J2 (accounting for the
% Earth's oblateness). 
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
    syms x y z J2 muu RE
    r = sqrt(x^2 + y^2 + z^2);
    U = (muu/r)*(1 - J2*(RE/r)^2*(1.5*(z/r)^2 - 0.5));
    dUdx = simplify(diff(U,x));
    dUdy = simplify(diff(U,y));
    dUdz = simplify(diff(U,z));
end