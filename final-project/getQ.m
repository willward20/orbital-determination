function Q = getQ(del_t, sigX, sigY, sigZ)
% Returns a process noise transition matrix.
%
% INPUTS
% 
% del_t: change in time
% sigX, sigY, sigZ: covariances
%
% OUTPUTS
%
% Q: PNTM
%
% +============================================================+
    Q = [(0.5*del_t*sigX)^2   0   0   0.5*del_t*sigX^2   0   0;
          0   (0.5*del_t*sigY)^2   0   0   0.5*del_t*sigY^2   0;
          0   0   (0.5*del_t*sigZ)^2   0   0   0.5*del_t*sigZ^2;
          0.5*del_t*sigX^2   0   0   sigX^2   0   0;
          0   0.5*del_t*sigY^2   0   0   sigY^2   0;
          0   0   0.5*del_t*sigZ^2   0   0   sigZ^2];
end