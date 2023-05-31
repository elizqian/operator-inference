function [dXdt,ind] = ddt(X,dt,scheme)
% Uses specified difference scheme to approximate dX/dt with uniform time 
% spacing of size dt
%
% INPUTS
% X         N-by-K data matrix where each column is the state at one time
% dt        time step
% scheme    specifies which scheme is used to approximate time derivative
% 
% OUTPUTS
% ind       M-by-1 vector of indices of state data in X that correspond to
%           calculated time derivative
% dXdt      N-by-M data matrix of state time derivatives corresponding to
%           the state data in X(:,ind)
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 11 July 2019

[~,K] = size(X);

switch scheme
    
    case {1,'FE'}   % forward differencing corresponding to forward Euler integration 
        dXdt = (X(:,2:end) - X(:,1:end-1))/dt;
        ind = 1:K-1;
        
    case {'BE'}     % backward differencing corresponding to backward Euler integration
        dXdt = (X(:,2:end) - X(:,1:end-1))/dt;
        ind = 2:K;
        
    case {2, '2c'}      % 2nd order central differencing
        dXdt = (X(:,3:end) - X(:,1:end-2))/(2*dt);
        ind = 2:K-1;
    
    case '2imp'         % 2nd order backward differencing (implicit)
        dXdt = (3*X(:,3:end) - 4*X(:,2:end-1) + X(:,1:end-2))/(2*dt);
        ind = 3:K;
        
    case '2ex'          % 2nd order forward differencing (explicit)
        dXdt = (-3*X(:,1:end-2) + 4*X(:,2:end-1) - X(:,3:end))/(2*dt);
        ind = 1:K-2;
    
    case {4,'4c'}       % 4th order central differencing
        % 4th order 5-point stencil from https://en.wikipedia.org/wiki/Five-point_stencil
        dXdt = (1/12*X(:,1:end-4) - 2/3*X(:,2:end-3) + 2/3*X(:,4:end-1) - 1/12*X(:,5:end))/dt;
        ind = 3:K-2;
        
    case '4imp'         % 4th order backward differencing
        dXdt = (25/12*X(:,5:end) - 4*X(:,4:end-1) + 3*X(:,3:end-2) - 4/3*X(:,2:end-3) + 1/4*X(:,1:end-4))/dt;
        ind = 5:K;
        
    case '4ex'          % 4th order forward differencing
        dXdt = (-25/12*X(:,1:end-4) + 4*X(:,2:end-3) - 3*X(:,3:end-2) + 4/3*X(:,4:end-1) - 1/4*X(:,5:end))/dt;
        ind = 1:K-4;
        
    otherwise
        error('Specified difference scheme not implemented in ddt.m')
end