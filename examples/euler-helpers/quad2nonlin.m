% converts specific volume state vector/matrix in which Euler equations are
% quadratic to conservative variables in which eqs are nonlinear

% INPUT
% quad_state    specific volume state data, [u; p; 1./rho] stacked

% OUTPUT
% nl_state  conservative state data with [rho; rho*u; E] stacked in vectors

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 12 June 2019

function nl_state = quad2nonlin(quad_state)

gamma = 1.4;

ntemp	= size(quad_state,1);
N    	= ntemp/3;
u = quad_state(1:N,:);
p = quad_state(N+1:2*N,:);
rho = 1./quad_state(2*N+1:3*N,:);
E = p/(gamma-1) + 0.5*rho.*u.^2;

nl_state = [rho; rho.*u; E];