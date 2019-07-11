% converts conservative state vector/matrix in which Euler equations are
% nonlinear to specific volume variables in which eqs are quadratic

% INPUT
% nl_state  conservative state data with [rho; rho*u; E] stacked in vectors

% OUTPUT
% quad_state    specific volume state data, [u; p; 1./rho] stacked

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 12 June 2019

function quad_state = nonlin2quad(nl_state)

gamma = 1.4;

ntemp	= size(nl_state,1);
N    	= ntemp/3;
rho     = nl_state(1:N,:);
E       = nl_state(2*N+1:end,:);
u       = nl_state(N+1:2*N,:)./rho;
p       = (gamma-1)*(E-0.5*rho.*u.^2);

quad_state = [u; p; 1./rho];