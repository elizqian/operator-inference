% Full-order model solves the 1D Euler equations on periodic domain [0,2)
% using left-difference solver (supersonic flows only)

% INPUTS
% xgrid     mesh on which to solve
% dt        timestep
% init      struct with fields {rho, u, p} that defines initial condition
% Tfinal    final time for integration

% OUTPUTS
% s_all     conservative state data
% init      struct, adds field s0 which gives the initial conservative state
% time      vector of time values corresponding to s_all state data
% runtime   time to execute FOM

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 12 June 2019

function [s_all,init,time,runtime] = FOM(xgrid,dt,init,Tfinal)
start = tic;

N = length(xgrid);
dx      = xgrid(2)-xgrid(1);

% time discretization
time = 0:dt:Tfinal;
K = length(time)-1;

% parameters
gamma = 1.4;    % specific heat ratio

% set up conservative state initial condition
s1 = init.rho;
s2 = s1.*init.u;
s3 = init.p/(gamma-1) + 0.5.*s2.^2./s1;
init.s0 = [s1; s2; s3];

% allocate space for snapshots to be stored
s_all = zeros(3*N,K);
s_all(:,1) = [s1; s2; s3];

% pull primitive variables to start
rho = s1;
u = s2./rho;
p = (s3 - 0.5*rho.*u.^2)*(gamma-1);

for i = 1:K
    
    % calculate fluxes
    dF1dx = (rho.*u    - rho([N,1:N-1]).*u([N,1:N-1]))/dx;
    dF2dx = (rho.*u.^2 + p - rho([N,1:N-1]).*u([N,1:N-1]).^2 - p([N,1:N-1]))/dx;
    dF3dx = ((s3+p).*u - (s3([N,1:N-1])+ p([N,1:N-1])).*u([N,1:N-1]))/dx;
    
    % timestep
    s1 = s1 - dt*dF1dx;
    s2 = s2 - dt*dF2dx;
    s3 = s3 - dt*dF3dx;
    
    % pull primitive variables after timestep
    rho = s1;
    u = s2./rho;
    p = (s3 - 0.5*rho.*u.^2)*(gamma-1);
    
    s_all(:,i) = [s1; s2; s3];
end
runtime = toc(start);